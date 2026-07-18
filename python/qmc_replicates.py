#!/usr/bin/env python3
"""qmc_replicates.py — replicate-scatter error estimate for MoCafe runs.

Randomized quasi-Monte Carlo (`par%launch_sequence = 'sobol'`) is a
*randomized* scheme: one scrambled net gives no ordinary error bar, but
independent scrambles (different `par%qmc_seed`) are independent unbiased
replicates, and their scatter *is* the error estimate.  This tool runs a
namelist several times with different seeds, reads the outputs, and reports
the replicate mean and RMS for the integrated quantities.

Run it for the pseudo-random mode too (varying `par%iseed` instead) and it
reports the **effective photon-number gain**: the factor by which an ordinary
Monte-Carlo run would have to increase `no_photons` to reach the RMS the
quasi-random run achieved, estimated from the variance ratio at fixed cost
(RMS scales as N^-1/2, so gain = (rms_random / rms_sobol)^2).

Examples
--------
    # four scrambles of a sobol run
    python3 qmc_replicates.py sed_point.in --mode sobol --n 4

    # compare against four pseudo-random replicates and print the gain
    python3 qmc_replicates.py sed_point.in --compare --n 4 --np 8

The namelist is copied for each replicate with the seed and the output name
rewritten, so the original input file is never modified.
"""
from __future__ import annotations

import argparse
import os
import re
import glob
import shutil
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from mocafe_io import MoCafeFile, find_mocafe_outputs  # noqa: F401
    HAVE_IO = True
except Exception:                                          # pragma: no cover
    HAVE_IO = False
import h5py


# --------------------------------------------------------------------------
def set_par(text: str, key: str, value: str) -> str:
    """Set (or insert) par%<key> in a namelist body."""
    pat = re.compile(rf"^\s*par%{re.escape(key)}\s*=.*$", re.IGNORECASE | re.MULTILINE)
    line = f" par%{key} = {value}"
    if pat.search(text):
        return pat.sub(line, text, count=1)
    # insert just after the &parameters line
    return re.sub(r"(&parameters\s*\n)", r"\1" + line + "\n", text, count=1)


def integrated(path: str) -> dict:
    """Integrated quantities of one output file."""
    out = {}
    with h5py.File(path, "r") as f:
        for name in ("Scattered", "Direct", "Direct0"):
            key = f"{name}/data"
            if key in f:
                out[name] = float(np.asarray(f[key]).sum())
        for grp in ("Direct", "Scattered"):
            if grp in f and "TOT_LUM" in f[grp].attrs:
                out["TOT_LUM"] = float(np.atleast_1d(f[grp].attrs["TOT_LUM"])[0])
                break
    if "Scattered" in out and "Direct" in out:
        out["Total"] = out["Scattered"] + out["Direct"]
    return out


def run_one(infile: str, workdir: str, tag: str, mode: str, seed: int,
            nproc: int, exe: str, rundir: str) -> dict:
    """Run one replicate and return its integrated quantities.

    The run happens in `rundir` so that relative paths inside the namelist
    (`kext_file`, `source_spectrum`, ...) resolve exactly as they do for a
    normal run of that input.
    """
    text = open(infile).read()
    #--- MoCafe writes the output under its base name in the working
    #--- directory, so name it there and collect it afterwards.
    out_h5 = os.path.join(rundir, f"{tag}.h5")
    text = set_par(text, "out_file", f"'{tag}.h5'")
    text = set_par(text, "file_format", "'hdf5'")
    if mode == "sobol":
        text = set_par(text, "launch_sequence", "'sobol'")
        text = set_par(text, "qmc_seed", str(seed))
    else:
        text = set_par(text, "launch_sequence", "'random'")
        text = set_par(text, "iseed", str(seed))

    nml = os.path.join(workdir, f"{tag}.in")
    with open(nml, "w") as fh:
        fh.write(text)

    cmd = ["mpirun", "-np", str(nproc), exe, nml]
    res = subprocess.run(cmd, cwd=rundir, capture_output=True, text=True)
    if res.returncode != 0:
        sys.stderr.write(res.stdout[-2000:] + res.stderr[-2000:])
        raise RuntimeError(f"run failed: {' '.join(cmd)}")
    if not os.path.exists(out_h5):
        #--- MoCafe may decorate the name (an "_obs" stem, a per-observer
        #--- "_NNN" suffix); take the observer file it actually wrote.
        cand = [p for p in sorted(glob.glob(os.path.join(rundir, tag + "*.h5")))
                if not p.endswith(("_tau.h5", "_jlam.h5", "_dustsed.h5",
                                   "_bwdust.h5", "_stokes.h5"))]
        if not cand:
            raise RuntimeError(f"no output produced for {tag} in {rundir}")
        out_h5 = cand[0]
    stats = integrated(out_h5)
    for leftover in glob.glob(os.path.join(rundir, tag + "*")):
        if leftover.endswith((".h5", ".fits", ".fits.gz")):
            shutil.move(leftover, os.path.join(workdir, os.path.basename(leftover)))
    return stats


def summarize(results: list[dict]) -> dict:
    """Mean and RMS (unbiased, over replicates) of each quantity."""
    keys = sorted(set().union(*(r.keys() for r in results)))
    stats = {}
    for k in keys:
        v = np.array([r[k] for r in results if k in r], dtype=float)
        if v.size < 2:
            continue
        stats[k] = (v.mean(), v.std(ddof=1))
    return stats


# --------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="MoCafe namelist (.in)")
    ap.add_argument("--mode", choices=("sobol", "random"), default="sobol",
                    help="launch mode of the replicate set (default: sobol)")
    ap.add_argument("--compare", action="store_true",
                    help="also run pseudo-random replicates and report the gain")
    ap.add_argument("--n", type=int, default=4, help="replicates (default 4; 8 is better)")
    ap.add_argument("--np", dest="nproc", type=int, default=4, help="MPI ranks")
    ap.add_argument("--seed0", type=int, default=1000, help="first seed")
    ap.add_argument("--exe", default=None, help="MoCafe.x (default: ../MoCafe.x)")
    ap.add_argument("--cwd", default=None,
                    help="working directory for the runs, so relative paths in the "
                         "namelist resolve (default: the input file's directory)")
    ap.add_argument("--keep", action="store_true", help="keep the replicate outputs")
    args = ap.parse_args()

    exe = args.exe or os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "MoCafe.x"))
    if not os.path.exists(exe):
        sys.exit(f"MoCafe.x not found: {exe} (use --exe)")

    rundir = os.path.abspath(args.cwd or os.path.dirname(os.path.abspath(args.input)))
    workdir = tempfile.mkdtemp(prefix="qmcrep_")
    try:
        modes = [args.mode] + (["random"] if args.compare and args.mode == "sobol" else [])
        allstats = {}
        for mode in modes:
            print(f"--- {mode}: {args.n} replicates, -np {args.nproc} ---")
            res = []
            for i in range(args.n):
                seed = args.seed0 + i
                r = run_one(args.input, workdir, f"{mode}_{i}", mode, seed,
                            args.nproc, exe, rundir)
                print(f"  seed {seed}: " +
                      "  ".join(f"{k}={v:.6g}" for k, v in sorted(r.items())))
                res.append(r)
            allstats[mode] = summarize(res)

        for mode, st in allstats.items():
            print(f"\n{mode}: replicate mean +/- RMS ({args.n} replicates)")
            for k, (m, s) in st.items():
                rel = s / abs(m) if m else float("nan")
                print(f"  {k:10s} {m:14.6e} +/- {s:10.3e}  ({100*rel:.3f}%)")

        if len(allstats) == 2:
            a, b = allstats["sobol"], allstats["random"]
            print("\nEffective photon-number gain (random RMS / sobol RMS)^2")
            print("  >1 means the quasi-random run matches a larger ordinary run")
            for k in sorted(set(a) & set(b)):
                ra, rb = a[k][1], b[k][1]
                bias = abs(a[k][0] - b[k][0])
                sig = np.hypot(ra, rb)
                if ra > 0:
                    print(f"  {k:10s} gain = {(rb/ra)**2:8.2f}"
                          f"   (mean difference {bias/sig if sig else 0:.2f} sigma)")
            print("\nA mean difference of more than ~3 sigma would indicate a bias;"
                  "\nlaunch-only RQMC should show none.")
        return 0
    finally:
        if args.keep:
            print(f"\nreplicate outputs kept in {workdir}")
        else:
            shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
