#!/usr/bin/env python3
"""Consolidated validation / regression checks for the MoCafe dust-emission
suite.  Run from examples/benchmarks after the bench_* runs (make_bench.sh):

    python validate.py

Checks, each with a PASS/FAIL tolerance:
  1. Energy conservation (Lucy): the intrinsic dust luminosity equals the
     absorbed power exactly (the emitted spectrum is normalized to the locally
     absorbed power cell by cell).  The emergent SED is smaller at high optical
     depth because the dust emission is partly reabsorbed -- that is physics,
     not an error, so it is not checked here.
  2. Two-method agreement: the absorbed = emitted luminosity from the Lucy
     path-length estimator and the Bjorkman & Wood immediate-reemission scheme
     agree across an optical-depth sweep.
  3. SEDust emission vs the SHG multi-code benchmark (CRT / DIRTY / SKIRT /
     DustEM / TRADING / MCFOST / DARTRAY), when the reference data is present.
"""
import glob
import os
import re
import sys
import numpy as np

sys.path.insert(0, os.path.abspath("../../python"))
import mocafe_io as mio

HERE = os.path.dirname(os.path.abspath(__file__))
results = []   # (name, ok, detail)


def check(name, ok, detail):
    results.append((name, ok, detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")


# ---- 1. energy conservation on every _dustsed run ------------------------
print("1. Energy conservation (Lucy):")
for f in sorted(glob.glob(os.path.join(HERE, "*_dustsed.h5"))):
    mf = mio.load_mocafe(f)
    d = mf.dust_sed()
    if d is None:
        continue
    sec = mf.section("SED_emergent")
    L_abs = float(sec.attr("L_ABS")) if sec and sec.attr("L_ABS") else None
    L_int = float(np.sum(d["intrinsic"]))
    tag = os.path.basename(f).replace("_dustsed.h5", "")
    if L_abs:
        r1 = abs(L_int - L_abs) / L_abs
        check(f"{tag} intrinsic dust L = absorbed L", r1 < 1e-3, f"rel = {r1:.1e}")


# ---- 2. Lucy vs Bjorkman & Wood across the tau sweep ---------------------
print("2. Lucy vs Bjorkman & Wood absorbed luminosity:")


def labs_from_log(path):
    if not os.path.isfile(path):
        return None
    txt = open(path).read()
    m = re.search(r"total cell emission \(absorbed\) L\s*:\s*([0-9.Ee+\-]+)", txt)
    if m:
        return float(m.group(1))
    m = re.search(r"total absorbed L\s*=\s*([0-9.Ee+\-]+)", txt)
    return float(m.group(1)) if m else None


for lucy in sorted(glob.glob(os.path.join(HERE, "bench_tau*_lucy.log"))):
    tau = re.search(r"tau([0-9.]+)_lucy", lucy).group(1)
    bw = lucy.replace("_lucy.log", "_bw01.log")
    Ll, Lb = labs_from_log(lucy), labs_from_log(bw)
    if Ll and Lb:
        r = abs(Ll - Lb) / Ll
        check(f"tau={tau} Lucy vs BW01", r < 0.02, f"|dL|/L = {r:.2%}")


# ---- 3. SEDust vs the SHG multi-code benchmark ---------------------------
print("3. SEDust emission vs SHG multi-code benchmark (incl. SKIRT / DIRTY):")
REF = os.path.expanduser("~/MoCafe/Grain/SHG_Benchmark/Results_FullSolution")
CODES = ["CRT", "DIRTY", "SKIRT", "DustEM", "TRADING", "MCFOST", "DARTRAY"]


def loadfix(f):
    rows = []
    for line in open(f):
        line = line.strip()
        if not line or line[0] == "#":
            continue
        vals = []
        for t in line.split():
            t2 = re.sub(r"^([+-]?\d*\.?\d+)([+-]\d{2,3})$", r"\1e\2", t)
            try:
                vals.append(float(t2))
            except ValueError:
                vals.append(0.0)
        rows.append(vals)
    return np.array(rows)


def norm(lam, e):
    i = e > 0
    A = np.trapz(e[i], np.log(lam[i]))
    return e / A if A > 0 else e


if os.path.isdir(REF):
    for Utag, key in [("1e-02", "1.E-02"), ("1e+00", "1.E+00"),
                      ("1e+02", "1.E+02"), ("1e+04", "1.E+04")]:
        ours = glob.glob(os.path.join(HERE, f"sedust_zubko_Mathis_U_*{key}*.dat"))
        if not ours:
            continue
        od = loadfix(ours[0]); ol, oe = od[:, 0], od[:, 1]
        curves = []
        for c in CODES:
            fp = f"{REF}/{c}_SHG_Mathis_U_{Utag}.dat"
            if os.path.exists(fp):
                d = loadfix(fp); curves.append((d[:, 0], d[:, 2:].sum(axis=1)))
        if not curves:
            continue
        lam = curves[0][0]; oe_i = np.interp(lam, ol, oe)
        med = np.median([norm(lam, c[1]) for c in curves], 0)
        on = norm(lam, oe_i)
        m = (lam > 3) & (lam < 1000) & (med > 0.02 * med.max())
        rel = np.median(np.abs(on[m] - med[m]) / med[m])
        check(f"SHG U={Utag}", rel < 0.06,
              f"median rel.diff vs {len(curves)}-code median = {rel:.1%}")
else:
    print(f"  [skip] reference not found at {REF}")


# ---- summary -------------------------------------------------------------
npass = sum(ok for _, ok, _ in results)
print(f"\n{npass}/{len(results)} checks passed.")
sys.exit(0 if npass == len(results) else 1)
