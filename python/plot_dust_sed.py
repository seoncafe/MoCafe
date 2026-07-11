#!/usr/bin/env python3
"""Plot the dust-emission SED from a MoCafe ``_dustsed`` output.

Overlays the intrinsic (unattenuated) and emergent (line-of-sight) dust SEDs as
lambda*L_lambda.  Works for Lucy-mode runs (par%use_dustemis with the 'lucy'
method).

Usage:
    python plot_dust_sed.py <run>            # <run> is a stem or a _dustsed file
    python plot_dust_sed.py <run> -o out.pdf
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import mocafe_io as mio


def resolve(run, suffix):
    """Accept a direct file path or a run stem; return the output file path."""
    if os.path.isfile(run):
        return run
    f = mio.find_mocafe_file(run, suffix)
    if f is None:
        raise SystemExit(f"no {suffix} output found for {run!r}")
    return f


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run", help="run stem or a _dustsed file")
    ap.add_argument("-o", "--out", default=None, help="output PDF (default: <stem>_sed.pdf)")
    ap.add_argument("--obs", type=int, default=0, help="observer index for the emergent SED")
    ap.add_argument("--emergent", action="store_true",
                    help="overlay the emergent SED shape (rescaled to the intrinsic "
                         "total; the emergent absolute flux needs the observer "
                         "distance, which the _dustsed file does not store)")
    args = ap.parse_args()

    path = resolve(args.run, "_dustsed")
    mf = mio.load_mocafe(path)
    d = mf.dust_sed()
    if d is None:
        raise SystemExit(f"{path} has no dust SED (SED_intrinsic / Wavelength)")

    lam = d["wavelength"]
    dlam = d["dwavelength"]
    intr = d["intrinsic"]                       # L_lam * dlam per bin [erg/s]
    lLl_intr = lam * np.where(dlam > 0, intr / dlam, 0.0)

    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax.plot(lam, lLl_intr, color="#333333", lw=2.0, label="intrinsic")

    emg = d["emergent"]
    if args.emergent and emg is not None:
        emg = np.atleast_2d(emg)
        k = min(args.obs, emg.shape[0] - 1)
        lLl_emg = lam * np.where(dlam > 0, emg[k] / dlam, 0.0)
        # rescale the shape to the intrinsic total (the emergent absolute flux
        # would need the observer distance, which is not in the _dustsed file).
        s_emg = np.trapz(lLl_emg, np.log(lam))
        s_int = np.trapz(lLl_intr, np.log(lam))
        if s_emg > 0:
            lLl_emg = lLl_emg * (s_int / s_emg)
        ax.plot(lam, lLl_emg, color="#d62728", lw=1.6,
                label=f"emergent (obs {k}, shape)")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(max(1.0, lam.min()), lam.max())
    ymax = max(lLl_intr.max(), 1e-30)
    ax.set_ylim(ymax * 1e-4, ymax * 2.0)
    ax.set_xlabel(r"$\lambda\ [\mu\mathrm{m}]$")
    ax.set_ylabel(r"$\lambda L_\lambda\ [\mathrm{erg\,s^{-1}}]$")
    ax.legend(frameon=False, loc="upper left")
    ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)

    out = args.out or (mio._strip_mocafe_ext(os.path.basename(path)) + "_sed.pdf")
    fig.tight_layout()
    fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
