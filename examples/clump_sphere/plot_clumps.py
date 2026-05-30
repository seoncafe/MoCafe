#!/usr/bin/env python3
"""Inspect and validate the MoCafe clump examples.

  * Prints a validation table comparing the single-centered-giant-clump run
    (clump_giant) against the homogeneous-sphere reference (smooth_ref):
    a correct clump ray-tracer reproduces the smooth result to MC noise.
  * Saves a montage (clump_montage.png) of the scattered-light images for all
    available cases.

Run after run.sh.  Usage:  python3 plot_clumps.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits

CASES = [
    ("smooth_ref",     "smooth sphere (ref)"),
    ("clump_giant",    "giant clump (= ref)"),
    ("clump_basic",    "uniform clumps"),
    ("clump_powerlaw", "density profile"),
    ("clump_cone",     "biconical"),
]


def load_section(stem, section):
    fn = stem + ".fits.gz"
    if not os.path.exists(fn):
        return None
    with fits.open(fn) as h:
        for hdu in h:
            if hdu.name == section and hdu.data is not None:
                return np.array(hdu.data, dtype=float)
    return None


def validate():
    a = {s: load_section("smooth_ref", s) for s in ("Scattered", "Direct", "Direct0")}
    b = {s: load_section("clump_giant", s) for s in ("Scattered", "Direct", "Direct0")}
    if a["Scattered"] is None or b["Scattered"] is None:
        print("validate: need smooth_ref and clump_giant outputs; skipping.")
        return
    print("=== Validation: centered giant clump vs homogeneous sphere ===")
    print(f"{'section':12s} {'smooth':>14s} {'giant':>14s} {'rel diff':>10s}")
    for s in ("Scattered", "Direct", "Direct0"):
        sa, sb = a[s].sum(), b[s].sum()
        rel = (sb - sa) / sa * 100 if sa != 0 else 0.0
        print(f"{s:12s} {sa:14.6e} {sb:14.6e} {rel:+9.3f}%")


def logimg(ax, img, title):
    m = img > 0
    vmax = img[m].max() if m.any() else 1.0
    vmin = vmax * 1e-4
    data = np.ma.masked_where(~m, img)
    cmap = plt.cm.inferno.copy()
    cmap.set_bad("white")
    ax.imshow(np.log10(np.clip(data, vmin, vmax)), origin="lower",
              cmap=cmap, interpolation="nearest")
    ax.set_title(title, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])


def montage():
    avail = [(s, t) for s, t in CASES if load_section(s, "Scattered") is not None]
    if not avail:
        print("montage: no outputs found; run run.sh first.")
        return
    n = len(avail)
    fig, axes = plt.subplots(1, n, figsize=(3.0 * n, 3.2))
    if n == 1:
        axes = [axes]
    for ax, (stem, title) in zip(axes, avail):
        logimg(ax, load_section(stem, "Scattered"), title)
    fig.suptitle("Scattered light (log scale)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig("clump_montage.png", dpi=120)
    print("montage: wrote clump_montage.png")


if __name__ == "__main__":
    validate()
    montage()
