#!/usr/bin/env python3
"""Validate and inspect the MoCafe AMR dust-sphere examples (run after run.sh).

Prints integrated-flux comparisons of the AMR runs against the homogeneous
Cartesian reference, and saves a scattered-light montage (amr_montage.png).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits

CASES = [("smooth_ref", "Cartesian (ref)"),
         ("amr_uniform", "AMR uniform L5"),
         ("amr_refined", "AMR refined L3-6")]


def sec(stem, s):
    fn = stem + ".fits.gz"
    if not os.path.exists(fn):
        return None
    with fits.open(fn) as h:
        for hd in h:
            if hd.name == s and hd.data is not None:
                return np.array(hd.data, float)
    return None


def validate():
    ref = {s: sec("smooth_ref", s) for s in ("Scattered", "Direct", "Direct0")}
    if ref["Scattered"] is None:
        print("validate: run smooth_ref first.")
        return
    for stem, _ in CASES[1:]:
        if sec(stem, "Scattered") is None:
            continue
        print(f"=== {stem} vs smooth_ref ===")
        for s in ("Scattered", "Direct", "Direct0"):
            a, b = ref[s].sum(), sec(stem, s).sum()
            rel = (b - a) / a * 100 if a != 0 else 0.0
            print(f"  {s:10s} ref={a:.5e} amr={b:.5e} rel={rel:+.3f}%")


def logimg(ax, img, title):
    m = img > 0
    vmax = img[m].max() if m.any() else 1.0
    data = np.ma.masked_where(~m, img)
    cmap = plt.cm.inferno.copy()
    cmap.set_bad("white")
    ax.imshow(np.log10(np.clip(data, vmax*1e-4, vmax)), origin="lower",
              cmap=cmap, interpolation="nearest")
    ax.set_title(title, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])


def montage():
    avail = [(s, t) for s, t in CASES if sec(s, "Scattered") is not None]
    if not avail:
        print("montage: no outputs; run run.sh first.")
        return
    fig, axes = plt.subplots(1, len(avail), figsize=(3.0*len(avail), 3.2))
    if len(avail) == 1:
        axes = [axes]
    for ax, (stem, title) in zip(axes, avail):
        logimg(ax, sec(stem, "Scattered"), title)
    fig.suptitle("AMR dust sphere: scattered light (log scale)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig("amr_montage.png", dpi=120)
    print("montage: wrote amr_montage.png")


if __name__ == "__main__":
    validate()
    montage()
