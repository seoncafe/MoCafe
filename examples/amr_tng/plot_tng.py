#!/usr/bin/env python3
"""Montage of the TNG50 dust-scattering images (run after run.sh / tng.in)."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits


def sec(fn, s):
    with fits.open(fn) as h:
        for hd in h:
            if hd.name == s and hd.data is not None:
                return np.array(hd.data, float)
    return None


def panel(ax, img, title):
    m = img > 0
    vmax = img[m].max() if m.any() else 1.0
    data = np.ma.masked_where(~m, img)
    cmap = plt.cm.inferno.copy(); cmap.set_bad("white")
    ax.imshow(np.log10(np.clip(data, vmax*1e-4, vmax)), origin="lower",
              cmap=cmap, interpolation="nearest")
    ax.set_title(title, fontsize=10); ax.set_xticks([]); ax.set_yticks([])


fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.4))
for a, (s, t) in zip(ax, [("Scattered", "scattered"), ("Direct", "direct"),
                          ("Direct0", "direct0 (unattenuated)")]):
    img = sec("tng_image.fits.gz", s)
    if img is not None:
        panel(a, img, "TNG50 " + t)
fig.suptitle("TNG50 galaxy dust scattering (converter to MoCafe AMR)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("tng_montage.png", dpi=120)
print("wrote tng_montage.png")
