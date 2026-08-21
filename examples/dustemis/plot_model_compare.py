#!/usr/bin/env python3
"""Compare the intrinsic dust-emission SED of the four SEDust grain models
(astrodust, DL07, MRN, Zubko) as computed by MoCafe.

Reads the `_dustsed` outputs of four otherwise-identical runs that differ only
in par%dust_model, and overlays lambda*L_lambda.  Naming the model sets both
halves of its dust physics: each run transports through that model's own
extinction and reemits with that model's own grains.  The curves therefore
differ in level as well as in shape -- the models do not absorb the same
stellar power out of an identical cloud, because tau is normalized at
par%lambda_ref and the models part company away from it -- on top of the
spectral differences (PAH bands, silicate features, FIR peak, submm slope).

Regenerate the inputs with model_compare_{astrodust,dl07,mrn,zubko}.in:
    for m in astrodust dl07 mrn zubko; do
        mpirun -np 4 ../../MoCafe.x model_compare_$m.in
    done
then run this script.
"""
import os
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

MODELS = [
    ("astrodust", "Astrodust (HD23)", "#1f77b4"),
    ("dl07",      "DL07",             "#d62728"),
    ("mrn",       "MRN (DL84)",       "#ff7f0e"),
    ("zubko",     "Zubko (ZDA BARE-GR-S)", "#2ca02c"),
]


def load_sed(model):
    """Return (lambda [um], lambda*L_lambda [erg/s]) from a _dustsed file."""
    path = os.path.join(HERE, f"sed_{model}_dustsed.h5")
    with h5py.File(path, "r") as f:
        lam = f["/Wavelength/data"][:]        # bin centers [um]
        dlam = f["/Dwavelength/data"][:]       # bin widths [um]
        sed = f["/SED_intrinsic/data"][:]      # L_lambda * dlam per bin [erg/s]
    Llam = np.where(dlam > 0, sed / dlam, 0.0)  # erg/s/um
    return lam, lam * Llam                       # lambda*L_lambda [erg/s]


LAM_LO, LAM_HI = 1.0, 2.0e3

fig, ax = plt.subplots(figsize=(7.0, 4.6))
curves = []
for model, label, color in MODELS:
    lam, lLl = load_sed(model)
    ax.plot(lam, lLl, color=color, lw=1.8, label=label)
    curves.append((lam, lLl))

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(LAM_LO, LAM_HI)
# The top of the axis comes from the data inside the plotted range, not from the
# autoscale: MRN has no PAHs, so its near-infrared bins hold essentially no
# emission, and a log autoscale that has to accommodate them spans hundreds of
# decades and leaves the far-infrared peaks in a sliver.
ymax = max(lLl[(lam >= LAM_LO) & (lam <= LAM_HI)].max() for lam, lLl in curves)
ax.set_ylim(ymax * 2e-4, ymax * 2.0)
ax.set_xlabel(r"$\lambda\ [\mu\mathrm{m}]$")
ax.set_ylabel(r"$\lambda L_\lambda\ [\mathrm{erg\,s^{-1}}]$")
ax.legend(frameon=False, loc="upper left")
ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)

out = os.path.join(HERE, "dust_model_compare.pdf")
fig.tight_layout()
fig.savefig(out)
print("wrote", out)
