#!/usr/bin/env python3
"""Compare the intrinsic dust-emission SED of the three SEDust grain models
(astrodust, DL07, Zubko) as computed by MoCafe.

Reads the `_dustsed` outputs of three otherwise-identical runs that differ only
in par%dust_model_sed, and overlays lambda*L_lambda.  Each run absorbs the same
stellar power, so the total emitted luminosity is identical and the curves
differ only in spectral shape (PAH bands, silicate features, FIR peak, submm
slope).

Regenerate the inputs with model_compare_{astrodust,dl07,zubko}.in:
    for m in astrodust dl07 zubko; do
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


fig, ax = plt.subplots(figsize=(7.0, 4.6))
for model, label, color in MODELS:
    lam, lLl = load_sed(model)
    ax.plot(lam, lLl, color=color, lw=1.8, label=label)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1.0, 2.0e3)
ymax = ax.get_ylim()[1]
ax.set_ylim(ymax * 1e-4, ymax)
ax.set_xlabel(r"$\lambda\ [\mu\mathrm{m}]$")
ax.set_ylabel(r"$\lambda L_\lambda\ [\mathrm{erg\,s^{-1}}]$")
ax.legend(frameon=False, loc="upper left")
ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)

out = os.path.join(HERE, "dust_model_compare.pdf")
fig.tight_layout()
fig.savefig(out)
print("wrote", out)
