#!/usr/bin/env python3
"""Validate the external-field SED runs.

For the absolute case, cross-check the derived luminosity against the closed
form pi * J * A_surface with A_surface = 4 pi (rmax * kpc)^2, and confirm the
emergent direct SED follows the input Planck shape.
"""
import numpy as np
import h5py

KPC = 3.0856776e21   # cm
RMAX = 1.0           # kpc
J_BAND = 1.0e-3      # erg/s/cm^2/sr (set in make_ext_spectrum.py)

A_surf = 4.0 * np.pi * (RMAX * KPC) ** 2
lum_expected = np.pi * J_BAND * A_surf
print(f"A_surface           = {A_surf:.6e} cm^2")
print(f"pi*J*A (expected L) = {lum_expected:.6e} erg/s")

for fn in ("sed_external.h5", "sed_external_abs.h5"):
    with h5py.File(fn, "r") as f:
        totlum = float(f["Direct/"].attrs["TOT_LUM"][0])
        wave   = f["Wavelength/data"][:]
        direc0 = f["Direct0/data"][:]
    sed = direc0.sum(axis=(1, 2))
    print(f"\n{fn}")
    print(f"  TOT_LUM (MoCafe)  = {totlum:.6e} erg/s")
    if fn == "sed_external_abs.h5":
        rel = abs(totlum - lum_expected) / lum_expected
        print(f"  relative error    = {rel:.3e}  (derived J -> luminosity)")
    print(f"  direct SED nonzero bins: {int((sed > 0).sum())}/{sed.size}, "
          f"peak at lambda = {wave[sed.argmax()]:.3f} um")
