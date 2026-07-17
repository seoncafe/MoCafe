#!/usr/bin/env python3
"""Write an ABSOLUTE external-field spectrum for the sed_external_abs case.

The file is in par%spectrum_type = 'per_um' units, but for an external field
column 2 is a mean-intensity density J_lambda [erg/s/cm^2/sr/um] (not a
luminosity density).  Its integral over the wavelength grid is the band mean
intensity J [erg/s/cm^2/sr], from which MoCafe derives the entering
luminosity pi * J * A_surface.  Here J_lambda is a 1e4 K blackbody scaled so
that the band integral over [LAM_MIN, LAM_MAX] equals J_TARGET.
"""
import numpy as np

LAM_MIN, LAM_MAX = 0.0912, 1000.0
NLAM     = 400
T        = 1.0e4        # K
J_TARGET = 1.0e-3       # erg/s/cm^2/sr, band mean intensity

hc_k = 1.43877687750393e4   # h c / k in [um K]
lam = np.logspace(np.log10(LAM_MIN), np.log10(LAM_MAX), NLAM)   # um
x = hc_k / (lam * T)
blam = np.where(x < 700.0, 1.0 / (lam**5 * np.expm1(np.minimum(x, 700.0))), 0.0)

jint = np.trapz(blam, lam)
jlam = blam * (J_TARGET / jint)

with open("ext_spectrum_per_um.dat", "w") as f:
    f.write("# External-field spectrum, 1e4 K blackbody shape.\n")
    f.write(f"# Scaled so band J over [{LAM_MIN}, {LAM_MAX}] um = {J_TARGET:.3e} erg/s/cm^2/sr.\n")
    f.write("# lambda[um]   J_lambda[erg/s/cm^2/sr/um]   (par%spectrum_type = 'per_um')\n")
    for l, j in zip(lam, jlam):
        f.write(f"  {l:.6e}   {j:.6e}\n")

print(f"ext_spectrum_per_um.dat: band J = {J_TARGET:.3e} erg/s/cm^2/sr over [{LAM_MIN}, {LAM_MAX}] um")
