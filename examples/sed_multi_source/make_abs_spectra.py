#!/usr/bin/env python3
"""Build ABSOLUTE source spectra for the sed_multi_source example.

Reads the FSPS shape spectra in ../../data/ (lambda[um], L_lambda[arb]) and
rescales each so that its integral over the RT wavelength range
[LAM_MIN, LAM_MAX] equals a target luminosity [erg/s].  The output files are
in par%spectrum_type = 'per_um' units: lambda[um], L_lambda[erg/s/um].
With par%src_lum unset, MoCafe derives each source luminosity from these
files (the derived values differ from the targets only by the wavelength-bin
discretization of the RT grid).
"""
import numpy as np

LAM_MIN, LAM_MAX = 0.0912, 2000.0   # match par%lambda_min/lambda_max

CASES = [
    ("../../data/spec_young_fsps.dat", "spec_young_per_um.dat", 3.0e36,
     "FSPS SSP 10 Myr (young), rescaled to L = 3.0e36 erg/s over [0.0912, 2000] um"),
    ("../../data/spec_old_fsps.dat",   "spec_old_per_um.dat",   1.0e36,
     "FSPS SSP 10 Gyr (old), rescaled to L = 1.0e36 erg/s over [0.0912, 2000] um"),
]

for src, out, ltarget, note in CASES:
    lam, flam = np.loadtxt(src, unpack=True)
    m = (lam >= LAM_MIN) & (lam <= LAM_MAX)
    lint = np.trapz(flam[m], lam[m])          # [arb erg/s]
    scale = ltarget / lint
    with open(out, "w") as f:
        f.write(f"# {note}\n")
        f.write("# lambda[um]   L_lambda[erg/s/um]   (par%spectrum_type = 'per_um')\n")
        for l, y in zip(lam, flam * scale):
            f.write(f"  {l:.6e}   {y:.6e}\n")
    print(f"{out}: L({LAM_MIN}-{LAM_MAX} um) = {ltarget:.3e} erg/s  (scale {scale:.6e})")
