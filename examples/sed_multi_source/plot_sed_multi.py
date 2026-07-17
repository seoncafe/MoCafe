#!/usr/bin/env python3
"""Validate the sed_multi_source example.

Checks that (1) the derived total luminosity equals the sum of the two file
integrals over the wavelength grid, and (2) the emergent unattenuated direct
SED (Direct0 summed over the image) reproduces the luminosity-weighted input
spectrum (Source_lum is written only for a single global spectrum, so the
reference here is built from the two files directly).
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt

FILES  = ["spec_young_per_um.dat", "spec_old_per_um.dat"]

with h5py.File("sed_multi.h5", "r") as f:
    wave   = f["Wavelength/data"][:]
    dwave  = f["Dwavelength/data"][:]
    direc0 = f["Direct0/data"][:]              # (nlam, ny, nx)
    totlum = float(f["Direct/"].attrs["TOT_LUM"][0])

sed_mc = direc0.sum(axis=(1, 2))
sed_mc /= sed_mc.sum()

# reference: bin-integrated (interp x dwave) sum of the two absolute files,
# following the code convention (log-lambda interpolation, zero outside file).
ref = np.zeros_like(wave)
ltot = 0.0
for fn in FILES:
    lam, flam = np.loadtxt(fn, unpack=True)
    y = np.interp(np.log(wave), np.log(lam), flam, left=np.nan, right=np.nan)
    y = np.where((wave < lam[0]) | (wave > lam[-1]), 0.0, np.nan_to_num(y))
    ref += y * dwave
    ltot += (y * dwave).sum()
ref_pdf = ref / ref.sum()

print(f"TOT_LUM (derived by MoCafe) : {totlum:.6e} erg/s")
print(f"sum of file bin-integrals   : {ltot:.6e} erg/s")
print(f"relative difference         : {abs(totlum-ltot)/ltot:.3e}")

fig, ax = plt.subplots(figsize=(6, 4))
ax.loglog(wave, ref_pdf/dwave, "k-", lw=1, label="input (young + old)")
ax.loglog(wave, sed_mc/dwave, "ro", ms=3, label=r"MC Direct0")
ax.set_xlabel(r"$\lambda$ [$\mu$m]")
ax.set_ylabel(r"normalized $L_\lambda$")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig("sed_multi.pdf")
print("wrote sed_multi.pdf")
