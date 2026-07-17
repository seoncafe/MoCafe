#!/usr/bin/env python3
"""Check the composed runs: the output total luminosity equals L_int + L_ext.

L_int and L_ext are echoed by MoCafe at run time; here we re-derive L_ext from
the closed form pi * J * A_surface (sphere) and confirm TOT_LUM = L_int + L_ext.
"""
import numpy as np
import h5py

def totlum(fn):
    with h5py.File(fn, "r") as f:
        return float(f["Direct/"].attrs["TOT_LUM"][0])

# compose_mono: distance_unit='' (code units), rmax=1, J=0.025, L_int=1.0
Lext_mono = np.pi * 0.025 * 4.0 * np.pi * (1.0) ** 2
print("compose_mono.h5")
print(f"  L_int + L_ext (expected) = {1.0 + Lext_mono:.6e}")
print(f"  TOT_LUM (MoCafe)         = {totlum('compose_mono.h5'):.6e}")

# compose_sed: distance_unit='kpc', rmax=1, J=1e-3, L_int=3e41
KPC = 3.0856776e21
Lext_sed = np.pi * 1.0e-3 * 4.0 * np.pi * (1.0 * KPC) ** 2
print("compose_sed.h5")
print(f"  L_int + L_ext (expected) = {3.0e41 + Lext_sed:.6e}")
print(f"  TOT_LUM (MoCafe)         = {totlum('compose_sed.h5'):.6e}")
