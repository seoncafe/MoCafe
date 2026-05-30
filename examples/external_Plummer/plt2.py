#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from read_L1642_data import *
from cirrus_lib import *

plt.rcParams.update({'font.size': 16})

u1 = fits.open('Plummer_u1_obs.fits.gz')
t1 = fits.open('Plummer_u1_obs_tau.fits.gz')
u2 = fits.open('Plummer_u2_obs.fits.gz')
t2 = fits.open('Plummer_u2_obs_tau.fits.gz')

t2a = 2.5*np.log10(np.exp(1.0))
A1  = t1[0].data * t2a
A2  = t2[0].data * t2a

#fscal = 0.2
fscal = 0.0
Iu1 = u1[0].data + (u1[1].data - u1[2].data) * fscal
Iu2 = u2[0].data + (u2[1].data - u2[2].data) * fscal

dist    = 124.0
radius  = 1.528
pc      = 3.0857e18

Lscale1 = 0.24
Lscale2 = 0.22
if fscal == 0.0:
  Lscale1 = 0.13
  Lscale2 = 0.18

wavl    = 3500.0
Iu      = ISRF(wavl)
Lum_u0  = 4.0*np.pi*(radius*pc)**2*(np.pi*Iu.value)
Lum_u   = 4.0*np.pi*(radius*pc)**2*(np.pi*Iu.value) * 1e9

#read L1642 data
a     = read_L1642_data()

ndat  = a.Au.size
for i in np.arange(ndat):
   w      = np.where(np.abs(A1 - a.Au[i]) <= a.Au_err[i])[0]
   u1_avg = np.mean(u1[0].data[w]*Lum_u0)
   u1_sig = np.std(u1[0].data[w]*Lum_u0)
   print('%5.2f %10.4e %10.4e %3d' % (a.Au[i], u1_avg, u1_sig, w.size))

# Wavelengths
# u,      b,      y,      385,    415
# 3500.0, 4670.0, 5550.0, 3840.0, 4160.0

fig, axs = plt.subplots(ncols=2,nrows=1, figsize=(14, 7))
x2 = 10.0

ax = axs[0]
ax.scatter(A1, Iu1*Lum_u*Lscale1, s=0.5, color='blue', label='isotropic')
ax.scatter(a.Au, a.u, color='black')
ax.legend(markerscale=5.0, frameon=False, handletextpad=0.05)
ax.set_xlim(0.0, x2)
ax.set_ylim(0.0, 24.0)
ax.set_ylabel(r'I(u) [10$^{-9}$ erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$ \AA$^{-1}$] ')
ax.set_xlabel(r'A$_u$ [mag]')

ax = axs[1]
ax.scatter(A2, Iu2*Lum_u*Lscale2, s=0.5, color='red',  label='anisotropic')
ax.scatter(a.Au, a.u, color='black')
ax.legend(markerscale=5.0, frameon=False, handletextpad=0.05)
ax.set_xlim(0.0, x2)
ax.set_ylim(0.0, 24.0)
ax.set_ylabel(r'I(u) [10$^{-9}$ erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$ \AA$^{-1}$] ')
ax.set_xlabel(r'A$_u$ [mag]')

plt.show()
