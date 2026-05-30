#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from read_L1642_data import *
from cirrus_lib import *

u1 = fits.open('Plummer_u1_obs.fits.gz')
t1 = fits.open('Plummer_u1_obs_tau.fits.gz')
u2 = fits.open('Plummer_u2_obs.fits.gz')
t2 = fits.open('Plummer_u2_obs_tau.fits.gz')

dist   = 124.0
radius = 1.528
pc     = 3.0857e18
t2a    = 2.5*np.log10(np.exp(1.0))

scale = 1.2e8
wavl  = 3500.0
Iu    = ISRF(wavl)
Lum_u = 4.0*np.pi*(radius*pc)**2*(np.pi*Iu.value) * scale

# Wavelengths
# u,      b,      y,      385,    415
# 3500.0, 4670.0, 5550.0, 3840.0, 4160.0
plt.scatter(t1[0].data*t2a, u1[0].data*Lum_u, s=0.5, color='blue')
plt.scatter(t2[0].data*t2a, u2[0].data*Lum_u, s=0.5, color='red')

a = read_L1642_data()
plt.scatter(a.Au, a.Iu, color='black')
plt.xlim(0.0, 13.0)
plt.ylim(0.0, 24.0)
plt.ylabel(r'I(u) [10$^{-9}$ erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$ \AA$^{-1}$] ')
plt.xlabel(r'A$_u$ [mag]')
plt.show()
