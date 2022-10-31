#!/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from   square_subplots import square_subplots
from   scipy import interpolate
from   scipy.ndimage import gaussian_filter
from   astropy.io import fits
from   matplotlib.ticker import AutoMinorLocator, MultipleLocator
from   mpl_toolkits.axes_grid1.inset_locator import inset_axes
from   radial_profile import radial_profile

np.seterr(divide='ignore', invalid='ignore')

fname = '../../../MoCafe_pol_v1.02_hybrid/examples/run_point/point_tau10_R'
hdu   = fits.open(fname+'_stokes.fits.gz')

Iim   = hdu[0].data
Qim   = hdu[1].data
Uim   = hdu[2].data
Vim   = hdu[3].data

polim  = np.sqrt(Qim**2 + Uim**2)/Iim
r, pol = radial_profile(polim)
r, I   = radial_profile(Iim)

#-------
fname = '../../../MoCafe_pol_v1.02/examples/run_point/point_tau10_R'
hdu   = fits.open(fname+'_stokes.fits.gz')

Iim   = hdu[0].data
Qim   = hdu[1].data
Uim   = hdu[2].data
Vim   = hdu[3].data

polim  = np.sqrt(Qim**2 + Uim**2)/Iim
r2, pol2 = radial_profile(polim)
r2, I2   = radial_profile(Iim)

#-------
fig, ax = plt.subplots(1,2,figsize=(8.0,3.5))

ax[0].set_xlim(0,1)
ax[0].set_ylim(0,0.25)
ax[0].set_xlabel('r')
ax[0].set_ylabel('polarization')
ax[0].plot(r, pol)
ax[0].plot(r2, pol2, "o", markersize=1.5)

ax[1].set_xlim(0,1)
ax[1].set_ylim(0,0.25)
ax[1].set_xlabel('r')
ax[1].set_ylabel('intensity')
ax[1].plot(r, I)
ax[1].plot(r2, I2, "o", markersize=1.5)

plt.subplots_adjust(wspace=0.4, hspace=0.2, left=0.1,right=0.9,bottom=0.1,top=0.9)
s = square_subplots(fig)
plt.show()
