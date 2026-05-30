#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

fname = 'uniform_taug100_r'
fname = 'test'
fname = 'M020_001_taug100_g'

a = fits.open(fname+'_obs.fits.gz')
fig, axs = plt.subplots(ncols=2,nrows=2,figsize=(9,9))
axs = axs.flat

ax = axs[0]
im = ax.imshow(a[0].data + a[1].data)
fig.colorbar(im,ax=ax)

factor = 4.0*np.pi
factor = 2.0

w           = np.where(a[2].data <= 0.0)
bkg         = np.median(a[2].data)*factor
im_cloud    = a[0].data + a[1].data*factor + bkg
im_cloud[w] = 2.0*bkg

ax = axs[1]
im = ax.imshow(im_cloud)
fig.colorbar(im,ax=ax)

im_cloud[:,:]  = im_cloud[:,:] - 2.0*bkg
np.mean(im_cloud[w])

ax = axs[2]
im = ax.imshow(im_cloud, vmin=0.0)
fig.colorbar(im,ax=ax)

plt.show()
