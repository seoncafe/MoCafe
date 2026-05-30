#!/usr/bin/env python
#fname = 'Plummer_u1_obs.fits.gz'
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

fname = ['Plummer_u1a_obs.fits.gz','Plummer_u1b_obs.fits.gz','Plummer_u1c_obs.fits.gz']
n     = len(fname)
nside = [128, 256, 512]

fig, ax = plt.subplots(ncols=3,figsize=(18,5))

for i in np.arange(n):
   a1 = fits.open(fname[i])
   b1 = fits.open('../../../MoCafe_v1.18/examples/external_Plummer/'+fname[i])

   plt.axes(ax[i])
   plt.scatter(a1[0].data, b1[0].data,s=2)
   lim2 = np.amax(a1[0].data)*1.05
   plt.plot([0.0,lim2],[0.0,lim2],color='red')
   plt.axis('scaled')
   plt.xlim(0.0, lim2)
   plt.ylim(0.0, lim2)
   plt.title(r'nside = %d' % nside[i])
   plt.xlabel(r'method = 2')
   plt.ylabel(r'method = 1')

plt.show()
