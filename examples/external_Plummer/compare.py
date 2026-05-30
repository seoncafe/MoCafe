#!/usr/bin/env python
#fname = 'Plummer_u1_obs.fits.gz'
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

#fname = 'Plummer_u2_obs.fits.gz'
fname = 'Plummer_u1_obs.fits.gz'
a1 = fits.open(fname)
b1 = fits.open('../../../MoCafe_v1.18/examples/external_Plummer/'+fname)

fig, (ax1,ax2) = plt.subplots(ncols=2,figsize=(12,5))

plt.axes(ax1)
plt.scatter(a1[0].data, b1[0].data,s=2)
lim2 = np.amax(a1[0].data)*1.02
plt.plot([0.0,lim2],[0.0,lim2],color='red')
plt.axis('scaled')
plt.xlim(0.0, lim2)
plt.ylim(0.0, lim2)

plt.axes(ax2)
plt.scatter(a1[1].data, b1[1].data,s=2)
lim2 = np.amax(a1[1].data)*1.02
plt.plot([0.0,lim2],[0.0,lim2],color='red')
plt.axis('scaled')
plt.xlim(0.0, lim2)
plt.ylim(0.0, lim2)

plt.show()
