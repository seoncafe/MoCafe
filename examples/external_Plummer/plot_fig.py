#!/usr/bin/env python
import matplotlib.pyplot as plt
from astropy.io import fits

a     = fits.open('Plummer_tau20_R_obs.fits.gz')
a_tau = fits.open('Plummer_tau20_R_obs_tau.fits.gz')
Iscatt1 = a[0].data
Idirec1 = a[1].data

b     = fits.open('Plummer_tau20_R_angular_obs.fits.gz')
b_tau = fits.open('Plummer_tau20_R_angular_obs_tau.fits.gz')
Iscatt2 = b[0].data
Idirec2 = b[1].data

plt.subplot(121)
plt.scatter(a_tau[0].data, Iscatt1,s=0.5, color='black', label='Isotropic, Scattered')
plt.scatter(b_tau[0].data, Iscatt2,s=0.5, color='red',   label='AnnIsotropic, Scattered')
plt.legend(frameon=False)
#plt.xlim(0.0, 13.0)
plt.xlim(0.0, 8.0)
plt.ylim(0.0, 0.015)

plt.subplot(122)
plt.scatter(a_tau[0].data, Idirec1,s=0.3, color='black', label='Isotropic, Direct')
plt.scatter(b_tau[0].data, Idirec2,s=0.3, color='red',   label='AnnIsotropic, Direct')
plt.legend(frameon=False)
#plt.xlim(0.0, 13.0)
plt.xlim(0.0, 8.0)
plt.ylim(0.0, 0.041)

plt.show()
