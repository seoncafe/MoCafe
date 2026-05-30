from radial_profile import *
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

fname = 'uniform_taug100_g_111'
a=fits.open(fname+'_obs.fits.gz')

im0, im1, im2 = a[0].data, a[1].data, a[2].data

r, p0 = radial_profile(im0)
r, p1 = radial_profile(im1)
r, p2 = radial_profile(im2)

plt.plot(r, p0+p1, color='black')
plt.plot(r, p0+p1-p2, color='black', linestyle='dashed')
plt.plot(r, p0, color='blue')
plt.plot(r, p1, color='red')
plt.plot(r, p2, color='red', linestyle='dotted')
