#!/usr/bin/env python
from astropy.io import fits
from astropy import units as u
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import glob
import os.path

flist = glob.glob('*_obs.fits.gz')
nf    = len(flist)

for idx in np.arange(nf):
   # peeled-off intensity image
   fname = flist[idx]
   hdu   = fits.open(fname)
   hdr   = hdu[0].header
   dx    = hdr['CD1_1']
   dy    = hdr['CD2_2']
   dist  = hdr['DISTANCE']
   try:
      dxfreq = hdr['DXFREQ']
   except:
      dxfreq = 1.0
   try:
      distance_cm = hdr['DIST_CM']
      dist        = dist * distance_cm
   except:
      distance_cm = 1.0

   try:
      source_geometry = hdr['SRC_GEOM']
   except:
      source_geometry = 'point'

   #if source_geometry == 'external_sph':
   #   rmax   = hdr['XMAX'] * distance_cm
   #   cosvt0 = rmax / dist
   #   lum_fac = (1.0+cosvt0)/(2.0*np.pi*3.0*cosvt0**2)*( 2.0*cosvt0**2 - (1.0-cosvt0)*(1.0 - np.sqrt(1.0-cosvt0**2)) )
   #else:
   #   lum_fac = 1.0

   im1   = hdu[0].data
   im2   = hdu[1].data
   im    = im1 + im2
   try:
      im_direc0 = hdu[2].data
   except:
      pass
   hdu.close()

   scale = 4.0*np.pi * (dist * np.tan(dx * np.pi/180.0)) * (dist * np.tan(dy * np.pi/180.0)) * dxfreq
   ftot  = im.sum()  * scale
   fsca  = im1.sum() * scale
   fdir  = im2.sum() * scale
   try:
      fdir0 = im_direc0.sum() * scale
   except:
      fdir0 = -999.9

   print('')
   print(fname)
   print('   TOTAL  : %10.7f' % ftot)
   print('   SCATT  : %10.7f' % fsca)
   print('   DIRECT : %10.7f' % fdir)
   print('   DIRECT0: %10.7f' % fdir0)
   #if (lum_fac != 1.0):
   #   print('   LUM_FAC: %10.7f' % lum_fac)

   # peeled-off Stokes I image
   jj     = fname.find('_obs.fits.gz')
   basis  = fname[0:jj]
   fname2 = basis+'_stokes.fits.gz'
   if os.path.exists(fname2):
      hdu   = fits.open(fname2)
      im    = hdu[0].data
      hdu.close()

      scale = 4.0*np.pi * (dist * np.tan(dx * np.pi/180.0)) * (dist * np.tan(dy * np.pi/180.0)) * dxfreq
      ftot  = im.sum()  * scale

      print('')
      print(fname2)
      print('   TOTAL : %10.7f' % ftot)
