#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from read_L1642_data import *
from cirrus_lib import *
from scipy.optimize import minimize

plt.rcParams.update({'font.size': 16})

#read L1642 data
dat   = read_L1642_data()
fname = 'Plummer_u1'

def model_out(x, fname, band):
   # cloud radius
   radius  = 1.528
   pc      = 3.0857e18

   f_diff = x[0]
   L_scal = x[1]
   mod    = fits.open(fname+'_obs.fits.gz')
   tau    = fits.open(fname+'_obs_tau.fits.gz')

   # Wavelengths
   # u,      b,      y,      385,    415
   # 3500.0, 4670.0, 5550.0, 3840.0, 4160.0
   if band == 'u':
      wavl     = 3500.0
   elif band == 'b':
      wavl     = 4670.0
   elif band == 'y':
      wavl     = 5550.0
   elif band == '385':
      wavl     = 3840.0
   elif band == '415':
      wavl     = 4160.0

   I_radiation = ISRF(wavl)
   Lum         = 4.0*np.pi*(radius*pc)**2*(np.pi*I_radiation.value) * 1e9 * L_scal

   tau2a = 2.5*np.log10(np.exp(1.0))
   Amod  = tau[0].data * tau2a
   Imod  = (mod[0].data + (mod[1].data - mod[2].data) * f_diff) * Lum

   return Amod, Imod

def chisqrt(x, arg_1, arg_2):
   # cloud radius
   radius  = 1.528
   pc      = 3.0857e18

   f_diff = x[0]
   L_scal = x[1]
   fname  = arg_1
   band   = arg_2
   mod    = fits.open(fname+'_obs.fits.gz')
   tau    = fits.open(fname+'_obs_tau.fits.gz')

   # Wavelengths
   # u,      b,      y,      385,    415
   # 3500.0, 4670.0, 5550.0, 3840.0, 4160.0
   if band == 'u':
      wavl     = 3500.0
      Aobs     = dat.Au
      Aobs_err = dat.Au_err
      Iobs     = dat.Iu
      Iobs_err = (dat.Iu_err1 + dat.Iu_err2)/2.0
      #Iobs_err = np.sqrt(dat.Iu_err1**2 + dat.Iu_err2**2)/2.0
   elif band == 'b':
      wavl     = 4670.0
      Aobs     = dat.Ab
      Aobs_err = dat.Ab_err
      Iobs     = dat.Ib
      Iobs_err = (dat.Ib_err1 + dat.Ib_err2)/2.0
      #Iobs_err = np.sqrt(dat.Ib_err1**2 + dat.Ib_err2**2)/2.0
   elif band == 'y':
      wavl     = 5550.0
      Aobs     = dat.Ay
      Aobs_err = dat.Ay_err
      Iobs     = dat.Iy
      Iobs_err = (dat.Iy_err1 + dat.Iy_err2)/2.0
      #Iobs_err = np.sqrt(dat.Iy_err1**2 + dat.Iy_err2**2)/2.0
   elif band == '385':
      wavl     = 3840.0
      Aobs     = dat.A385
      Aobs_err = dat.A385_err
      Iobs     = dat.I385
      Iobs_err = (dat.I385_err1 + dat.I385_err2)/2.0
      #Iobs_err = np.sqrt(dat.I385_err1**2 + dat.I385_err2**2)/2.0
   elif band == '415':
      wavl     = 4160.0
      Aobs     = dat.A415
      Aobs_err = dat.A415_err
      Iobs     = dat.I415
      Iobs_err = (dat.I415_err1 + dat.I415_err2)/2.0
      #Iobs_err = np.sqrt(dat.I415_err1**2 + dat.I415_err2**2)/2.0

   I_radiation = ISRF(wavl)
   Lum         = 4.0*np.pi*(radius*pc)**2*(np.pi*I_radiation.value) * 1e9 * L_scal

   tau2a = 2.5*np.log10(np.exp(1.0))
   Amod  = tau[0].data * tau2a
   Imod  = (mod[0].data + (mod[1].data - mod[2].data) * f_diff) * Lum

   Imod_avg = np.zeros(Aobs.size)
   Imod_sig = np.zeros(Aobs.size)
   for i in np.arange(Aobs.size):
      #w           = np.where(np.abs(Amod - Aobs[i]) <= Aobs_err[i])[0]
      w           = np.where(np.abs(Amod - Aobs[i]) <= Aobs_err[i]*2.0)[0]
      Imod_avg[i] = np.mean(Imod[w])
      Imod_sig[i] = np.std(Imod[w])
      #print('%5.2f %10.4e %10.4e %3d' % (Aobs[i], Imod_avg[i], Imod_sig[i], w.size))

   #chi2 = np.sum((Iobs - Imod_avg)**2/(Iobs_err**2 + Imod_sig**2)) / (Iobs.size - 2.0)
   #chi2 = np.sum((Iobs - Imod_avg)**2/(Iobs_err**2 + Imod_sig**2*0.01)) / (Iobs.size - 2.0)
   chi2 = np.sum((Iobs - Imod_avg)**2/Iobs_err**2) / (Iobs.size - 2.0)

   #chi2 = np.sum((Iobs - Imod_avg)**2/Imod_sig**2) / (Iobs.size - 2.0)
   #chi2 = np.sum((Iobs - Imod_avg)**2/(Iobs_err**2*0.5 + Imod_sig**2)) / (Iobs.size - 2.0)
   #-----
   w    = np.where(Aobs < 3.0)[0]
   chi2 = np.sum((Iobs[w] - Imod_avg[w])**2/Iobs_err[w]**2) / (Iobs[w].size - 2.0)
   return chi2

#x0   = [0.2, 0.2]
#bnds = ((0.0, 0.5), (0.0, 0.8))
#bnds = ((0.0, None), (0.0, None))
x0   = [0.2, 0.2]
bnds1 = ((0.0, 0.2), (0.0, None))
bnds2 = ((0.0, None), (0.0, None))

fig, axs = plt.subplots(ncols=2,nrows=1, figsize=(14, 7))
x2 = 13.0

#-----------------------------------
method='Nelder-Mead'
#method='SLSQP'
#method='Powell'
#-----------------------------------
ax = axs[0]
result = minimize(chisqrt, x0, args = ('Plummer_u1', 'u'), bounds=bnds1, method=method)
Amod, Imod = model_out(result.x, 'Plummer_u1', 'u')
ax.scatter(Amod, Imod, s=0.5, color='red', label='scattering only')

#x0 = result.x
result = minimize(chisqrt, x0, args = ('Plummer_u1', 'u'), bounds=bnds2, method=method)
Amod, Imod = model_out(result.x, 'Plummer_u1', 'u')
ax.scatter(Amod, Imod, s=0.5, color='green', label='scatt. + direc')
print(result)

ax.errorbar(dat.Au, dat.Iu, [dat.Iu_err1, dat.Iu_err2], fmt='o', color='black')
ax.legend(markerscale=5.0, frameon=False, handletextpad=0.05)

ax.set_title('isotropic radiation')
ax.set_xlim(0.0, x2)
ax.set_ylim(0.0, 24.0)
ax.set_ylabel(r'I(u) [10$^{-9}$ erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$ \AA$^{-1}$] ')
ax.set_xlabel(r'A$_u$ [mag]')

#-----------------------------------
ax = axs[1]
result = minimize(chisqrt, x0, args = ('Plummer_u2', 'u'), bounds=bnds1, method=method)
Amod, Imod = model_out(result.x, 'Plummer_u2', 'u')
ax.scatter(Amod, Imod, s=0.5, color='red', label='scattering only')

#x0 = result.x
result = minimize(chisqrt, x0, args = ('Plummer_u2', 'u'), bounds=bnds2, method=method)
Amod, Imod = model_out(result.x, 'Plummer_u2', 'u')
ax.scatter(Amod, Imod, s=0.5, color='green', label='scatt. + direc')
print(result)

#ax.scatter(dat.Au, dat.Iu, color='black')
ax.errorbar(dat.Au, dat.Iu, [dat.Iu_err1, dat.Iu_err2], fmt='o', color='black')
ax.legend(markerscale=5.0, frameon=False, handletextpad=0.05)

ax.set_title('anisotropic radiation')
ax.set_xlim(0.0, x2)
ax.set_ylim(0.0, 24.0)
ax.set_ylabel(r'I(u) [10$^{-9}$ erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$ \AA$^{-1}$] ')
ax.set_xlabel(r'A$_u$ [mag]')

plt.show()
