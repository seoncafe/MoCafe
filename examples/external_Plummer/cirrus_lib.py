#!/usr/bin/env python
# 2024.04.16, K.I. Seon
import numpy as np
from astropy import units as u
from astropy.modeling import models

def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(3).mean(1)

def ISRF(wavl, per_Hz=False):
   try:
      wavl.value
   except:
      wavl = wavl * u.AA

   T_arr = [3000.0, 4000.0, 7500.0]
   W_arr = [7e-13, 1.65e-13, 1e-14]
   n     = len(T_arr)
   for i in np.arange(n):
      bb = models.BlackBody(T_arr[i]*u.K)
      if i==0:
         flux = bb(wavl) * W_arr[i]
      else:
         flux = flux + bb(wavl) * W_arr[i]

   # the output unit: erg cm^-2 s^-1 sr^-1 Hz^-1
   if per_Hz == False:
      # unit: erg cm^-2 s^-1 sr^-1 A^-1
      freq = wavl.to(u.Hz, equivalencies=u.spectral())
      flux = flux * freq / wavl
   return flux
