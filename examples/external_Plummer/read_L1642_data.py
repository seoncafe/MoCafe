#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

class read_L1642_data:
   def __init__(self, band=None):
      pos, u, u_err1, u_err2, b, b_err1, b_err2, y, y_err1, y_err2, I385, I385_err1, I385_err2, I415, I415_err1, I415_err2 \
                = np.loadtxt('Final2020JUL9N.txt',skiprows=2,usecols=(0,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19), unpack=True, max_rows=36)
      w         = np.argsort(pos)
      pos1      = np.int32(pos[w])
      self.Iu        = u[w]
      self.Iu_err1   = u_err1[w]
      self.Iu_err2   = u_err2[w]
      self.Ib        = b[w]
      self.Ib_err1   = b_err1[w]
      self.Ib_err2   = b_err2[w]
      self.Iy        = y[w]
      self.Iy_err1   = y_err1[w]
      self.Iy_err2   = y_err2[w]
      self.I385      = I385[w]
      self.I385_err1 = I385_err1[w]
      self.I385_err2 = I385_err2[w]
      self.I415      = I415[w]
      self.I415_err1 = I415_err1[w]
      self.I415_err2 = I415_err2[w]
   
      pos, tau200, Au, Au_err, A385, A385_err, A415, A415_err, Ab, Ab_err, Ay, Ay_err \
            = np.loadtxt('L1642_AV+5colour+titleNEW.txt',skiprows=5,usecols=(0,11,17,18,19,20,21,22,23,24,25,26), unpack=True, max_rows=36)
      w             = np.argsort(pos)
      pos2          = np.int32(pos[w])
      self.tau200   = tau200[w]
      self.Au       = Au[w]
      self.Au_err   = Au_err[w]
      self.Ab       = Ab[w]
      self.Ab_err   = Ab_err[w]
      self.Ay       = Ay[w]
      self.Ay_err   = Ay_err[w]
      self.A385     = A385[w]
      self.A385_err = A385_err[w]
      self.A415     = A415[w]
      self.A415_err = A415_err[w]

      if band == 'u':
         self.Aobs      = self.Au
         self.Aobs_err  = self.Au_err
         self.Iobs      = self.Iu
         self.Iobs_err1 = self.Iu_err1
         self.Iobs_err2 = self.Iu_err2
      elif band == 'b':
         self.Aobs      = self.Ab
         self.Aobs_err  = self.Ab_err
         self.Iobs      = self.Ib
         self.Iobs_err1 = self.Ib_err1
         self.Iobs_err2 = self.Ib_err2
      elif band == 'y':
         self.Aobs      = self.Ay
         self.Aobs_err  = self.Ay_err
         self.Iobs      = self.Iy
         self.Iobs_err1 = self.Iy_err1
         self.Iobs_err2 = self.Iy_err2
      elif band == '385':
         self.Aobs      = self.A385
         self.Aobs_err  = self.A385_err
         self.Iobs      = self.I385
         self.Iobs_err1 = self.I385_err1
         self.Iobs_err2 = self.I385_err2
      elif band == '415':
         self.Aobs      = self.A415
         self.Aobs_err  = self.A415_err
         self.Iobs      = self.I415
         self.Iobs_err1 = self.I415_err1
         self.Iobs_err2 = self.I415_err2
