#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

class read_L1642_data():
   def __init__(self):
      pos, u, u_err1, u_err2, b, b_err1, b_err2, y, y_err1, y_err2, I385, I385_err1, I385_err2, I415, I415_err1, I415_err2 \
                = np.loadtxt('Final2020JUL9N.txt',skiprows=2,usecols=(0,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19), unpack=True, max_rows=36)
      w         = np.argsort(pos)
      pos1      = np.int32(pos[w])
      self.u         = u[w]
      self.u_err1    = u_err1[w]
      self.u_err2    = u_err2[w]
      self.b         = b[w]
      self.b_err1    = b_err1[w]
      self.b_err2    = b_err2[w]
      self.y         = y[w]
      self.y_err1    = y_err1[w]
      self.y_err2    = y_err2[w]
      self.I385      = I385[w]
      self.I385_err1 = I385_err1[w]
      self.I385_err2 = I385_err2[w]
      self.I415      = I415[w]
      self.I415_err1 = I415_err1[w]
      self.I415_err2 = I415_err2[w]
   
      pos, tau200, Au, Au_err, A384, A384_err, A415, A415_err, Ab, Ab_err, Ay, Ay_err \
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
      self.A384     = A384[w]
      self.A384_err = A384_err[w]
      self.A415     = A415[w]
      self.A415_err = A415_err[w]
