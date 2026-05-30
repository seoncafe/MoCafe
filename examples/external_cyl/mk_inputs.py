#!/usr/bin/env python
import os
import socket
import glob
import numpy as np

def tau_band(taug=1.0,band='r', data_dir='../../data/'):
  cext_ref  = np.loadtxt(data_dir+'mueller_SDSS_g.dat',       usecols=1,skiprows=1,max_rows=1)
  cext_band = np.loadtxt(data_dir+'mueller_SDSS_'+band+'.dat',usecols=1,skiprows=1,max_rows=1)
  return taug*cext_band/cext_ref

def get_dust_params(band='r',data_dir='../../data/'):
  albedo = np.loadtxt(data_dir+'mueller_SDSS_'+band+'.dat',usecols=2,skiprows=1,max_rows=1)
  hgg    = np.loadtxt(data_dir+'mueller_SDSS_'+band+'.dat',usecols=3,skiprows=1,max_rows=1)
  return albedo, hgg

def write_input(mach=1.0,taug=1.0,band='r',centering=0,nsample=1, no_photons=1e8):
  source_geometry = 'external_rec'
  #use_stokes      = '.true.'
  use_stokes      = '.false.'
  sightline_tau   = '.true.'
  scatt_mat_file  = '../../data/mueller_SDSS_'+band+'.dat'
  no_print        = no_photons/10.0

  if mach == 0:
     density_file = ''
     file_in      = "uniform_taug%03d_%s" % (taug*100, band) + ".in"
  else:
     density_file = '../../b040_fbm3d_128/M%03db040_%03d.fits.gz' % (mach*10, nsample)
     file_in      = "M%03d_%03d_taug%03d_%s" % (mach*10, nsample, taug*100, band) + ".in"

  tau         = tau_band(taug=taug,band=band)
  albedo, hgg = get_dust_params(band=band)
  #hgg    = 0.53689
  #albedo = 0.67312

  xmax = 1.0
  ymax = 1.0
  zmax = 1.0
  nx   = 128
  ny   = 128
  nz   = 128
  nxim = 128
  nyim = 128
  distance = xmax*1e3

  f = open(file_in,'w')
  f.write("&parameters\n")
  f.write(" par%%no_photons     = %0.2e\n" % no_photons)
  f.write(" par%%albedo         = %7.5f\n" %albedo)
  f.write(" par%%hgg            = %7.5f\n" %hgg)
  f.write(" par%%scatt_mat_file = '%s'\n" %scatt_mat_file)
  f.write(" par%%use_stokes     = %s\n" %use_stokes)
  f.write(" par%%sightline_tau  = %s\n" %sightline_tau)
  f.write(" par%%tauhomo        = %0.2e\n" % tau)
  if mach > 0:
     f.write(" par%%density_file    = '%s'\n" % density_file)
  f.write(" par%%centering       = %d\n" % centering)
  f.write(" par%%source_geometry = '%s'\n" % source_geometry)
  f.write("\n")
  f.write(" par%distance_unit   = ''\n")
  f.write(" par%%xmax    = %5.1f\n" %xmax)
  f.write(" par%%ymax    = %5.1f\n" %ymax)
  f.write(" par%%zmax    = %5.1f\n" %zmax)
  f.write(" par%%nx = %d\n" %nx)
  f.write(" par%%ny = %d\n" %ny)
  f.write(" par%%nz = %d\n" %nz)
  f.write("\n")
  f.write(" par%%nxim = %d\n" %nxim)
  f.write(" par%%nyim = %d\n" %nyim)
  f.write(" par%inclination_angle = 0.0\n")
  f.write(" par%position_angle    = 0.0\n")
  f.write(" par%phase_angle       = 0.0\n")
  f.write(" par%%distance          = %0.1e\n" % distance)
  f.write("\n")
  f.write(" par%use_master_slave  = .true.\n")
  f.write(" par%%no_print = %0.2e\n" % no_print)
  f.write("/\n")
  f.close()
  return file_in

#---------------------------------------------
#mach_list = np.array([0.5, 1.0, 2.0, 3.0, 4.0, 6.0])
#taug_list = np.array([0.1, 0.5, 0.7, 0.8, 1.0, 2.0])
mach_list = np.array([0.0, 2.0])
taug_list = np.array([1.0])
band_list = ['u','g','r','i','z']
nsample_list = np.arange(1)+1
#nsample_list = np.arange(10)+1
#---------------------------------------------
hostname = socket.gethostname()
num_cores = 72
if hostname == 'mocafe': num_cores = 88
#---------------------------------------------
file_sh  = 'run.sh'
fsh      = open(file_sh,'w')
fsh.write("#!/bin/bash\n")
fsh.write('exec < /dev/null 2>&1\n')
fsh.write('trap "" HUP\n')
fsh.write("\n")
fsh.write("EXEC=../../MoCafe.x\n")
fsh.write("HOST=%s\n" % hostname)
fsh.write("num_cores=%d\n" % num_cores)
fsh.write("\n")
for nsample in nsample_list:
  for mach in mach_list:
    for taug in taug_list:
      for band in band_list:
        file_in = write_input(mach=mach,taug=taug,band=band,nsample=nsample)
        fsh.write("mpirun -hosts $HOST -ppn $num_cores $EXEC %s\n" % (file_in))
fsh.close()
os.chmod(file_sh, 0o755)
