#!/usr/bin/env python
import os
import socket
import glob
import numpy as np

nsample   = 1
band      = 'K'
centering = 1
#source_geometry = 'external'
source_geometry = 'uniform'

no_photons = 1e8
use_stokes = '.true.'
peel_tau   = '.true.'
no_print   = no_photons/10.0

if band == 'K':
   hgg     = 0.13095
   albedo  = 0.43864
   scatt_mat_file  =  '../data/mueller_K.dat'
if band == 'V':
   hgg     = 0.53689
   albedo  = 0.67312
   scatt_mat_file='../data/mueller_V.dat'

xmax = 1.0
ymax = 1.0
zmax = 1.0
rmax = 1.0
nx   = 128
ny   = 128
nz   = 128
nxim = 128
nyim = 128
distance = xmax*1e3

mach_list = np.array([0.0, 0.5, 1.0, 2.0, 3.0])
tau_list  = np.array([0.01, 0.1, 1.0])

hostname = socket.gethostname()
#-----------
file_log = 'log.txt'
file_sh  = 'run.sh'
NTHREADS = 160
fsh      = open(file_sh,'w')
fsh.write("#!/bin/bash\n")
fsh.write("\n")
fsh.write("EXEC=../MoCafe.x\n")
fsh.write("HOST=%s_only\n" % hostname)
fsh.write("\n")
#-----------
for mach in mach_list:
  for tau in tau_list:
    file_in      = "M%03d_%03d%s_t%03d" % (mach*10, nsample, band, tau*100) + ".in"
    if mach == 0:
       density_file = ''
    else:
       density_file = '../b040_fbm3d_128/M%03db040_%03d.fits.gz' % (mach*10, nsample)

    f = open(file_in,'w')
    f.write("&parameters\n")
    a = int(np.log10(no_photons))
    b = no_photons/10.0**a
    f.write(" par%%no_photons     = %0.1fe%d\n" % (b,a))
    f.write(" par%%hgg            = %7.5f\n" %hgg)
    f.write(" par%%albedo         = %7.5f\n" %albedo)
    f.write(" par%%scatt_mat_file = '%s'\n" %scatt_mat_file)
    f.write(" par%%use_stokes     = %s\n" %use_stokes)
    f.write(" par%%peel_tau       = %s\n" %peel_tau)
    a = int(np.log10(tau))
    b = tau/10.0**a
    f.write(" par%%taumax          = %0.1fe%d\n" % (b,a))
    f.write(" par%%density_file    = '%s'\n" % density_file)
    f.write(" par%%centering       = %d\n" % centering)
    f.write(" par%%source_geometry = '%s'\n" % source_geometry)
    f.write("\n")
    f.write(" par%distance_unit   = ''\n")
    f.write(" par%%xmax    = %5.1f\n" %xmax)
    f.write(" par%%ymax    = %5.1f\n" %ymax)
    f.write(" par%%zmax    = %5.1f\n" %zmax)
    f.write(" par%%rmax    = %5.1f\n" %rmax)
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
    f.write(" par%iseed    = 0\n")
    a = int(np.log10(no_print))
    b = no_print/10.0**a
    f.write(" par%%no_print = %0.1fe%d\n" % (b,a))
    f.write("/\n")
    f.close()

    fsh.write("mpirun -machinefile $HOST $EXEC %s\n" % (file_in))
fsh.close()

hostfile = open(hostname+'_only','w')
if hostname == 'mocafe':
   hostfile.write("mocafe:88\n")
if hostname == 'lart1':
   hostfile.write("lart1:72\n")
if hostname == 'lart2':
   hostfile.write("lart2:72\n")
if hostname == 'lart3':
   hostfile.write("lart3:72\n")
hostfile.close()
