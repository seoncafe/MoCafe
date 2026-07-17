#!/bin/bash
#--- External illumination in SED mode.
#---   sed_external.in     : Planck field shape, absolute scale from ext_intensity
#---   sed_external_abs.in : absolute field spectrum file, J (hence luminosity) derived
#--- Outputs: sed_external.h5, sed_external_abs.h5  (checked by check_sed_external.py)
EXEC=../../MoCafe.x
UNAME=`uname`
case $UNAME in
'Linux')  NTHREADS=`nproc --all` ;;
'Darwin') NTHREADS=`sysctl -n hw.ncpu` ;;
*)        NTHREADS=8 ;;
esac
echo "Running $EXEC with $NTHREADS processes."

python3 make_ext_spectrum.py
mpirun -np $NTHREADS $EXEC sed_external.in
mpirun -np $NTHREADS $EXEC sed_external_abs.in
python3 check_sed_external.py
