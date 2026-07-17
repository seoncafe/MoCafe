#!/bin/bash
#--- Compose an internal source with an isotropic external field, in one run.
#---   compose_mono.in : monochromatic (scattering-only)
#---   compose_sed.in  : SED (wavelength-resolved, hot star + cool external field)
#--- Outputs: compose_mono.h5, compose_sed.h5  (checked by check_compose.py)
EXEC=../../MoCafe.x
UNAME=`uname`
case $UNAME in
'Linux')  NTHREADS=`nproc --all` ;;
'Darwin') NTHREADS=`sysctl -n hw.ncpu` ;;
*)        NTHREADS=8 ;;
esac
echo "Running $EXEC with $NTHREADS processes."
mpirun -np $NTHREADS $EXEC compose_mono.in
mpirun -np $NTHREADS $EXEC compose_sed.in
python3 check_compose.py
