#!/bin/bash
#--- Monochromatic multiple internal sources (two point sources).
#--- Output: mono_multi.h5  (TOT_LUM = sum of src_lum)
EXEC=../../MoCafe.x
UNAME=`uname`
case $UNAME in
'Linux')  NTHREADS=`nproc --all` ;;
'Darwin') NTHREADS=`sysctl -n hw.ncpu` ;;
*)        NTHREADS=8 ;;
esac
echo "Running $EXEC with $NTHREADS processes."
mpirun -np $NTHREADS $EXEC mono_multi.in
