#!/bin/bash
#--- Multiple stellar sources with absolute ('per_um') spectrum files:
#--- the source luminosities are derived from the file integrals.
#--- Output: sed_multi.h5 (validated by plot_sed_multi.py)
EXEC=../../MoCafe.x
UNAME=`uname`
case $UNAME in
'Linux')  NTHREADS=`nproc --all` ;;
'Darwin') NTHREADS=`sysctl -n hw.ncpu` ;;
*)        NTHREADS=8 ;;
esac
echo "Running $EXEC with $NTHREADS processes."

python3 make_abs_spectra.py
mpirun -np $NTHREADS $EXEC sed_multi.in
python3 plot_sed_multi.py
