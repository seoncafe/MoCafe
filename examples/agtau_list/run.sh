#!/bin/bash
#--- Run the three (a,g)/tau scan demo cases.
#--- Outputs: ag_only.h5, tau_only.h5, both.h5  (analyzed by show_results.ipynb)
EXEC=../../MoCafe.x
UNAME=`uname`
case $UNAME in
'Linux')  NTHREADS=`nproc --all` ;;
'Darwin') NTHREADS=`sysctl -n hw.ncpu` ;;
*)        NTHREADS=8 ;;
esac
echo "Running $EXEC with $NTHREADS processes."

mpirun -np $NTHREADS $EXEC ag_only.in
mpirun -np $NTHREADS $EXEC tau_only.in
mpirun -np $NTHREADS $EXEC both.in
mpirun -np $NTHREADS $EXEC single_tau.in   # ground truth for the tau-scan validation
