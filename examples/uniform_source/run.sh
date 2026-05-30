#!/bin/bash

EXEC=../../MoCafe.x
HOST=`hostname`
UNAME=`uname`

case $UNAME in
'Linux')
   NTHREADS=`nproc --all`
   ;;
'Darwin')
   NTHREADS=`sysctl -n hw.ncpu`
   ;;
*)
   NTHREADS=16
   ;;
esac
echo ""
echo "Running $EXEC on $HOST with $NTHREADS threads."


mpirun -np $NTHREADS $EXEC uniform_tau10_hgg02_a05.in
