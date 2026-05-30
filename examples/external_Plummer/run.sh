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


#mpirun -np $NTHREADS $EXEC Plummer_u1a.in
#mpirun -np $NTHREADS $EXEC Plummer_u1b.in
#mpirun -np $NTHREADS $EXEC Plummer_u1c.in
mpirun -np $NTHREADS $EXEC Plummer_u1.in
#mpirun -np $NTHREADS $EXEC Plummer_u2.in
#mpirun -np $NTHREADS $EXEC Plummer_tau20_R.in
#mpirun -np $NTHREADS $EXEC Plummer_tau20_R_angular.in

#---- We need to use the following command if more CPUs are neeeded.
#---- However, for the present purpose, using a single machine is fast enough.
