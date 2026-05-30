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


mpirun -np $NTHREADS $EXEC point_tau20_hgg05_a05.in
mpirun -np $NTHREADS $EXEC point_tau01_hgg00_a05.in
mpirun -np $NTHREADS $EXEC rec_tau50_hgg05_a05.in

#mpirun -np $NTHREADS $EXEC point_tau10_R.in
#mpirun -np $NTHREADS $EXEC point_tau10_hgg00_a10.in

#mpirun -np $NTHREADS $EXEC point_tau10_hgg02_a10.in
#mpirun -np $NTHREADS $EXEC point_tau10_hgg02_a05.in

#mpirun -np $NTHREADS $EXEC b040_point_tau50_hgg02_a05.in

#---- We need to use the following command if more CPUs are neeeded.
#---- However, for the present purpose, using a single machine is fast enough.
#HOST=all_hosts
#mpirun -machinefile $HOST $EXEC point_tau10_g02_a05.in
