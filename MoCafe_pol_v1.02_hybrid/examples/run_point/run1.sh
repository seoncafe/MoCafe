#!/bin/bash

EXEC=../../MoCafe_pol.x
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

NTHREADS=1

mpirun -np $NTHREADS $EXEC point_tau10_R.in
mpirun -np $NTHREADS $EXEC point_tau10_hgg00_a10.in
#mpirun -np $NTHREADS $EXEC point_tau10_hgg02_a10.in
#mpirun -np $NTHREADS $EXEC point_tau10_hgg02_a05.in

#mpirun -np $NTHREADS $EXEC b040_point_tau50_hgg02_a05.in

#mpirun -np $NTHREADS $EXEC rec111_rot000_tau50_hgg00_a10.in
#mpirun -np $NTHREADS $EXEC rec111_rot045_tau50_hgg00_a10.in
#mpirun -np $NTHREADS $EXEC rec111_rot090_tau50_hgg00_a10.in

#---- We need to use the following command if more CPUs are neeeded.
#---- However, for the present purpose, using a single machine is fast enough.
#HOST=all_hosts
#mpirun -machinefile $HOST $EXEC point_tau10_g02_a05.in
