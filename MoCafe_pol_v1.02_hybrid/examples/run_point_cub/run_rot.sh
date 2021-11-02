#!/bin/bash

EXEC=../MoCafe_pol.x
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

mpirun -np $NTHREADS $EXEC cub111_rot000_tau50_hgg00_a10.in
mpirun -np $NTHREADS $EXEC cub111_rot045_tau50_hgg00_a10.in
mpirun -np $NTHREADS $EXEC cub111_rot090_tau50_hgg00_a10.in
mpirun -np $NTHREADS $EXEC cub111_rot135_tau50_hgg00_a10.in
mpirun -np $NTHREADS $EXEC cub111_rot180_tau50_hgg00_a10.in
