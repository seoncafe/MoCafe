#!/bin/bash

HOST_FILE='all_hosts'
EXEC=../../MoCafe_pol.x
HOST=`hostname`
UNAME=`uname`

$EXEC point_tau10_R.in
$EXEC point_tau10_hgg00_a10.in

#mpirun -machinefile $HOST_FILE $EXEC point_tau10_hgg02_a10.in
#mpirun -machinefile $HOST_FILE $EXEC point_tau10_hgg02_a05.in

#mpirun -machinefile $HOST_FILE $EXEC b040_point_tau50_hgg02_a05.in

#mpirun -machinefile $HOST_FILE $EXEC rec111_rot000_tau50_hgg00_a10.in
#mpirun -machinefile $HOST_FILE $EXEC rec111_rot045_tau50_hgg00_a10.in
#mpirun -machinefile $HOST_FILE $EXEC rec111_rot090_tau50_hgg00_a10.in

#---- We need to use the following command if more CPUs are neeeded.
#---- However, for the present purpose, using a single machine is fast enough.
#HOST=all_hosts
#mpirun -machinefile $HOST $EXEC point_tau10_g02_a05.in
