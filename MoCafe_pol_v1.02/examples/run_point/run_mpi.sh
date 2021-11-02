#!/bin/bash

EXEC=../../MoCafe_pol.x

#---- We need to use the following command if more CPUs are neeeded.
#---- However, for the present purpose, using a single machine is fast enough.
HOST=all_hosts
mpirun -machinefile $HOST $EXEC point_tau10_R.in
#mpirun -machinefile $HOST $EXEC point_tau10_hgg00_a10.in
#mpirun -machinefile $HOST $EXEC point_tau10_hgg02_a10.in
#mpirun -machinefile $HOST $EXEC point_tau10_hgg02_a05.in
#b040_point_tau50_hgg02_a05.in
#point_tau10_R.in
#point_tau10_hgg00_a10.in
#point_tau10_hgg02_a05.in
#point_tau10_hgg02_a10.in
