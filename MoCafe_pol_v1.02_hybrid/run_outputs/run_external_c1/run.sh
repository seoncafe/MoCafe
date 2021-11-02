#!/bin/bash

EXEC=../MoCafe.x
HOST=lart2_only

mpirun -machinefile $HOST $EXEC M000_001K_t001.in
mpirun -machinefile $HOST $EXEC M000_001K_t010.in
mpirun -machinefile $HOST $EXEC M000_001K_t100.in
mpirun -machinefile $HOST $EXEC M005_001K_t001.in
mpirun -machinefile $HOST $EXEC M005_001K_t010.in
mpirun -machinefile $HOST $EXEC M005_001K_t100.in
mpirun -machinefile $HOST $EXEC M010_001K_t001.in
mpirun -machinefile $HOST $EXEC M010_001K_t010.in
mpirun -machinefile $HOST $EXEC M010_001K_t100.in
mpirun -machinefile $HOST $EXEC M020_001K_t001.in
mpirun -machinefile $HOST $EXEC M020_001K_t010.in
mpirun -machinefile $HOST $EXEC M020_001K_t100.in
mpirun -machinefile $HOST $EXEC M030_001K_t001.in
mpirun -machinefile $HOST $EXEC M030_001K_t010.in
mpirun -machinefile $HOST $EXEC M030_001K_t100.in
