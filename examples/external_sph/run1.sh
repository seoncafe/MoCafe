#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

EXEC=../../MoCafe.x
HOST=lart4
num_cores=72

mpirun -hosts $HOST -ppn $num_cores $EXEC test_taug100_r.in
mpirun -hosts $HOST -ppn $num_cores $EXEC test_taug100_r0.in
