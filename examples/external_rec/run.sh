#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

EXEC=../../MoCafe.x
HOST=mocafe
num_cores=88

mpirun -hosts $HOST -ppn $num_cores $EXEC uniform_taug100_g_cub111.in

mpirun -hosts $HOST -ppn $num_cores $EXEC uniform_taug100_g.in
mpirun -hosts $HOST -ppn $num_cores $EXEC uniform_taug100_z.in
mpirun -hosts $HOST -ppn $num_cores $EXEC uniform_taug100_u.in
mpirun -hosts $HOST -ppn $num_cores $EXEC uniform_taug100_r.in
mpirun -hosts $HOST -ppn $num_cores $EXEC uniform_taug100_i.in

mpirun -hosts $HOST -ppn $num_cores $EXEC M020_001_taug100_u.in
mpirun -hosts $HOST -ppn $num_cores $EXEC M020_001_taug100_g.in
mpirun -hosts $HOST -ppn $num_cores $EXEC M020_001_taug100_r.in
mpirun -hosts $HOST -ppn $num_cores $EXEC M020_001_taug100_i.in
mpirun -hosts $HOST -ppn $num_cores $EXEC M020_001_taug100_z.in
