#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

EXEC=../../MoCafe.x
#HOST=mocafe
#num_cores=88
HOST=lart4
num_cores=72

#mpirun -hosts $HOST -ppn $num_cores $EXEC uniform_taug100_u.in
#mpirun -hosts $HOST -ppn $num_cores $EXEC uniform_taug100_g.in
#mpirun -hosts $HOST -ppn $num_cores $EXEC uniform_taug100_r.in
#mpirun -hosts $HOST -ppn $num_cores $EXEC uniform_taug100_i.in
#mpirun -hosts $HOST -ppn $num_cores $EXEC uniform_taug100_z.in
#
#mpirun -hosts $HOST -ppn $num_cores $EXEC M020_001_taug100_u.in
#mpirun -hosts $HOST -ppn $num_cores $EXEC M020_001_taug100_g.in
#mpirun -hosts $HOST -ppn $num_cores $EXEC M020_001_taug100_r.in
#mpirun -hosts $HOST -ppn $num_cores $EXEC M020_001_taug100_i.in
#mpirun -hosts $HOST -ppn $num_cores $EXEC M020_001_taug100_z.in

mpirun -hosts $HOST -ppn $num_cores $EXEC M100_001_taug100_a100_g090.in
