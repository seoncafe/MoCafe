#!/bin/bash
#--- Clumpy-dust-medium examples for MoCafe (AMR_CLUMPS_PLAN.md Part B).
#
#   smooth_ref      homogeneous dust sphere (Cartesian grid) -- reference
#   clump_giant     one centered clump filling the sphere (file-loaded);
#                   must reproduce smooth_ref (the ray-tracer validation)
#   clump_basic     uniform random non-overlapping clump population
#   clump_powerlaw  clumps with a centrally-concentrated dust-opacity profile
#   clump_cone      biconical clump distribution
#
# After running, compare/inspect with:  python3 plot_clumps.py
EXEC=../../MoCafe.x
HOST=`hostname`
UNAME=`uname`
case $UNAME in
'Linux')  NTHREADS=`nproc --all` ;;
'Darwin') NTHREADS=`sysctl -n hw.ncpu` ;;
*)        NTHREADS=8 ;;
esac
if [ "$NTHREADS" -gt 16 ]; then NTHREADS=16; fi
echo "Running $EXEC on $HOST with $NTHREADS processes."

#--- generate the single centered giant clump loaded by clump_giant.in
python3 ../../python/make_clumps.py --single --rmax 1.0 --clump-radius 1.0 \
        --tau0 1.0 --out clump_giant.fits

mpirun -np $NTHREADS $EXEC smooth_ref.in
mpirun -np $NTHREADS $EXEC clump_giant.in
mpirun -np $NTHREADS $EXEC clump_basic.in
mpirun -np $NTHREADS $EXEC clump_powerlaw.in
mpirun -np $NTHREADS $EXEC clump_cone.in
