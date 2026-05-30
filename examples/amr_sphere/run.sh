#!/bin/bash
#--- AMR octree dust-sphere examples for MoCafe (AMR_CLUMPS_PLAN.md Part A).
#
#   smooth_ref    homogeneous dust sphere on the Cartesian grid (32^3) -- ref
#   amr_uniform   uniform-density sphere on a uniform-level-5 octree (= 32^3);
#                 must reproduce smooth_ref (validates the octree ray-tracer)
#   amr_refined   same uniform sphere on a radially-refined octree (levels
#                 3-6); exercises ray traversal across refinement boundaries
#
# Inspect/compare with:  python3 plot_amr.py
EXEC=../../MoCafe.x
case `uname` in
'Linux')  NTHREADS=`nproc --all` ;;
'Darwin') NTHREADS=`sysctl -n hw.ncpu` ;;
*)        NTHREADS=8 ;;
esac
if [ "$NTHREADS" -gt 16 ]; then NTHREADS=16; fi
echo "Running $EXEC with $NTHREADS processes."

#--- build the generic AMR grids read by the inputs
python3 ../../python/AMR_grid/make_amr_sphere.py --level 5 --out amr_sphere_uniform.fits
python3 ../../python/AMR_grid/make_amr_sphere.py --level-min 3 --level-max 6 --out amr_sphere_refined.fits

mpirun -np $NTHREADS $EXEC smooth_ref.in
mpirun -np $NTHREADS $EXEC amr_uniform.in
mpirun -np $NTHREADS $EXEC amr_refined.in
