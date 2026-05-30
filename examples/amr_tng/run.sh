#!/bin/bash
#--- Illustris-TNG dust-scattering example (dust-only AMR octree).
#
# Pipeline:  TNG gas cutout (HDF5) --convert--> generic AMR grid --> MoCafe.
# Requires a TNG PartType0 gas cutout HDF5.  Pass its path as $1, or set the
# CUTOUT environment variable, or download one via the TNG API (see README).
#
#   ./run.sh /path/to/cutout_TNG50-1_snap99_subNNN.hdf5
#
# Inspect afterwards:  python3 plot_tng.py  (or see tng_montage.png).
EXEC=../../MoCafe.x
CONV=../../python/AMR_grid/convert_illustris_to_generic.py
CUTOUT=${1:-${CUTOUT:-cutout.hdf5}}
case `uname` in
'Linux')  NTHREADS=`nproc --all` ;;
'Darwin') NTHREADS=`sysctl -n hw.ncpu` ;;
*)        NTHREADS=8 ;;
esac
if [ "$NTHREADS" -gt 16 ]; then NTHREADS=16; fi

if [ ! -f "$CUTOUT" ]; then
  echo "ERROR: TNG cutout '$CUTOUT' not found."
  echo "  Pass the cutout path as the first argument, or download one with"
  echo "  convert_illustris_to_generic.py --api-key KEY --simulation TNG50-1 ..."
  exit 1
fi

#--- (1) convert the cutout to generic AMR grids (uniform and adaptive)
python3 $CONV "$CUTOUT" -o tng_uniform.fits --grid-type uniform --level-max 5 \
        --boxsize 60 --emit-ndust
python3 $CONV "$CUTOUT" -o tng_amr.h5       --grid-type amr --level-min 3 \
        --level-max 6 --boxsize 60 --emit-ndust

#--- (2) run MoCafe (central source illuminating the galaxy dust)
mpirun -np $NTHREADS $EXEC tng.in        # uniform grid, dust_model=laursen09
mpirun -np $NTHREADS $EXEC tng_amr.in    # adaptive grid, dust_model=global_dgr
