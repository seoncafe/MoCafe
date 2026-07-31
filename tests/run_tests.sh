#!/bin/bash
# Build and run the wavelength-sampling checks.
#
#   cd tests && ./run_tests.sh
#
# These cover the two places a photon's wavelength is drawn, both of which used
# to hand back the center of a wavelength bin:
#
#   test_spectrum_sampler   the source spectra (spectrum_sampler_mod): a
#                           piecewise-linear luminosity density, its analytic
#                           cumulative distribution, and the exact inversion.
#                           Also reports how far the old bin-center scheme was
#                           from the exact band integral of a blackbody.
#
#   test_power_law_bin      the draw inside one bin used by the dust emission
#                           (sample_power_law_bin, logarithmic_mean): a spectrum
#                           read as a power law across the bin, its integral, and
#                           the closed-form inversion, including the degenerate
#                           ends (equal values, a zero end, a ratio near one).
#
# Both are self-checking: every line prints the quantity and what it must equal.
# The objects in ../src are built with -ipo and cannot be linked here, so the two
# modules these need are compiled fresh into this directory.
set -e
cd "$(dirname "$0")"

FC=${FC:-mpiifort}
if ! command -v $FC >/dev/null; then FC=mpif90; fi
FLAGS="-O2 -cpp -DMPI"
case "$FC" in
   mpiifort|mpiifx|ifort|ifx) MODFLAG="-module ." ;;
   *)                         MODFLAG="-J . -I ." ;;
esac

echo "compiling with $FC ..."
$FC $FLAGS $MODFLAG -c ../src/define.f90                -o define.o
$FC $FLAGS $MODFLAG -c ../src/spectrum_sampler_mod.f90  -o spectrum_sampler_mod.o

for t in test_spectrum_sampler test_power_law_bin; do
   $FC $FLAGS $MODFLAG $t.f90 spectrum_sampler_mod.o define.o -o $t.x
done

for t in test_spectrum_sampler test_power_law_bin; do
   echo
   ./$t.x
done

rm -f *.o *.mod
