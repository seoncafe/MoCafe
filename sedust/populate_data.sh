#!/bin/bash
# Populate the vendored SEDust data files needed by MoCafe's dust-emission
# stage.  These are bulky (~55 MB) optics tables, so they are kept out of git
# (see .gitignore) and copied here from the canonical SEDust tree instead.
# Re-run once after a fresh checkout on a machine where SEDust is available.
set -e
SED=${SEDUST_SRC:-/nfs/mocafe/kiseon/MoCafe/Grain/SEDust}
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$HERE/data/dielectric" "$HERE/data/release" "$HERE/tmatrix/output"
cp "$SED"/data/dielectric/{DH21_aeff,index_CpaD03,index_CpeD03,index_silD03,PAHion.31,PAHneu.31,q_D16graphite.dat,qlib_gra_D16MGemt_1.400} "$HERE/data/dielectric/"
cp "$SED"/data/kext_astrodust_MW.dat "$HERE/data/"
cp "$SED"/data/release/{extinction.dat,kext_albedo_WD_MW_3.1_60_D03.all_2003,scattering.dat,size_distribution.dat} "$HERE/data/release/"
cp "$SED"/tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat "$HERE/tmatrix/output/"
echo "Populated vendored SEDust data under $HERE"
