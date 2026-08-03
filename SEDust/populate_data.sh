#!/bin/bash
# Populate SEDust optics data that is NOT committed to git.
#
# The minimal data set the default 'astrodust'+PAH dust-emission path needs is
# already in the repository (see .gitignore), so a fresh checkout runs
# dust emission out of the box -- you normally do NOT need this script.
#
# The only file kept out of git is the 17 MB D16 spheroid Q-library
# (qlib_gra_D16MGemt_1.400), read only by q_graphite_d16_mod, which nothing uses
# by default.  Run this script to fetch it if you enable that non-default
# graphite path.  It also re-copies the rest, so it doubles as a full refresh
# from a canonical SEDust tree.
#
# The sources under sed/src/ and sed/build_lib.sh are byte-identical to SEDust
# v1.00 revision 04b6ed8 (2026-08-03), so `cmp` against that tree is the whole
# check that this copy is current.  The one file of the library that reaches
# into the T-matrix, src/euv_astrodust_tmatrix.f90, is deliberately not copied:
# it serves the astrodust EUV band, which only a lam_min shorter than 0.0912 um
# asks for, and MoCafe never passes one.  Without it the archive links with no
# T-matrix at all, and euv_tmatrix = .true. is refused (status 6) instead of
# being answered with a different grain.
#
# >>> EDIT the path below to your own SEDust location, or run with
# >>>   SEDUST_SRC=/your/path/to/SEDust ./populate_data.sh
set -e
SED=${SEDUST_SRC:-/home/kiseon/MoCafe/Grain/SEDust_v1.00}
if [ ! -d "$SED" ]; then
  echo "ERROR: SEDust source tree not found at: $SED" >&2
  echo "       Set it to your own location, e.g.:" >&2
  echo "         SEDUST_SRC=/path/to/SEDust $0" >&2
  exit 1
fi
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$HERE/data/dielectric" "$HERE/data/release" "$HERE/data/zubko" "$HERE/tmatrix/output"
cp "$SED"/data/dielectric/{DH21_aeff,index_CpaD03,index_CpeD03,index_silD03,index_DH21Ad_P0.20_0.00_1.400,PAHion.31,PAHneu.31,q_D16graphite.dat,qlib_gra_D16MGemt_1.400} "$HERE/data/dielectric/"
cp "$SED"/data/{kext_astrodust_MW.dat,kext_astrodust_MW_euv.dat,kext_dl07_MW.dat,kext_dl07_MW_euv.dat,kext_zubko_BARE_GR_S.dat} "$HERE/data/"
cp "$SED"/data/release/{extinction.dat,kext_albedo_WD_MW_3.1_60_D03.all_2003,scattering.dat,size_distribution.dat} "$HERE/data/release/"
# The astrodust Q table, cut at the Lyman limit.  Its wavelength axis becomes
# the model's own grid (SEDust's euv_extended_lambda_grid), so the table's
# extent sets both the band the transport can be given cross sections over and
# how much work the emission solver does per cell.  Upstream's table runs to
# 1.0e-4 um (12.4 keV); the 633 wavelengths it adds below 0.0912 um carry no
# photons at MoCafe's par%lambda_min and cost about a factor 2.4 in
# dust-emission wall time, so they are cut here.  The cut is a pure truncation
# on whole a_eff blocks -- no value is recomputed -- which is why taking
# par%lambda_min below the Lyman limit needs nothing but the untruncated
# upstream file copied in over this one.
awk '/^#/ { if ($0 ~ /lambda stride/) sub(/of 1762/, "of 1129"); print; next }
     (NF >= 8 && $1+0 >= 0.0912) { print }' \
   "$SED"/tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat \
   > "$HERE/tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat"
# Zubko (ZDA 2004 BARE-GR-S) optics/calorimetry/config for par%dust_model='zubko'
cp "$SED"/data/zubko/* "$HERE/data/zubko/"
echo "Populated SEDust data under $HERE (from $SED)"
