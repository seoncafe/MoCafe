#!/bin/bash
# Refresh this SEDust copy's data from a canonical SEDust tree.
#
# What a MoCafe run needs is already in the repository, so a fresh checkout runs
# dust emission out of the box -- you normally do NOT need this script.  Run it
# to bring this copy up to date with a newer SEDust, or to fetch the files kept
# out of git: the 17 MB D16 spheroid Q-library (read only when
# qpah_graphite_source is 'd16_spheroid', a sensitivity setting no MoCafe path
# selects), and the text products described below.
#
# WHAT THE LAYOUT IS.  SEDust keeps one directory per dust model, holding
# everything that model owns, and keeps what the models SHARE beside those
# directories:
#
#   data/<model>/     one per par%dust_model: astrodust, dl07, mrn, zubko
#   data/dielectric/  optical constants, PAH cross sections
#   data/release/     published tables + the size distribution
#
# WHAT IS COMMITTED, AND WHAT IS NOT.  A model's whole optics -- its wavelength
# axis, its (lambda, a_eff) cross sections and its size-integrated extinction
# curve -- is ONE file, data/<model>/sedust_<model>.h5, and that is what
# build_dust reads and what this tree ships.  SEDust writes the same numbers as
# text as well (the q_*.dat Q tables and the kext_*.dat curves, ~250 MB for the
# four models), for a tree built without HDF5; those are regenerable by the
# programs that wrote the product, nothing here opens them, and they are not
# committed.  This script copies them anyway when the source tree has them, so
# that an HDF5=0 build has what it needs.
#
# The distributed ZDA optics tables under zubko/ are out of git for the same
# reason: the product carries them as its default /qtable, so a run that has the
# product opens none of them.
#
# What is NOT optics still comes from beside the product and IS committed: the
# size distribution under release/, the dielectric functions and PAH cross
# sections under dielectric/, and the ZDA config, size distributions and
# calorimetry under zubko/.
#
# The sources under sed/src/ and sed/build_lib.sh are byte-identical to their
# SEDust v1.00 counterparts, so `cmp` against that tree is the whole check that
# this copy is current.  The one file of the library that reaches into the
# T-matrix, src/euv_astrodust_tmatrix.f90, is deliberately not copied: it serves
# the astrodust EUV band, which only a lam_min shorter than 0.0912 um asks for,
# and MoCafe never passes one -- its ionizing band, when it wants one, is inside
# the stored axis and comes off the product like everything else.  Without that
# file the archive links with no T-matrix at all, and euv_tmatrix = .true. is
# refused (status 8 out of build_dust) instead of being answered with a
# different grain.  The
# T-matrix driver's own inputs (DH21_*, q_DH21Ad_*.gz under astrodust/) are not
# copied for the same reason.
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
mkdir -p "$HERE/data"

# Shared material data.  A dielectric function belongs to no one model -- DL07
# and Zubko read the same D03 astrosilicate -- so it sits beside the model
# directories rather than inside one of them.
mkdir -p "$HERE/data/dielectric" "$HERE/data/release"
cp "$SED"/data/dielectric/{index_CpaD03,index_CpeD03,index_silD03,index_DH21Ad_P0.20_0.00_1.400,PAHion.31,PAHneu.31,q_D16graphite.dat,qlib_gra_D16MGemt_1.400} "$HERE/data/dielectric/"
cp "$SED"/data/release/{extinction.dat,kext_albedo_WD_MW_3.1_60_D03.all_2003,kext_albedo_MRN,scattering.dat,size_distribution.dat} "$HERE/data/release/"

# One directory per model: the product, the text products behind it, and -- for
# zubko, whose model definition IS a set of files -- that definition.
for m in astrodust dl07 mrn zubko; do
  mkdir -p "$HERE/data/$m"
  cp "$SED"/data/"$m"/sedust_"$m".h5 "$HERE/data/$m/"
  for f in "$SED"/data/"$m"/q_*.dat "$SED"/data/"$m"/kext_*.dat; do
    [ -e "$f" ] && cp "$f" "$HERE/data/$m/"
  done
done
# The distributed ZDA optics tables.  They ARE the published BARE-GR-S model,
# and sedust_zubko.h5 now carries them as its default /qtable, so a run with the
# product opens none of them -- measured: with all three absent, the output of
# examples/dustemis/model_compare_zubko.in is identical dataset by dataset.
# They are therefore not committed, and are fetched here for a tree built
# without HDF5, which has no other source for this model.
cp "$SED"/data/zubko/{Gra_121_1201.dat,suvSil_121_1201.dat,PAH_28_1201_neu.dat} "$HERE/data/zubko/"
cp "$SED"/data/zubko/{Graphitic_Calorimetry_1000.dat,Silicate_Calorimetry_1000.dat} "$HERE/data/zubko/"
cp "$SED"/data/zubko/{ZDA_BARE_GR_S_Config.dat,ZDA_BARE_GR_S_SzDist_Gra.dat,ZDA_BARE_GR_S_SzDist_PAH.dat,ZDA_BARE_GR_S_SzDist_Sil.dat,zubko_descriptor.txt} "$HERE/data/zubko/"

echo "Populated SEDust data under $HERE (from $SED)"
