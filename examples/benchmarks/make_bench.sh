#!/bin/bash
# Generate tau-sweep benchmark inputs for Mode 1 (single-Teq) and Mode 2 (B&W).
mk() {  # mk <tau> <method> <extra>
  tau=$1; method=$2; out="bench_tau${tau}_${method}"
  cat > ${out}.in <<EOF
&parameters
 par%no_photons      = 2.0e6
 par%use_sed         = .true.
 par%use_dustemis    = .true.
 par%dust_emission_method = '${method}'
 par%dust_single_teq = .true.
 par%save_jlam       = .true.
 par%luminosity      = 3.828e33
 par%nlambda         = 120
 par%lambda_min      = 0.0912
 par%lambda_max      = 3000.0
 par%lambda_ref      = 0.55
 par%kext_file       = '../../data/kext_astrodust_MW.dat'
 par%tstar           = 1.0e4
 par%taumax          = ${tau}
 par%source_geometry = 'point'
 par%distance_unit   = 'pc'
 par%rmax = 1.0
 par%nx = 17
 par%ny = 17
 par%nz = 17
 par%nxim = 33
 par%nyim = 33
 par%obsx = 0.0
 par%obsy = 0.0
 par%obsz = 1.0
 par%distance = 1.0e5
 par%use_master_slave = .false.
 par%iseed = 1234
 par%no_print = 2.0e6
 par%file_format = 'hdf5'
 par%out_file = '${out}.h5'
/
EOF
}
for tau in 0.1 1.0 5.0 20.0; do
  mk $tau lucy
  # bw variant: drop single_teq line (bw ignores it) and lucy jlam requirement
  sed 's/dust_emission_method = .lucy./dust_emission_method = '"'"'bw01'"'"'/' bench_tau${tau}_lucy.in | \
    sed '/dust_single_teq/d; /save_jlam/d; /sed_qtable/d; /sed_sizedist/d' > bench_tau${tau}_bw01.in
  sed -i "s/bench_tau${tau}_lucy.h5/bench_tau${tau}_bw01.h5/" bench_tau${tau}_bw01.in
done
echo "generated bench inputs"
ls bench_*.in
