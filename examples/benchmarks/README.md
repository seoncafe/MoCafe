# MoCafe v2.00 dust-emission benchmarks (Stage 5)

## 1. SHG dust-emission benchmark (Camps et al. 2015, A&A 580, A87)

Validates the **emission engine** (the SEDust library) against the
published multi-code stochastic-heating-of-grains (SHG) benchmark, using the
same Zubko et al. (2004) BARE-GR-S dust model the benchmark defines.

- `shg_bench.f90` links `SEDust/sed/lib/libsedust.a`, builds the Zubko
  BARE-GR-S model through `build_dust` from `SEDust/data`, and computes the dust
  emission for the Mathis ISRF scaled by U = 1e-2 … 1e6.  Build/run (the archive
  reads its optics as HDF5, so it comes before the HDF5 libraries):
  ```
  H=/data/opt/hdf5_intel
  ifort -O2 -qopenmp -I../../SEDust/sed/lib shg_bench.f90 \
        ../../SEDust/sed/lib/libsedust.a \
        $H/lib/libhdf5_fortran.a $H/lib/libhdf5.a -lsz -ldl -lz -lm \
        -o shg_bench.x
  OMP_NUM_THREADS=8 ./shg_bench.x
  ```
  It builds on the whole 1201-point axis (`include_euv = .true.`), the range the
  benchmark's own reference spectra cover; the non-ionizing view would cut it at
  the Lyman limit, which is right for a transport run but not for this.
- **What it measures is the library as delivered.**  Whether SEDust's `zubko`
  is made of the Zubko et al. tables as distributed or of SEDust's own
  recomputation of them is SEDust's decision, and reproducing the published
  model is SEDust's business; this program goes through the plain entry point
  and reports how far the result is from the seven-code median.  It matters
  which: the two differ by up to 36% over the bins above a tenth of the peak
  (measured here), because the recomputation reproduces the silicate to 2.4e-6
  but the graphite only to 8% in `Q_abs` and implements the PAH component
  differently.  SEDust now makes the choice an argument, `zubko_optics`, whose
  default `'zda'` is the distributed tables the seven codes read; this program
  leaves it unset, so what it measures is the library's own default, and a
  MoCafe run with `par%dust_model = 'zubko'` gets the same model.
- `cmp_shg.py` / `plot_shg.py` compare the output against the reference
  spectra of 7 codes (CRT, DIRTY, SKIRT, DustEM, TRADING, MCFOST, DARTRAY)
  under `~/MoCafe/Grain/SHG_Benchmark/Results_FullSolution/`.  Those files are
  in the benchmark's own units, ~7 decades below ours at every U, so both
  scripts normalize before comparing: the benchmark constrains the **shape**,
  while the level is set by the absorbed power and is checked by the energy
  conservation test below.  `plot_shg.py` used to overlay them raw, which piled
  every envelope into one clump at the bottom of the axes instead of putting it
  on the curve it belongs to; the figure has been regenerated.

**Result:** the SEDust emission agrees with the 7-code median to a **median
relative difference of 2.5–4.4%** across U = 1e-2…1e4 (3.8 / 4.2 / 4.4 / 2.5% at
U = 1e-2 / 1e0 / 1e2 / 1e4), inside the multi-code envelope at every U, and the
emission peak wavelength matches the code median to within one wavelength bin
(288/308, 142/143, 60/62, 34/34 µm) — i.e. within the inter-code scatter.
See `shg_benchmark.pdf`.

Measured against SEDust `ff4ff47` plus the `zubko_optics` revision on top of it.
The agreement had gone to 1.6–5.9% for as long as the library served the
recomputation for `zubko`; making the distributed tables the default moved the
spectra back by 1.3–2.7% (median over the bins above a tenth of the peak, ~36%
at worst) and the agreement to the figures above.  Before that, `04b6ed8`'s
`calc_P` fix — which stopped the emission solvers absorbing photons the
radiation field does not carry — had moved the spectra by 0.4% (median over the
bins above a thousandth of the peak, 3–6% at worst) and the agreement by at
most 0.3 percentage points.

## 2. Internal mode cross-check (Mode 1 Lucy vs Mode 2 B&W)

`make_bench.sh` generates a τ_V sweep (0.1, 1, 5, 20) for a point source in a
uniform astrodust sphere, run in both modes.  The absorbed (= emitted) energy
fraction agrees between the equilibrium single-T Lucy solve and the Bjorkman
& Wood immediate-reemission mode:

| τ_V | Lucy (single-Teq, niter=1) | B&W 2001 |
|-----|---------------------------|----------|
| 0.1 | 0.0686 | 0.0686 |
| 1.0 | 0.454  | 0.455  |
| 5.0 | 0.861  | 0.864  |
| 20  | 0.989  | 0.9998 |

At τ_V = 20 the non-iterative Lucy (niter=1) misses dust self-absorption,
which B&W captures through immediate reemission; running Lucy with
`par%dust_niter > 1` closes the gap (the two agree to <0.3% at τ_V ≤ 5).

These fractions were not remeasured for SEDust `04b6ed8` and do not need to be.
What they report is the absorbed energy, which the transport fixes: at
`niter = 1` nothing is reemitted before it is counted, and B&W reemits without
the SEDust solver at all (`bw_mod` does not call it).  The update leaves the
absorbed luminosity of every `examples/dustemis/model_compare_*.in` run
unchanged to every digit; it moves how that energy is distributed over
wavelength, which is section 1 above.

## External RT benchmarks (deferred — need reference data)

The 3-D transport benchmarks in the plan (DUSTY 1-D, Pascucci 2-D disk,
Gordon et al. 2017 TRUST 3-D slab, and a direct Hyperion cross-check) require
their published reference solutions, which are not included here.  Add them
under this directory when available.
