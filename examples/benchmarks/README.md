# MoCafe v2.00 dust-emission benchmarks (Stage 5)

## 1. SHG dust-emission benchmark (Camps et al. 2015, A&A 580, A87)

Validates the **emission engine** (the SEDust library) against the
published multi-code stochastic-heating-of-grains (SHG) benchmark, using the
same Zubko et al. (2004) BARE-GR-S dust model the benchmark defines.

- `shg_bench.f90` links `SEDust/sed/lib/libsedust.a`, builds the Zubko
  BARE-GR-S model (`build_zubko`) from the optics committed under
  `SEDust/data/zubko/`, and computes the dust emission for the Mathis ISRF
  scaled by U = 1e-2 … 1e6.  Build/run:
  ```
  ifort -O2 -qopenmp -I../../SEDust/sed/lib shg_bench.f90 \
        ../../SEDust/sed/lib/libsedust.a -o shg_bench.x
  OMP_NUM_THREADS=8 ./shg_bench.x
  ```
  It prints one warning about not finding `../data/kext_zubko_BARE_GR_S.dat`,
  because that path is relative to `SEDust/sed/` and the benchmark runs from
  here.  Nothing is missing: the benchmark exercises `dust_emission` only, and
  never asks the model for its extinction curve.
- `cmp_shg.py` / `plot_shg.py` compare the output against the reference
  spectra of 7 codes (CRT, DIRTY, SKIRT, DustEM, TRADING, MCFOST, DARTRAY)
  under `~/MoCafe/Grain/SHG_Benchmark/Results_FullSolution/`.

**Result:** the SEDust emission agrees with the 7-code median to a **median
relative difference of 2.5–4.4%** across U = 1e-2…1e4, and the emission peak
wavelength matches the code median to within one wavelength bin
(288/308, 142/143, 60/62, 34/34 µm) — i.e. within the inter-code scatter.
See `shg_benchmark.pdf`.

Remeasured against SEDust revision `04b6ed8`, whose `calc_P` fix stopped the
emission solvers absorbing photons the radiation field does not carry.  That
moved these spectra by 0.4% (median over the bins above a thousandth of the
peak, 3–6% at worst) and the agreement with the code median by at most 0.3
percentage points, from the 3.7 / 4.2 / 4.2 / 2.2% measured before it to
3.8 / 4.2 / 4.4 / 2.5% at U = 1e-2 / 1e0 / 1e2 / 1e4.

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
