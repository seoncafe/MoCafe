# Quasi-Random Photon Launching: Applicability to MoCafe v2.00

> **Status (2026-07-18): Stage Q1 implemented.**  `src/qmc_mod.f90` (ported
> from the MoCHII reference implementation; raw net verified against
> `scipy.stats.qmc.Sobol`), `par%launch_sequence = 'random'|'sobol'`,
> `par%qmc_seed`; inverse-CDF wavelength samplers next to the alias tables;
> fixed 7-dimension layout (origin, component, wavelength, direction,
> external-sphere entry) consumed in `gen_photon`; the global photon number
> is carried as `photon%id` (LaRT convention).  Validated: default path
> bit-identical; sobol `Direct0` bit-identical across MPI task counts;
> wavelength-bin deviations O(1) vs binomial (max 1.9 vs 139 at N=1e5);
> scattered field unbiased (0.7 sigma over 4+4 replicates).  Excluded (setup
> error): Stokes, `(a,g)`/tau scans, external rec/cyl,
> `radiation_angular_PDF_file`.  Dust-emission packets stay on the Mersenne
> Twister (step 6 below not taken).

This note analyzes whether the randomized quasi-Monte Carlo (RQMC) launch
scheme proposed for MoCHII (`MoCHII/docs/QUASI_RANDOM_LAUNCH.md`) is worth
adopting in MoCafe v2.00, and what would have to change.  The short answer:
**yes, as an opt-in variance-reduction feature, with a smaller expected
payoff than in MoCHII and two MoCafe-specific obstacles (alias-table
wavelength sampling and rejection-based source geometries) that shape the
design.**

## 1. What the MoCHII proposal does

MoCHII replaces the pseudo-random draws of the *launch* variables (source
component, frequency bin, direction, entry point) of each stellar packet by an
Owen-scrambled Sobol point indexed by the global photon number, while keeping
the Mersenne Twister for every conditional, variable-length draw after launch
(scattering, Russian roulette).  Key requirements: a fixed meaning for every
Sobol dimension, global (MPI-independent) photon indexing, a monotonic inverse
CDF for discrete choices, and a fixed scramble key across nonlinear
iterations with several independent scrambles for error estimation.

The scheme is especially attractive in MoCHII because its ionizing transport
is absorption-only with *analytic* attenuation: given the launch variables,
the direct contribution is deterministic, so RQMC addresses a fixed,
low-dimensional integral.

## 2. How MoCafe differs, and what that implies

| Aspect | MoCHII | MoCafe v2.00 | Consequence for RQMC |
|---|---|---|---|
| Post-launch transport | Analytic attenuation (deterministic) | Stochastic scattering random walk (forced first scattering, HG, peel-off) | Only the **direct** component is deterministic given the launch variables; the scattered field keeps its `N^(-1/2)` behavior |
| Direct estimator | Analytic rates | `peeling_direct_photon`, called once at birth in `gen_photon` | The `Direct`/`Direct0` images and the external direct peel inherit the full RQMC benefit unchanged |
| Frequency/wavelength choice | Monotonic band CDF (`ion_cdf`) | **Alias tables** (`sed_mod`: `sed_src_alias`; `sources_mod`: `src_alias`; external field: `sed_ext_alias`) | Alias sampling scrambles the unit interval; the wavelength dimension needs a monotonic inverse-CDF sampler first (the MoCHII note flags exactly this caveat for the SED path) |
| Source position | Point or simple surface mappings | Point, but also disks/bulges that use **rejection** (`exp_spiral`, `sersic` box rejection, `boxy`/`bar`/`xbar`) and stateful samplers (`rand_gauss`, `rand_r1exp`) | Rejection loops consume a variable number of uniforms and cannot occupy fixed Sobol dimensions; those geometries keep Mersenne Twister positions |
| Global photon index | `ip` strided loop | Present in **both** runners: equal-share (`run_simulation_mod.f90:134`) and master-slave (`ip = numreceived - num_send_at_once + ii`, line 66) | Random-access Sobol indexed by `ip` works under either scheduler and makes the launch set MPI-count independent |
| Nonlinear iteration | Ionization/thermal loop | Lucy dust-emission iteration (`dust_niter`) | The common-random-numbers argument carries over: a fixed scramble across Lucy iterations stabilizes the convergence differences |
| Regression gates | Fixed-state gamma tests | Bit-identity gates for the `(a,g)`/tau scans and the shared peel routines | RQMC must be strictly opt-in (`'random'` default); the scans should reject the QMC mode initially |

The structural conclusion: MoCafe cannot expect the near-analytic gains MoCHII
anticipates, because scattering keeps a stochastic history after launch.  What
remains is still substantial, and it concentrates exactly where MoCafe v2.00
runs spend their photon budget.

## 3. Where RQMC launch would help MoCafe (ranked)

1. **Wavelength-bin stratification in SED runs.**  With plain Monte Carlo the
   packet counts in the wavelength bins fluctuate binomially, and every
   wavelength-resolved product (`Direct`/`Scattered` cubes, `J_lambda`,
   the absorbed power that drives dust emission) inherits that noise bin by
   bin.  Inverting one stratified Sobol coordinate through the luminosity CDF
   makes the count in a bin of probability `p` deviate from `pN` by order one.
   This is the cheapest, most certain gain, and it applies to every SED run
   regardless of optical depth.
2. **`Direct` / `Direct0` images and the emergent direct SED.**  The direct
   peel is evaluated at birth from the launch variables alone, so it is the
   MoCafe analogue of MoCHII's analytic rates: launch-only RQMC transfers to
   it without any transport change.  This includes the external direct peels.
3. **External-field entry uniformity.**  The isotropic external illumination
   (`external_illumination_sph/cyl/rec`) uses closed-form mappings of 4-5
   uniforms; low-discrepancy points remove the surface clumping of entry
   positions and incidence angles that currently adds large-scale noise to
   external and composed runs.
4. **Origin and component allocation.**  The compose mode draws
   internal-vs-external per packet by the luminosity ratio, and the
   multi-source engine draws the emitting component by CDF; stratifying these
   discrete choices removes the shot noise on each subpopulation's photon
   budget (equivalent to near-deterministic splitting, with no weight
   bookkeeping).
5. **Lucy dust-iteration stability.**  Reusing one fixed scrambled launch set
   through the `dust_niter` loop correlates successive radiation fields, so
   the iteration-to-iteration change measures front/temperature motion rather
   than resampling noise.  Independent scrambles then give the replicate
   error estimate, as in the MoCHII policy.

## 4. Where it would not help

- The **scattered** images and SED at moderate-to-high optical depth: the
  variable-length random walk after launch dominates and stays pseudo-random.
- The Stokes path (Mueller-matrix scattering; alias-sampled phase angles) —
  out of scope, as in MoCHII.
- Sharp local features: clump shadow edges, single-pixel statistics.
- The `(a,g)`/tau scans: they are variance-reduction schemes of their own with
  bit-identity guarantees tied to the current RNG stream; combining them with
  RQMC needs a separate design and should be rejected at setup initially.

## 5. MoCafe-specific design

### 5.1 Dimension layout (superset, fixed for every packet)

| Dimension | Variable | Notes |
|---|---|---|
| 1 | Origin: internal vs external (compose mode) | Unused when not composing |
| 2 | Source component (`src_lum_cdf`) | Unused when `nsource = 1` |
| 3 | Wavelength bin | Inverse CDF, not the alias table |
| 4 | Polar direction `mu` (or cosine-weighted incidence) | |
| 5 | Azimuth `phi` | |
| 6-8 | Entry/position coordinates | External sph: 2; rec: 3 (face + 2 coordinates); cyl: 3; point: unused |

Unused-but-fixed dimensions are harmless; what must never happen is a
geometry-dependent reassignment of meanings.  Source geometries whose position
samplers are rejection-based or stateful (`gaussian`, `exponential`, `sech`,
`exp_spiral`, `sersic`, `boxy`, `bar`, `xbar`) keep Mersenne Twister draws for
the position only, while still taking wavelength and direction from the Sobol
point; this preserves fixed dimension semantics because the Mersenne Twister
draws live outside the Sobol stream.  (A later refinement could replace
`rand_zexp`/`rand_sech2` by their closed-form inverses, which exist, but
`rand_r1exp` and the bulge samplers are genuinely rejection-based.)

### 5.2 Inverse-CDF wavelength sampler

Keep the alias tables for the default path and add cumulative arrays
(`sed_src_cdf`, `src_cdf(:,is)`, `sed_ext_cdf`) built from the already
normalized probabilities, sampled by binary search (about `log2(nlambda) = 7`
steps; negligible).  Only the QMC path uses them, so the default sampling
stream is untouched.

### 5.3 Generator

As in the MoCHII note: a self-contained random-access Sobol generator
(Joe & Kuo direction numbers, at most 8 dimensions) with nested-uniform
scrambling by a keyed hash (Laine & Karras 2011; Burley 2020), indexed by the
global photon number and keyed by `(qmc_seed, dimension, replicate)`.  No
external dependency; a few tens of integer operations per packet.  Namelist:

```fortran
par%launch_sequence = 'random'   ! 'random' (default) or 'sobol'
par%qmc_seed        = 12345
```

The scramble key stays fixed across Lucy iterations within one run.
Dust-emission packets (`gen_dust_photon`) keep Mersenne Twister sampling
initially: their emitting-cell CDF changes between iterations and the
histories are scattering-dominated, mirroring the "diffuse emission later"
recommendation in the MoCHII note.

### 5.4 Reproducibility

Indexing by global `ip` makes the launch set independent of the MPI task
count and of the master-slave batching.  As in MoCHII, full bitwise
MPI-independence of a run is not achieved (the scattered histories still come
from rank-local generators); the claim is only about the launch set and the
direct-peel contributions.

## 6. Expected benefit, honestly stated

- For **scattering-only SED runs at low-to-moderate optical depth** (the
  configuration of `examples/sed_point`, `sed_external`, `sed_multi_source`),
  the direct component dominates most bins, so the emergent-SED noise should
  drop markedly; the scattered residual sets the floor.
- For **dust-emission runs**, the absorbed-power field is fed by direct and
  scattered pathlengths; only the direct share improves, but the wavelength
  and component stratification still improves the energy budget entering
  SEDust, and the fixed scramble stabilizes the Lucy iteration.
- For **optically thick, scattering-dominated** cases the gain will be small.
- No universal speedup factor should be assumed; report the effective
  photon-number gain (ordinary-MC packets needed to match the RQMC RMS) from
  replicate scatter, as the MoCHII note prescribes.

## 7. Validation plan (adapted to MoCafe)

- **V0 - net quality.**  Projections of the generated coordinates; CDF
  agreement of bin counts; isotropy; identical launch sets for different MPI
  task counts; distinct but balanced nets for different `qmc_seed`.
- **V1 - direct field.**  `examples/sed_point` at small `taumax`: compare
  ordinary MC vs several scrambles for the RMS of the `Direct0` cube and of
  bin counts vs `sed_lum`; the direct peel must be exactly the launch set.
- **V2 - unbiasedness of the scattered field.**  `taumax = 1` sphere: the
  scattered images from QMC-launch runs must agree with ordinary MC within
  replicate scatter (launch-only RQMC does not bias the transport estimator,
  but this must be demonstrated).
- **V3 - multi-source and compose.**  `sed_multi_source`, `compose_int_ext`:
  component/origin counts vs the luminosity ratios; total-luminosity checks
  (`TOT_LUM = L_int + L_ext`) must hold identically.
- **V4 - dust emission.**  One `examples/dustemis` benchmark: replicate RMS of
  the absorbed and emitted luminosity, PASS on the existing energy-conservation
  checks, and Lucy-iteration convergence traces with fixed vs varied scramble.
- **Gate.**  The default `'random'` path must remain bit-identical throughout;
  QMC + `(a,g)`/tau scan and QMC + Stokes are rejected at setup.

## 8. Verdict and recommended order

Adopting the MoCHII design in MoCafe is **worthwhile but second-order**: the
strongest MoCHII argument (deterministic analytic transport after launch)
does not carry over, so the payoff concentrates in (i) wavelength-bin
stratification of every SED product, (ii) the direct images/SED, (iii) the
external-field entry uniformity and origin/component allocation, and (iv)
Lucy-iteration stability.  The implementation cost is genuinely small (a
self-contained generator plus an inverse-CDF wavelength sampler, both outside
the default stream), and the design constraints (fixed dimensions, global
indexing, opt-in gating) are already worked out in the MoCHII note.

Suggested order:

1. Random-access scrambled Sobol generator + `par%launch_sequence` switch.
2. Single internal point source in SED mode: wavelength (inverse CDF), `mu`,
   `phi`; pass V0-V2.
3. Multi-source and compose layouts (dimensions 1-2); pass V3.
4. External entry mappings (dimensions 6-8); rerun V3 on `sed_external`.
5. Fixed-scramble policy through the Lucy iteration; pass V4.
6. Only then consider dust-emission packets or scattering-order dimensions,
   each with its own validation campaign.
