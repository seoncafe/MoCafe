# Quasi-Random Photon Launching: Applicability to MoCafe v1.20

> **Status (2026-07-18): implemented** (backported from the validated v2.00
> Stage Q1).  `src/qmc_mod.f90`, `par%launch_sequence = 'random'|'sobol'`,
> `par%qmc_seed`, `photon%id` (LaRT convention); dimensions 4-5 (direction)
> and 6-7 (external-sphere entry) consumed, 1-3 reserved for v2.00-layout
> compatibility.  Validated: default path bit-identical (point, external_sph,
> tau scan; equal-share runners); sobol point-source `Direct0` bit-identical
> across MPI task counts; scattered totals unbiased.  Excluded (setup error):
> Stokes, `(a,g)`/tau scans, external rec/cyl, `radiation_angular_PDF_file`.
> Note: the external direct image rides the `peeling_direct_external_sph2`
> estimator, which draws its own pseudo-random line of sight (by design, not
> a defect), so the external `Direct0` stays MPI-count dependent; its total
> is conserved.

This note analyzes whether the randomized quasi-Monte Carlo (RQMC) launch
scheme proposed for MoCHII (`MoCHII/docs/QUASI_RANDOM_LAUNCH.md`) is worth
adopting in MoCafe v1.20, the monochromatic scattering code.  A companion
analysis for v2.00 lives at
`MoCafe_v2.00/docs/QUASI_RANDOM_LAUNCH_MOCAFE.md`.  The short answer:
**the case is markedly weaker than for v2.00 — the two strongest v2.00 gains
(wavelength-bin stratification and dust-iteration stability) do not exist in
a monochromatic code — but extended-source and externally illuminated runs
would still benefit through the direct image and the first-scattering
distribution.  Recommendation: implement and validate in v2.00 first;
backport to v1.20 only if a production need appears.**

## 1. The MoCHII proposal in one paragraph

Replace the pseudo-random draws of the *launch* variables of each packet
(source component, frequency, direction, entry point) by an Owen-scrambled
Sobol point indexed by the global photon number; keep the Mersenne Twister
for every conditional, variable-length draw after launch.  Requirements:
fixed meaning for every Sobol dimension, global (MPI-independent) photon
indexing, monotonic inverse CDFs for discrete choices, and independent
scrambles for error estimation.  MoCHII expects large gains because its
ionizing transport is absorption-only with analytic attenuation, so the
entire post-launch estimator is deterministic.

## 2. What v1.20 offers RQMC — and what it does not

v1.20 transports a **single wavelength**: there is no frequency dimension, no
multi-source engine, no external-field composition, and no dust-emission
iteration.  The launch space is therefore small:

| Source type | Launch variables (closed-form?) | Dimensions |
|---|---|---|
| `point` | direction `mu`, `phi` (yes) | 2 |
| `uniform` sphere | radius `r = rmax * u^(1/3)`, `cos(theta)`, `phi` (yes) + direction (yes) | 5 |
| `uniform_xy`, `exponential`, `gaussian` | position (closed-form; `gaussian` uses the stateful Box-Muller sampler) + direction | 4-5 |
| `external_sph` / `external_rec` / `external_cyl` | entry point + cosine-weighted incidence (yes) | 4-5 |

Both runners expose a global photon index `ip` (equal-share:
`run_simulation_mod.f90:126`; master-slave: `ip = numreceived -
num_send_at_once + ii`, line 65), so a random-access Sobol point indexed by
`ip` is feasible under either scheduler, exactly as in MoCHII.

Two v1.20-specific observations shape the verdict.

**(a) For a point source, the direct image gains nothing.**  The direct
peel-off is evaluated at birth from the (fixed) source position with unit
weight; it does not depend on the launch direction at all.  The monochromatic
point-source `Direct`/`Direct0` images are already noise-free, so RQMC can
only act *indirectly*, by distributing the forced-first-scattering sites more
uniformly.  Since the first scattering order often dominates the scattered
image at the optical depths v1.20 is used for, this is a real but partial
effect: later orders stay pseudo-random.

**(b) Extended and external sources are the genuine v1.20 use case.**  For a
`uniform`/disk source the pixel-to-pixel noise of the direct surface
brightness comes from position sampling, and for external illumination it
comes from the entry-point and incidence sampling; both are fixed-dimension,
closed-form mappings that inherit the full RQMC benefit without touching the
transport.  Externally illuminated clouds (the `external_sph` family used for
the LDN 1642 work) are the configuration where launch stratification would be
most visible.

## 3. Interaction with the v1.20 production modes

- **`(a,g)` and tau scans.**  The scans reuse one forward history for the
  whole scan grid, so launch stratification would reduce noise coherently
  across every `(a, g, tau)` slice at once — structurally compatible, since
  the scans consume no launch draws.  However, the scans' bit-identity gate
  (a scan slice must reproduce a stand-alone run bit for bit) must then be
  re-validated *within* the QMC mode, and the default `'random'` path must
  stay untouched.  Initially the QMC mode should be gated off when
  `use_ag_list`/`use_tau_list` is set, as in the v2.00 plan.
- **Stokes / Mueller path.**  Phase-angle sampling uses a continuous
  alias-linear table (`scattering_car.f90:124`) — post-launch, so it simply
  stays pseudo-random.  Launch-only RQMC is compatible with Stokes transport
  (the birth Stokes basis follows deterministically from the direction), but
  the scattered polarization gains nothing.
- **ISRF angular distribution.**  The optional external angular distribution
  (`radiation_angular_PDF_file`) is sampled by `rand_alias_linear`
  (`external_radiation.f90:103`).  Alias sampling breaks the monotonicity a
  stratified coordinate needs; a QMC path would replace it with interpolated
  inverse-CDF sampling for that dimension (cheap, table is small).

## 4. Expected benefit, honestly stated

| Configuration | Expected gain |
|---|---|
| Point source, monochromatic | Small: only the first-scattering distribution improves; direct image already exact |
| Extended source (uniform/disk) | Moderate: smoother direct surface brightness + first-scattering stratification |
| External illumination | Moderate: uniform entry coverage removes large-scale launch noise from both direct and scattered fields |
| Deep (scattering-dominated) clouds | Small: the multiple-scattering random walk dominates and stays pseudo-random |
| `(a,g)`/tau scan production runs | Same as the underlying configuration, applied coherently to all slices — but gated initially |

There is no wavelength dimension to stratify and no nonlinear iteration to
stabilize, which removed the two highest-value items of the v2.00 analysis.
The generator cost itself is negligible (a self-contained random-access
Sobol + hash scrambling, at most 5 dimensions here), and the namelist
interface (`par%launch_sequence = 'random'|'sobol'`, `par%qmc_seed`) would
mirror the MoCHII/v2.00 design with the default `'random'` guaranteed
bit-identical.

## 5. Verdict and recommendation

RQMC launching in v1.20 is **technically straightforward** (small fixed
launch dimension, global photon index available, closed-form mappings for
the geometries that matter) but **strategically secondary**: v1.20 is the
maintenance line, its point-source runs barely benefit, and the strongest
arguments for RQMC live in v2.00 (wavelength stratification, multi-source
and compose allocation, Lucy-iteration stability).

Recommended course:

1. Implement and validate the feature in v2.00 first, following
   `MoCafe_v2.00/docs/QUASI_RANDOM_LAUNCH_MOCAFE.md`.
2. Backport to v1.20 only if an extended-source or externally illuminated
   production campaign (an LDN 1642-type analysis) shows a photon-budget
   bottleneck in the launch-dominated components.  The backport is then
   small: the direction/position/entry mappings above, the inverse-CDF
   replacement for the ISRF angular table, and the same validation gates
   (net quality, direct-field RMS, scattered-field unbiasedness,
   MPI-count-independent launch sets) restricted to the monochromatic case.
3. Keep the scans and the Stokes path excluded until the plain scattering
   mode has passed its gates; then re-validate the scan slice-identity
   property within the QMC mode before allowing the combination.
