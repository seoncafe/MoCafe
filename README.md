# MoCafe v2.00

MoCafe is a Monte-Carlo radiative-transfer code (Fortran 90 + MPI) for **dust
scattering and thermal emission**.  It follows photon packets emitted by a
source through a dusty medium, scatters them off dust grains (Henyey–Greenstein
or a tabulated Mueller matrix), and "peels off" a contribution toward each
observer to build a FITS or HDF5 image.  **v2.00** adds panchromatic dust
thermal emission: multi-wavelength transport, a radiation-field tally in each cell,
and UV-to-FIR SEDs of dusty sources, external galaxies, and the Milky Way.

The scattered-light image is in units of `(luminosity unit) cm^-2 sr^-1`.

## Features

- **Three media**, selected by `par%grid_type`:
  - `'car'` — uniform/structured **Cartesian** grid (point/extended/external
    sources; sphere, cylinder, Plummer, power-law, exponential geometries; or a
    3-D density file).
  - `'clump'` — **clumpy dust**: N non-overlapping spherical dust clumps in a
    sphere (RSA placement, CSR/DDA ray-tracing), generated internally or loaded
    from a file.
  - `'amr'` — **adaptive octree** read from a generic AMR file, fed by RAMSES /
    Illustris-TNG converters.
- **Single-run scans** (one Monte-Carlo run → many images):
  - `(albedo, asymmetry)` scan → `scatt(x,y,a,g)` (Seon 2010).
  - polychromatic `optical-depth` scan → adds a `tau` axis (Jonsson 2006);
    composes with the `(a,g)` scan into `scatt(x,y,a,g,tau)`.
- **Stokes polarization** via a Mueller-matrix phase function.
- **HDF5 and FITS** output/input through a format-agnostic facade (`par%file_format`).
- **MPI** parallelism (master–slave or equal-share), with MPI-3 shared memory
  for read-only grid/clump/octree arrays (one copy per node).
- **Dust thermal emission (v2.00)** — panchromatic (SED) transport, a `J_λ` tally in
  each cell, and two emission methods:
  - **Lucy (1999) + [SEDust](SEDust/)** — equilibrium + stochastically heated
    grains + PAHs; iterable for dust self-absorption; fast equilibrium
    (`dust_single_teq`) and table (`dust_fast_table`) options.
  - **Bjorkman & Wood (2001)** — immediate reemission, no iteration.
  - Plus **multiple stellar populations**, a **HEALPix all-sky interior
    observer** (the Milky-Way case), and a **Modified Random Walk** for high τ.
  - Runs on **both the Cartesian grid and the AMR octree**.
  - See **[README_HOWTO.md](README_HOWTO.md#dust-thermal-emission-v200)** and
    `MoCafe_v2.00_PLAN.md`.  Validated against the Camps et al. (2015) SHG
    dust-emission benchmark.

The clumpy and AMR media carry **dust only** (no gas velocity, temperature, or
Lyman-α radiative transfer).  The dust-emission modes run on the Cartesian grid
and the AMR octree (not the clumpy medium).

## References

- [Seon & Draine 2016, ApJ, 833, 201](https://ui.adsabs.harvard.edu/abs/2016ApJ...833..201S/abstract) & [Seon & Kim 2020, ApJS, 250, 9](https://ui.adsabs.harvard.edu/abs/2020ApJS..250....9S/abstract) — dust/Monte-Carlo radiative transfer by the author (general background; the shared Monte-Carlo core is ported from LaRT).
- [Seon 2010, PKAS, 25, 177](https://koreascience.kr/journal/CMHHCI/v25n4.do) — the `(a,g)` single-run scan.
- [Jonsson 2006, MNRAS, 372, 2](https://ui.adsabs.harvard.edu/abs/2006MNRAS.372....2J/abstract) — the polychromatic optical-depth scan (SUNRISE).
- [Lucy 1999, A&A, 344, 282](https://ui.adsabs.harvard.edu/abs/1999A%26A...344..282L/abstract) — the iterative dust-emission / mean-intensity estimator.
- [Bjorkman & Wood 2001, ApJ, 554, 615](https://ui.adsabs.harvard.edu/abs/2001ApJ...554..615B/abstract) — immediate reemission with temperature correction.
- [Camps et al. 2015, A&A, 580, A87](https://ui.adsabs.harvard.edu/abs/2015A%26A...580A..87C/abstract) — the stochastic-heating dust-emission benchmark.
- [Robitaille 2010, A&A, 520, A70](https://ui.adsabs.harvard.edu/abs/2010A%26A...520A..70R/abstract) & [Min et al. 2009, A&A, 497, 155](https://ui.adsabs.harvard.edu/abs/2009A%26A...497..155M/abstract) — the Modified Random Walk.
- Amanatides & Woo 1987, *A Fast Voxel Traversal Algorithm for Ray Tracing*, Eurographics '87

## Documentation

- **[README_HOWTO.md](README_HOWTO.md)** — build, run, every input parameter,
  the scans, clumpy/AMR modes, the converters, and the Python tools.
- **`docs/MoCafe_v2.00_UserGuide.pdf`** — the full v2.00 user guide (install,
  run, complete input reference, dust-emission workflow, outputs, examples,
  algorithms; LaTeX source `docs/MoCafe_v2.00_UserGuide.tex`).
- `docs/MoCafe_Geometry.pdf` — observer geometry and image conventions.
- `docs/MoCafe_agtau_scan.pdf`, `docs/MoCafe_clump.pdf`, `docs/MoCafe_amr.pdf`
  — algorithm memos for the scans, the clumpy medium, and the AMR grid.

## Author

Kwang-il Seon (KASI).

---

Last updated: 2026-07-11 11:51 KST
