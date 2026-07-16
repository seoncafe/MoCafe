# MoCafe v1.20

MoCafe is a Monte-Carlo radiative-transfer code (Fortran 90 + MPI) for **dust
scattering**.  It follows photon packets emitted by a source through a dusty
medium, scatters them off dust grains (Henyey–Greenstein or a tabulated
Mueller matrix), and "peels off" a contribution toward each observer at every
scattering site to build a FITS or HDF5 image of the source seen through the
dust.  Author: Kwang-il Seon (KASI).

The output image is in units of `(luminosity unit) cm^-2 sr^-1`.

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

The clumpy and AMR media are **dust only** at this stage (no gas velocity,
temperature, or Lyman-α machinery).

## References

- [Seon & Draine 2016, ApJ, 833, 201](https://ui.adsabs.harvard.edu/abs/2016ApJ...833..201S/abstract) & [Seon & Kim 2020, ApJS, 250, 9](https://ui.adsabs.harvard.edu/abs/2020ApJS..250....9S/abstract) — dust/Monte-Carlo radiative transfer by the author (general background; the shared MC machinery is ported from LaRT).
- [Seon 2010, PKAS, 25, 177](https://koreascience.kr/journal/CMHHCI/v25n4.do) — the `(a,g)` single-run scan.
- [Jonsson 2006, MNRAS, 372, 2](https://ui.adsabs.harvard.edu/abs/2006MNRAS.372....2J/abstract) — the polychromatic optical-depth scan (SUNRISE).
- Amanatides & Woo 1987, *A Fast Voxel Traversal Algorithm for Ray Tracing*, Eurographics '87

## Documentation

- **[README_HOWTO.md](README_HOWTO.md)** — build, run, every input parameter,
  the scans, clumpy/AMR modes, the converters, and the Python tools.
- `docs/MoCafe_Geometry.pdf` — observer geometry and image conventions.
- `docs/output_definitions.pdf` — output image definitions.
- `docs/MoCafe_agtau_scan.pdf`, `docs/MoCafe_clump.pdf`, `docs/MoCafe_amr.pdf`
  — algorithm memos for the scans, the clumpy medium, and the AMR grid.
- `UPGRADE_NOTES.md` — change history; `CLAUDE.md` — architecture overview.

## Author

Kwang-il Seon (KASI/UST)

---

Last updated: 2026-07-17 08:35 KST
