# MoCafe User Guide

MoCafe is a Monte-Carlo dust radiative-transfer code (Fortran 90 + MPI).  This
guide covers building, running, every input parameter, the single-run scans,
the clumpy and AMR media, the RAMSES/Illustris-TNG converters, and the Python
analysis tools.  All physics conventions (rotation angles, image geometry,
output normalization) are in `docs/MoCafe_Geometry.pdf` and
`docs/output_definitions.pdf`.

## How to Compile and Run

### Prerequisites

1. **Fortran/MPI compiler** — Fortran 2008 (Intel oneAPI `mpiifort`/`mpiifx`, or
   GNU `mpif90`).  MoCafe is MPI-only.
2. **CFITSIO** — <https://heasarc.gsfc.nasa.gov/fitsio/> (required).
3. **HDF5** (optional, default on) — 1.14+ with Fortran bindings built by the
   same compiler.  Build with `make HDF5=1` (default) or `make HDF5=0` to skip.
4. **Python 3** (for the analysis tools / converters): `numpy`, `astropy`,
   `h5py`, and `scipy` (Illustris converter only).

### Build

```bash
make                                # HDF5=1 (default) -> MoCafe.x
make HDF5=0                         # CFITSIO only (no HDF5 link)
make HDF5_PREFIX=/usr/local/hdf5    # override HDF5 location
make F90=mpif90                     # force a compiler
make DEBUG=1                        # checked/traceback flags
```

The compiler is auto-selected `mpiifort → mpiifx → mpif90`.  The build always
passes `-cpp -DMPI` and (default) `-DHDF5`.  When adding a `.f90` file, append
it to `OBJSB` in the `Makefile` in dependency order.

### Run

```bash
mpirun -np <N> ./MoCafe.x <input>.in
```

The single argument is a Fortran namelist file.  Sample inputs and `run.sh`
scripts are under `examples/*/`.

### Reference

See `params_type` in `src/define.f90` for the default value of every parameter.
All parameters take the `par%` prefix in the namelist (`&parameters / par%... /`).

---

## Input Parameters

### Simulation control

| Parameter | Default | Description |
|-----------|---------|-------------|
| `no_photons` | 1e5 | Number of photon packets |
| `no_print` | 1e7 | Progress print interval (photons) |
| `iseed` | 0 | RNG seed (0 → seeded from time/pid) |
| `launch_sequence` | `'random'` | `'sobol'` = quasi-random (Owen-scrambled Sobol) launch of the emission direction and the external entry point (`sph`, `rec`, `cyl`), indexed by the global photon number; transport draws stay on the Mersenne Twister.  Not combined with `use_stokes`, the `(a,g)`/tau scans, or `radiation_angular_PDF_file`.  See `docs/QUASI_RANDOM_LAUNCH_MOCAFE.md` |
| `qmc_seed` | 12345 | Scramble seed for `'sobol'`; a different seed is an independent replicate |
| `luminosity` | 1.0 | Source luminosity (output is per unit luminosity) |
| `use_master_slave` | `.true.` | Master–slave MPI (vs. equal-share when `.false.`) |
| `num_send_at_once` | 10000 | Photon batch size in master–slave mode |
| `use_reduced_wgt` | `.true.` | Forced first scattering + reduced-weight immortal photons |

### Grid geometry (Cartesian)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `nx`, `ny`, `nz` | 1, 1, 11 | Cells per axis |
| `nr` | -999 | Shorthand: sets `nx=ny=nz=nr` |
| `xmax`, `ymax`, `zmax` | 1.0 | Half-extent per axis; box is `(-xmax, xmax)` etc. |
| `rmax` | -999 | Sphere radius; if `>0` and not a cylinder, the box becomes a cube of side `rmax` and density is 0 for `r > rmax` |
| `rmin` | 0.0 | Inner cavity radius (shell / clump cavity) |
| `geometry` | `''` | `'sphere'`, `'cylinder'`, `'Plummer'` (else a uniform/file box) |
| `xyz_symmetry` | `.false.` | Use 1/8 of the box (auto-disabled if any of nx/ny/nz = 1, and not for peel-off) |
| `xy_periodic` | `.false.` | Periodic slab in x,y |
| `z_symmetry` | `.false.` | Mirror at z = 0 |

The system center is always at `(0, 0, 0)`.

### Grid type

| Parameter | Default | Description |
|-----------|---------|-------------|
| `grid_type` | `'car'` | `'car'` (Cartesian), `'clump'` (clumpy medium), `'amr'` (octree) |
| `use_clump_medium` | `.false.` | Legacy switch, maps to `grid_type='clump'` |
| `use_amr_grid` | `.false.` | Legacy switch, maps to `grid_type='amr'` |

### Distance unit

| Parameter | Default | Description |
|-----------|---------|-------------|
| `distance_unit` | `'kpc'` | `'kpc'`, `'pc'`, `'au'`, or `''` (dimensionless code units) |
| `distance2cm` | (from unit) | cm per code unit |

For a dimensionless model with no density file, `distance_unit` is forced to
`''` and the dust column is set by `taumax` / `tauhomo` (see below).

### Source geometry

| Parameter | Default | Description |
|-----------|---------|-------------|
| `source_geometry` | `'point'` | `'point'`, `'uniform'`, `'uniform_xy'`, `'gaussian'`, `'exponential'`, or `'external_{sph,cyl,rec}'` (and `external_*` variants) |
| `xs_point`, `ys_point`, `zs_point` | 0.0 | Point-source location |
| `source_zscale` | NaN | Scale height for `'gaussian'`/`'exponential'` sources |
| `radiation_angular_PDF_file` | `''` | Angular PDF for external illumination |

The scalar `source_geometry` above describes a **single** internal source or an
external-only field.  Two further configurations are available.

#### Multiple internal sources

| Parameter | Default | Description |
|-----------|---------|-------------|
| `nsource` | 1 | Number of internal sources (`> 1` activates the multi-source path) |
| `src_geometry(i)` | `'point'` | Geometry of each source: `'point'`, `'uniform'`, `'uniform_xy'`, `'gaussian'`, `'exponential'` |
| `src_lum(i)` | -999 | Luminosity of each source (equal split of `luminosity` if unset) |
| `src_x/y/z(i)` | 0.0 | Position of each `'point'` source |
| `src_zscale(i)` | -999 | Scale height of a `'gaussian'`/`'exponential'` source |

Sources are drawn in proportion to their luminosity and `luminosity` becomes
their sum.  Incompatible with `use_stokes` (the multi-source generator carries
no Stokes vector).

#### Composition: internal source(s) + external field

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ext_geometry` | `''` | Boundary the external field enters through when it **composes** with internal sources: `'sph'`, `'cyl'`, `'rec'` (blank auto-selects from `geometry`/`rmax`) |
| `ext_intensity` | -999 | Band mean intensity `J` of the isotropic external field; `> 0` switches the composition on |

Composition activates when internal sources (`nsource >= 1`) coexist with
`ext_intensity > 0` and `source_geometry` is *not* one of the `external_*`
external-only shorthands.  Each photon is then drawn internal-vs-external in
proportion to

```
L_tot = L_int + L_ext = Σ src_lum(i) + π · J · A_surface,
```

with `A_surface` the illuminated boundary area (`sph`: `4π(rmax·distance2cm)²`;
`rec`: the six box faces; `cyl`: the barrel plus the two end caps), and the
output is normalized to `L_tot`.  Because energy is additive, a composed run
equals the internal-only run plus the external-only run.  Incompatible with
`use_stokes` and with the `(a,g)`/tau scans.

### Dust & scattering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `albedo` | 0.3253 | Dust single-scattering albedo |
| `hgg` | 0.6761 | Henyey–Greenstein asymmetry `g` |
| `cext_dust` | 1.6059e-21 | Dust extinction per H nucleus [cm²] |
| `scatt_mat_file` | `''` | Tabulated Mueller matrix (`data/mueller_*.dat`); enables Stokes |
| `use_stokes` | `.true.` | Stokes polarization (auto-off if no `scatt_mat_file`) |

### Optical depth & density

| Parameter | Default | Description |
|-----------|---------|-------------|
| `taumax` | -999 | Radial optical depth, center → outer edge (sets the density scale) |
| `tauhomo` | -999 | Optical depth of the volume-equivalent homogeneous medium |
| `density_file` | `''` | 3-D density cube (FITS/HDF5); `rhokap = density·cext_dust·distance2cm` |
| `reduce_factor` | 1 | Down-sampling factor for the density file |
| `centering` | 0 | Re-centering mode for the density file |
| `density_zscale`, `density_rscale`, `density_powerlaw_index` | NaN | Analytic density modulations for `'sphere'`/`'cylinder'`/`'Plummer'` |

### Peel-off observers

If `nxim`,`nyim` > 0, an image is peeled off toward each observer.  Up to
`MAX_OBSERVERS = 181` observers per run.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `nobs` | 1 | Number of observers |
| `nxim`, `nyim` | 129 | Image pixels |
| `dxim`, `dyim` | NaN | Pixel size (auto to cover the grid if unset) |
| `distance` | NaN | Observer distance; renormalizes `(obsx,obsy,obsz)` |
| `inclination_angle`, `position_angle`, `phase_angle` | NaN | Observer angles [deg] (arrays) |
| `alpha`, `beta`, `gamma` | NaN | Euler rotation angles [deg] (alternative) |
| `obsx`, `obsy`, `obsz` | NaN | Explicit observer direction (alternative) |

See `docs/MoCafe_Geometry.pdf` for the angle conventions.

### Output control

| Parameter | Default | Description |
|-----------|---------|-------------|
| `file_format` | `'hdf5'` | `'hdf5'` or `'fits'` — authoritative for the output container |
| `out_file` | `''` | Output base name (derived from the input name if empty) |
| `out_bitpix` | -32 | FITS bitpix (-32 = float32, -64 = float64) |
| `save_direc0` | `.false.` | Also write the unattenuated direct image `Direct0` |
| `output_normalization` | `'luminosity'` | Output normalization |
| `sightline_tau` | `.false.` | Also write a per-pixel sight-line optical-depth map |

---

## Single-Run Scans

One Monte-Carlo run yields the scattered image for a whole grid of parameters.
Both scans are restricted to the no-Stokes Henyey–Greenstein path
(`use_stokes = .false.`).  References:
[Seon 2010, PKAS, 25, 177](https://ui.adsabs.harvard.edu/abs/2010PKAS...25..177S/abstract)
and [Jonsson 2006, MNRAS, 372, 2](https://ui.adsabs.harvard.edu/abs/2006MNRAS.372....2J/abstract);
the algorithm memo is `docs/MoCafe_agtau_scan.pdf`.

### (albedo, asymmetry-factor) scan — Seon 2010

| Parameter | Default | Description |
|-----------|---------|-------------|
| `use_ag_list` | `.false.` | Enable the `(a,g)` scan → `scatt(x,y,albedo,hgg)` |
| `albedo_list` | NaN | Albedo grid (empty → 0.1 … 1.0) |
| `hgg_list` | NaN | Asymmetry grid (empty → 0.0 … 0.9) |

Photons are tracked at the single simulated `g0 = par%hgg` (use 0.4–0.5);
per-`(a,g)` reweighting fills the extra axes with no extra random numbers.  The
`Direct`/`Direct0` images stay 2-D.  The slice at `(a*, g0)` is bit-identical to
a single run.

### Polychromatic optical-depth (tau) scan — Jonsson 2006

| Parameter | Default | Description |
|-----------|---------|-------------|
| `use_tau_list` | `.false.` | Enable the `tau` scan → adds a `tau` axis |
| `tau_list` | NaN | Target `taumax` values (empty → √2-spaced `tau/taumax` in [0.5,2]; the reference `par%taumax` is auto-inserted) |

Photons are tracked at the reference `tau0 = par%taumax`; a closed-form
reweighting (Jonsson eq. 31) produces every `tau_list` slice.  Composes with the
`(a,g)` scan:

| | image |
|---|---|
| neither | `scatt(x,y)` |
| `use_ag_list` | `scatt(x,y,a,g)` |
| `use_tau_list` | `scatt(x,y,tau)` + `direc(x,y,tau)` |
| both | `scatt(x,y,a,g,tau)` + `direc(x,y,tau)` |

Unlike `(a,g)`, the **direct** beam is `tau`-dependent (it gains a `tau` axis;
`Direct0` stays 2-D).  Variance grows as `tau/taumax` departs from 1 — keep
`tau_list` within ~2–3× of `par%taumax` and put `par%taumax` near the middle.
Works with internal point/extended sources **and** external illumination
(`external_{sph,cyl,rec}`) — the external direct beam is tau-aware.

---

## Clumpy Medium

`par%grid_type = 'clump'` (or `use_clump_medium = .true.`).  A sphere of radius
`par%rmax` holds N discrete, non-overlapping spherical dust clumps; the
inter-clump space is vacuum.  Dust only (no velocity/temperature).  Algorithm
memo: `docs/MoCafe_clump.pdf`; examples in `examples/clump_sphere/`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `rmax` | -999 | Outer sphere radius (required) |
| `clump_radius` | -1 | Clump radius (required) |
| `clump_f_cov` / `clump_f_vol` / `clump_N_clumps` | -1 | Clump count — specify exactly one |
| `clump_tau0` | -1 | Dust optical depth of one clump, center → surface (`κ = tau0/clump_radius`) |
| `clump_ndust` / `clump_nH` | -1 | Per-clump dust / H density (`κ = n·cext_dust·distance2cm`) |
| `rmin` | 0.0 | Inner cavity radius |
| `clump_fully_inside` | `.true.` | Require each clump entirely inside the shell |
| `cone_opening` | 0.0 | Biconical placement half-angle [deg] (0 → full sphere) |
| `clump_input_file` | `''` | Load a clump population instead of generating one |
| `save_clump_info` | `.false.` | Write the generated clump table to `<base>_clumps.*` |

If no per-clump opacity is given, `par%taumax` or `par%tauhomo` back-solves it.
RSA placement is reliable for volume filling factor ≲ 0.35; MoCafe aborts with
guidance if it cannot place the requested clumps.

**Spatially-varying profiles** (geometry only): `clump_radius_profile`,
`clump_density_profile`, `clump_number_profile` ∈
`constant | powerlaw | gaussian | exponential | file`, each with `*_alpha`,
`*_r0`; tabulated input via `clump_profile_file`.

**Generate a clump file** with the Python tool:

```bash
python python/make_clumps.py --single --rmax 1.0 --clump-radius 1.0 \
       --tau0 1.0 --out clump.fits          # one centered giant clump
python python/make_clumps.py --rmax 1.0 --clump-radius 0.1 --f-cov 1.0 \
       --tau0 1.0 --out clumps.fits          # uniform random population
```

---

## AMR Mode

`par%grid_type = 'amr'` (or `use_amr_grid = .true.`).  An adaptive octree is
read from a generic AMR file produced by the converters or the builder.  Dust
only.  Algorithm memo: `docs/MoCafe_amr.pdf`; examples in `examples/amr_sphere/`,
`examples/amr_tng/`, `examples/amr_ramses/`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `amr_type` | `'generic'` | Only `'generic'`; `'ramses'` is rejected (convert first) |
| `amr_file` | `''` | Generic AMR file (`.h5`, `.hdf5`, `.fits`, `.fits.gz`, `.dat`) |
| `dust_model` | `'global_dgr'` | Per-leaf dust opacity model (below) |
| `DGR` | 1e-2 | Dust-to-gas ratio for `global_dgr` |
| `Z_global` | 0.0134 | Uniform metallicity for `laursen09` when no `metallicity` column |
| `Z_ref` | 0.0134 | Reference solar metallicity |
| `f_ion_dust` | 0.01 | Dust survival fraction in ionized gas (Laursen+09) |

The dust scale is set by normalizing to `par%taumax` (radial pole sightline) or
`par%tauhomo` (volume average); for an asymmetric model prefer `tauhomo`.

### Generic AMR file format

A binary table (FITS/HDF5) or text file with **mandatory** columns (read by
name; an index fallback handles legacy files):

| Column | Unit | Description |
|--------|------|-------------|
| x, y, z | distance_unit | Leaf-cell centers |
| level | — | Octree refinement level |
| nH | cm⁻³ | Hydrogen number density (alias `dens`/`gasDen`) |

**Optional** columns (by name): `metallicity` (Z), `xHI` (neutral fraction),
`ndust` (dust pseudo-density).  Header keywords: `BOXLEN`, `ORIGINX/Y/Z`,
`NAXIS2` (= number of leaves).  Any `T`/`vx`/`vy`/`vz` columns are ignored.

### AMR dust models

| `dust_model` | Needs (columns) | Per-leaf opacity |
|--------------|-----------------|------------------|
| `global_dgr` | `nH` | `nH·cext_dust·DGR·distance2cm` |
| `from_file` | `ndust` | `ndust·cext_dust·distance2cm` |
| `laursen09` | `metallicity` **and** `xHI` | `(Z/Z_ref)(nHI + f_ion·nHII)·cext_dust·distance2cm` |

`laursen09` requires an explicit `xHI` column (the `T→xHI` CIE path is deferred
with temperature).

---

## Converters and Tools (`python/AMR_grid/`)

All converters are **dust only**: gas velocity and temperature are not read or
written.  See `AMR_CONVERTERS_PLAN.md`.

### Octree builder / validation grid

| Tool | Description |
|------|-------------|
| `AMR_grid.py` | In-memory `AMRGrid` octree builder: uniform / density-gradient refinement, `set_density/metallicity/neutral_fraction/dust_density`, FITS/HDF5/text writers (module-level `write_leaves` for flat leaf arrays). |
| `make_amr_sphere.py` | Generic AMR dust sphere (uniform level or radially refined) for validation. |

### RAMSES converter

`convert_ramses_to_generic.py` parses a native RAMSES output
(`amr_*`/`hydro_*` Fortran-unformatted files, directly — no `yt`).

```bash
python python/AMR_grid/convert_ramses_to_generic.py output_00042 \
       -o sim.fits --output-unit kpc --metallicity-index 6 --emit-ndust
```

| Flag | Default | Description |
|------|---------|-------------|
| `repository` | (positional) | RAMSES dir or `.../output_00042` |
| `-o, --output` | (required) | Output file (`.fits`/`.fits.gz`/`.h5`/`.dat`) |
| `-s, --snapnum` | inferred | Snapshot number |
| `--output-unit` | `kpc` | `cm`/`kpc`/`pc`/`au` |
| `--density-index` | from descriptor | 1-based hydro var for density |
| `--metallicity-index` | from descriptor | 1-based hydro var for Z (or `--metallicity Z` uniform) |
| `--ionization-index` | — | 1-based hydro var for `xHI` |
| `--emit-ndust` | off | Precompute `ndust` via Laursen+09 (`--Z-ref`, `--f-ion`) |

### Illustris-TNG converter

`convert_illustris_to_generic.py` reads a local TNG gas cutout (HDF5,
`PartType0`) or downloads one via the TNG.org API, converts comoving+h units to
physical kpc and `nH`, builds an octree, and assigns `nH`, `metallicity`
(`GFM_Metallicity`), `xHI` (`NeutralHydrogenAbundance`) by Voronoi
nearest-neighbor.

```bash
# local cutout
python python/AMR_grid/convert_illustris_to_generic.py cutout.hdf5 \
       -o galaxy.h5 --grid-type amr --level-max 7 --boxsize 60 --emit-ndust
# or download + convert
python python/AMR_grid/convert_illustris_to_generic.py \
       --api-key YOUR_KEY --simulation TNG50-1 --snap 99 --subhalo-id 0 \
       -o galaxy.h5
```

| Flag | Default | Description |
|------|---------|-------------|
| `input_file` | (positional/optional) | Local cutout HDF5 (omit when using `--api-key`) |
| `-o, --output` | (required) | Output file |
| `--grid-type` | `amr` | `amr` (adaptive) or `uniform` |
| `--level-min` / `--level-max` | 3 / 7 | Octree levels |
| `--dens-threshold` | 0.3 | Density-gradient refinement criterion |
| `--match-resolution` | off | Refine to the local Voronoi cell size (`--resolution-factor`) |
| `--center X Y Z` | density peak | Recenter [physical kpc] |
| `--boxsize` | enclose all | Octree box side [kpc] |
| `--metallicity Z` | — | Uniform Z when `GFM_Metallicity` absent |
| `--emit-ndust` | off | Precompute `ndust` via Laursen+09 |
| `--sfr-mask {keep,drop}` | keep | Density-based handling of star-forming cells (no temperature) |
| `--api-key`, `--simulation`, `--snap`, `--subhalo-id` | — | TNG API download |

### Python analysis (`python/`)

| Tool | Description |
|------|-------------|
| `mocafe_io.py` | Format-agnostic reader/converter for MoCafe outputs. CLI: `python mocafe_io.py {info,peek,convert} <file>`. API: `load_mocafe('out_obs.h5')` → `MoCafeFile` with `.scattered/.direct/.direct0/.total/.stokes('I')`, `.scattered_at(a=,g=,tau=)`, `.direct_at(tau=)`, axis accessors, and `find_mocafe_outputs(stem)`. |

```python
import mocafe_io
mf  = mocafe_io.load_mocafe('out_obs.h5')
img = mf.scattered_at(a=0.5, g=0.3, tau=1.5)   # nearest scan slice
```

---

## Output Formats

`par%file_format` is authoritative.  If `par%out_file`'s extension disagrees,
MoCafe rewrites it (and warns) so the main file and all derived files
(peel-off Stokes split, sight-line tau, per-observer images) stay in one format.

```
file_format='hdf5'  →  _obs.h5,        _obs_tau.h5,        _stokes.h5
file_format='fits'  →  _obs.fits.gz,   _obs_tau.fits.gz,   _stokes.fits.gz
```

HDF5 mirrors the FITS HDU layout (each image HDU → a group named after
`EXTNAME` with a `data` dataset; FITS keywords → group attributes), in the same
section order.  Convert between them with:

```bash
python python/mocafe_io.py convert out_obs.h5 out_obs.fits.gz
```

---

## Examples

| Directory | Shows |
|-----------|-------|
| `examples/point_source/`, `uniform_source/`, `external_*/` | Cartesian grids, point/extended/external sources |
| `examples/agtau_list/` | `(a,g)` and `tau` single-run scans (+ a notebook) |
| `examples/clump_sphere/` | Clumpy medium (uniform / profile / cone / file-loaded) + validation |
| `examples/amr_sphere/` | AMR octree sphere vs. the Cartesian reference |
| `examples/amr_tng/` | Illustris-TNG cutout → converter → MoCafe AMR |
| `examples/amr_ramses/` | RAMSES → converter → MoCafe AMR (template) |

Each directory's `run.sh` shows the typical invocation
(`mpirun -np $(nproc) ../../MoCafe.x <input>.in`).
