# Illustris-TNG dust-scattering example (dust-only AMR)

End-to-end demonstration of the TNG → MoCafe pipeline:

```
TNG gas cutout (HDF5, PartType0)
   │  python/AMR_grid/convert_illustris_to_generic.py   (dust-only)
   ▼
generic AMR octree file (x,y,z,level,nH[,metallicity,xHI,ndust])
   │  MoCafe   par%grid_type='amr'
   ▼
scattered-light image
```

Velocity and temperature are **not** used (dust only): the converter reads the
`Density` → `nH`, `GFM_Metallicity` → `metallicity`, and
`NeutralHydrogenAbundance` → `xHI`, and (with `--emit-ndust`) precomputes
`ndust` via Laursen+09.  MoCafe's dust opacity then comes from
`par%dust_model`:
- `global_dgr` — `nH·cext·DGR` (needs only `nH`),
- `laursen09`  — `(Z/Z_ref)(nHI+f_ion·nHII)·cext` (uses `metallicity`+`xHI`),
- `from_file`  — `ndust·cext` (needs `--emit-ndust`).

## Getting a cutout

Either supply a local TNG gas cutout HDF5, or download one (needs a free
tng-project.org API key):

```
python ../../python/AMR_grid/convert_illustris_to_generic.py \
    --api-key YOUR_KEY --simulation TNG50-1 --snap 99 --subhalo-id 448830 \
    -o tng_uniform.fits --grid-type uniform --level-max 5 --boxsize 60
```

## Run

```
./run.sh /path/to/cutout_TNG50-1_snap99_sub448830.hdf5
```

`run.sh` converts the cutout to a uniform-level-5 grid (`tng_uniform.fits`) and
an adaptive level-3–6 octree (`tng_amr.h5`), then runs `tng.in` and
`tng_amr.in`.  The dust is normalized to a volume-averaged optical depth
`par%tauhomo=1` (robust for an asymmetric galaxy; the pole sightline is not
representative).  `tng_montage.png` shows the scattered light tracing the
galaxy's dust.

## Converter options (see `--help`)

`--grid-type {amr,uniform}`, `--level-min/--level-max`, `--dens-threshold`,
`--match-resolution`, `--center`, `--boxsize`, `--metallicity Z` (uniform
fallback), `--emit-ndust`, `--sfr-mask {keep,drop}` (density-based, no
temperature), and the API-download flags.
