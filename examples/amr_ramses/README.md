# RAMSES dust-scattering example (dust-only AMR)

Pipeline:  RAMSES output → generic AMR octree → MoCafe.

```
python ../../python/AMR_grid/convert_ramses_to_generic.py \
    /path/to/output_00042 -o ramses.fits --output-unit kpc \
    --metallicity-index 6 --emit-ndust
mpirun -np $(nproc --all) ../../MoCafe.x ramses.in
```

`convert_ramses_to_generic.py` parses the native RAMSES octree
(`amr_*`/`hydro_*` Fortran-unformatted files, directly — no `yt`) and writes,
per leaf, `x,y,z,level,nH` plus optional `metallicity` (from the metal hydro
variable, `--metallicity-index`, or a uniform `--metallicity Z`), `xHI`
(`--ionization-index`), and `ndust` (`--emit-ndust`).  **Velocity and
temperature are not read** (the velocity/energy/pressure hydro variables are
skipped) — this is the dust-only converter.

The hydro variable indices come from `hydro_file_descriptor.txt` when present
(density and metallicity are detected by name); otherwise pass
`--density-index` / `--metallicity-index` explicitly.

MoCafe's dust opacity is chosen by `par%dust_model` in `ramses.in`:
`global_dgr` (nH only), `laursen09` (metallicity + xHI), or `from_file`
(`--emit-ndust`).  See `ramses.in` for a template and `docs/MoCafe_amr.pdf`
for the algorithm.

> Note: this directory ships only the template and converter; provide your own
> RAMSES `output_NNNNN/` directory.  The converter's Fortran-record reader is a
> faithful port of LaRT's (validated on real RAMSES); its generic-file *output*
> path is validated against MoCafe in this tree.
