#!/usr/bin/env python3
"""Build a generic AMR octree dust sphere for MoCafe (dust-only validation).

Writes a FITS/HDF5 binary table that MoCafe reads via ``par%amr_file``
(grid_type='amr'), with the schema expected by read_generic_amr.f90:

  columns : x, y, z, level, nH, T, vx, vy, vz   (mandatory, by index)
  header  : NAXIS2 (= nleaf), BOXLEN, ORIGINX/Y/Z

The medium is a uniform-density sphere of radius = box half-size inscribed in
the box; cells with center outside the sphere get nH = 0.  The octree can be
uniform (every leaf at --level) or radially refined (finer toward the center,
which exercises the octree's level-crossing ray traversal while keeping the
physical medium a uniform sphere -- so it must still reproduce the homogeneous
Cartesian sphere).  T and velocities are written for format compatibility but
are ignored by the dust-only transport.

This is a compact stand-in for the full python/AMR_grid/AMR_grid.py builder;
it covers the validation sphere.  See AMR_CLUMPS_PLAN.md Part A / F4.

Examples
--------
    python make_amr_sphere.py --level 5 --out amr_sphere_uniform.fits
    python make_amr_sphere.py --level-min 3 --level-max 6 --out amr_sphere_refined.fits
"""
import argparse
import numpy as np


def build_leaves(half, level_min, level_max, r_sphere):
    """Recursively subdivide the root cube [-half,half]^3 into an octree.
    A cell at `level` is split until it reaches the target level for its
    center radius (finer toward the center).  Returns arrays of leaf centers
    and levels.  Octant convention matches the Fortran amr_build_tree:
    child center = parent center + (2*bit-1)*child_half."""
    cx, cy, cz, lv = [], [], [], []

    def target_level(r):
        if level_max == level_min:
            return level_min
        # finer toward the center: full level_max for r < half/2, ramping to
        # level_min at the box corner.
        t = min(1.0, r / (0.75 * half))
        return int(round(level_max - t * (level_max - level_min)))

    def rec(x, y, z, h, level):
        r = np.sqrt(x*x + y*y + z*z)
        if level < target_level(r):
            hc = 0.5 * h
            for b in range(8):
                ix = b & 1
                iy = (b >> 1) & 1
                iz = (b >> 2) & 1
                rec(x + (2*ix - 1)*hc, y + (2*iy - 1)*hc, z + (2*iz - 1)*hc,
                    hc, level + 1)
        else:
            cx.append(x); cy.append(y); cz.append(z); lv.append(level)

    rec(0.0, 0.0, 0.0, half, 0)
    return (np.array(cx), np.array(cy), np.array(cz), np.array(lv, dtype=np.int32))


def main(argv=None):
    ap = argparse.ArgumentParser(description="Generate a generic AMR dust sphere.")
    ap.add_argument("--out", required=True, help="output file (.fits/.fits.gz/.h5)")
    ap.add_argument("--boxlen", type=float, default=2.0, help="box side length [code units]")
    ap.add_argument("--level", type=int, default=None, help="uniform refinement level")
    ap.add_argument("--level-min", type=int, default=3)
    ap.add_argument("--level-max", type=int, default=6)
    ap.add_argument("--n0", type=float, default=1.0, help="uniform density inside the sphere")
    ap.add_argument("--temperature", type=float, default=1.0e4)
    args = ap.parse_args(argv)

    half = 0.5 * args.boxlen
    r_sphere = half
    if args.level is not None:
        lmin = lmax = args.level
    else:
        lmin, lmax = args.level_min, args.level_max

    x, y, z, lv = build_leaves(half, lmin, lmax, r_sphere)
    n = x.size
    r = np.sqrt(x*x + y*y + z*z)
    nH = np.where(r < r_sphere, args.n0, 0.0)
    T  = np.full(n, args.temperature)
    vz = np.zeros(n);  vx = np.zeros(n);  vy = np.zeros(n)

    keys = dict(BOXLEN=args.boxlen, ORIGINX=-half, ORIGINY=-half, ORIGINZ=-half)

    if args.out.endswith((".h5", ".hdf5")):
        import h5py
        with h5py.File(args.out, "w") as f:
            g = f.create_group("AMR_GRID")
            for nm, arr in (("x", x), ("y", y), ("z", z), ("level", lv.astype("i4")),
                            ("nH", nH), ("T", T), ("vx", vx), ("vy", vy), ("vz", vz)):
                g.create_dataset(nm, data=arr)
            for k, v in keys.items():
                g.attrs[k] = v
            g.attrs["NLEAF"] = n
    else:
        from astropy.io import fits
        cols = [
            fits.Column(name="x",     format="1D", array=x),
            fits.Column(name="y",     format="1D", array=y),
            fits.Column(name="z",     format="1D", array=z),
            fits.Column(name="level", format="1J", array=lv.astype("i4")),
            fits.Column(name="nH",    format="1D", array=nH),
            fits.Column(name="T",     format="1D", array=T),
            fits.Column(name="vx",    format="1D", array=vx),
            fits.Column(name="vy",    format="1D", array=vy),
            fits.Column(name="vz",    format="1D", array=vz),
        ]
        tab = fits.BinTableHDU.from_columns(cols)
        for k, v in keys.items():
            tab.header[k] = v
        fits.HDUList([fits.PrimaryHDU(), tab]).writeto(args.out, overwrite=True)

    print(f"make_amr_sphere: wrote {n} leaves to {args.out} "
          f"(levels {lv.min()}-{lv.max()}, boxlen={args.boxlen}, "
          f"f_inside={np.mean(nH>0):.3f})")


if __name__ == "__main__":
    main()
