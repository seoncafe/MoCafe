#!/usr/bin/env python3
"""Standalone clumpy-dust population generator for MoCafe (dust-only).

Writes a clump file (FITS binary table or HDF5) that MoCafe reads via
``par%clump_input_file``.  The layout mirrors the Fortran ``write_clumps_info``
routine in ``src/clump_mod.f90``:

  * HDU 0 / group 0 : empty primary
  * HDU 2 / group 2 : binary table with columns ``X, Y, Z`` (clump centers, in
    code units) and, when the population is non-uniform, the per-clump columns
    ``R_CLUMP`` (radius) and ``RHOKAP`` (grey dust opacity per code length).
    Constant radius / opacity are stored in the header keywords ``CL_RAD`` /
    ``RHOKAP`` instead (MoCafe falls back to those when a column is absent).

Header keywords: ``N_CLUMPS, SPHERE_R, RMIN, CL_RAD, F_VOL, F_COV, TAU0,
RHOKAP, RMAX, DISTUNIT, DIST_CM``.

This is the dust-only slim of LaRT's ``make_clumps.x`` (the velocity /
temperature / Voigt columns and keywords are dropped); see
``AMR_CLUMPS_PLAN.md`` Part B / B4.

Examples
--------
Single giant clump filling the sphere (validation reference -- reproduces a
homogeneous dust sphere of radial optical depth tau0):

    python make_clumps.py --single --rmax 1.0 --clump-radius 1.0 \
        --tau0 1.0 --out clump_giant.fits

Uniform random non-overlapping population at a covering factor:

    python make_clumps.py --rmax 1.0 --clump-radius 0.1 --f-cov 1.0 \
        --tau0 1.0 --out clumps.fits
"""
import argparse
import sys
import numpy as np


# --------------------------------------------------------------------------
def opacity_from_inputs(args):
    """Return the grey dust opacity per code length (cl_rhokap) for one clump.

    Priority: tau0 (dimensionless, center->surface) > ndust > nH.  The ndust /
    nH branches follow MoCafe's convention rhokap = density * cext_dust *
    distance2cm, matching grid_create / clump_mod.
    """
    if args.tau0 is not None and args.tau0 > 0:
        return args.tau0 / args.clump_radius
    if args.ndust is not None and args.ndust > 0:
        return args.ndust * args.cext_dust * args.dist_cm
    if args.nh is not None and args.nh > 0:
        return args.nh * args.cext_dust * args.dist_cm
    raise SystemExit("ERROR: specify one of --tau0 / --ndust / --nh")


# --------------------------------------------------------------------------
def n_clumps_from_inputs(args):
    """Derive the clump count from --n-clumps / --f-vol / --f-cov."""
    R, r, rmin = args.rmax, args.clump_radius, args.rmin
    if args.n_clumps is not None and args.n_clumps > 0:
        return int(args.n_clumps)
    if args.f_vol is not None and args.f_vol > 0:
        return int(round(args.f_vol * (R**3 - rmin**3) / r**3))
    if args.f_cov is not None and args.f_cov > 0:
        return int(round((4.0 / 3.0) * args.f_cov *
                         (R**2 + R * rmin + rmin**2) / r**2))
    raise SystemExit("ERROR: specify one of --n-clumps / --f-vol / --f-cov")


# --------------------------------------------------------------------------
def place_rsa(n, R, r, rmin, fully_inside, rng):
    """Random Sequential Addition of n non-overlapping spheres of radius r in
    the shell [rmin, R].  Cell-hash accelerated; aborts on jamming."""
    if fully_inside:
        rmax_c, rmin_c = R - r, rmin + r
        if rmax_c <= rmin_c:
            raise SystemExit("ERROR: --fully-inside but no room "
                             "(rmin + 2*clump_radius > rmax).")
    else:
        rmax_c, rmin_c = R, rmin
    min_sep2 = (2.0 * r) ** 2

    # cell-hash grid (cell >= 2r so the 27-neighbor scan is complete)
    ng = max(1, int(2.0 * R / max(2.0 * r, 1e-30)))
    cell = 2.0 * R / ng
    from collections import defaultdict
    grid = defaultdict(list)

    def cidx(p):
        return tuple(np.clip(((p + R) / cell).astype(int), 0, ng - 1))

    pos = np.empty((n, 3))
    placed, stall, MAX_STALL = 0, 0, 5_000_000
    while placed < n:
        stall += 1
        if stall > MAX_STALL:
            raise SystemExit(
                f"ERROR: RSA jamming after 5e6 attempts at clump {placed+1}; "
                "lower the filling factor or clump_radius (f_vol < ~0.35).")
        p = (2.0 * rng.random(3) - 1.0) * rmax_c
        d2 = float(p @ p)
        if d2 > rmax_c**2 or d2 < rmin_c**2:
            continue
        ci, cj, ck = cidx(p)
        ok = True
        for ii in range(max(0, ci-1), min(ng, ci+2)):
            for jj in range(max(0, cj-1), min(ng, cj+2)):
                for kk in range(max(0, ck-1), min(ng, ck+2)):
                    for q in grid[(ii, jj, kk)]:
                        dq = pos[q] - p
                        if dq @ dq < min_sep2:
                            ok = False
                            break
                    if not ok:
                        break
                if not ok:
                    break
            if not ok:
                break
        if not ok:
            continue
        pos[placed] = p
        grid[(ci, cj, ck)].append(placed)
        placed += 1
        stall = 0
    return pos


# --------------------------------------------------------------------------
def write_clumps(fname, xyz, radius, rhokap, args, fmt):
    """Write the clump table.  Constant radius/opacity -> header keywords only;
    otherwise per-clump R_CLUMP / RHOKAP columns are written too."""
    n = xyz.shape[0]
    R, rmin = args.rmax, args.rmin
    f_vol = float(np.sum(radius**3) / max(R**3 - rmin**3, 1e-30))
    f_cov = float(0.75 * np.sum(radius**2) /
                  max(R**2 + R * rmin + rmin**2, 1e-30))
    const_tol = 1e-3
    write_radius = radius.ptp() > const_tol * max(abs(radius.mean()), 1e-30)
    write_rhokap = rhokap.ptp() > const_tol * max(abs(rhokap.mean()), 1e-30)

    keys = dict(
        N_CLUMPS=n, SPHERE_R=R, RMIN=rmin, CL_RAD=float(radius.max()),
        F_VOL=f_vol, F_COV=f_cov,
        TAU0=(args.tau0 if args.tau0 is not None else -1.0),
        RHOKAP=float(rhokap[0]), RMAX=R,
        DISTUNIT=args.distunit, DIST_CM=args.dist_cm)

    if fmt == "fits":
        from astropy.io import fits
        cols = [fits.Column(name="X", format="1E", array=xyz[:, 0].astype("f4")),
                fits.Column(name="Y", format="1E", array=xyz[:, 1].astype("f4")),
                fits.Column(name="Z", format="1E", array=xyz[:, 2].astype("f4"))]
        if write_radius:
            cols.append(fits.Column(name="R_CLUMP", format="1E",
                                    array=radius.astype("f4")))
        if write_rhokap:
            cols.append(fits.Column(name="RHOKAP", format="1E",
                                    array=rhokap.astype("f4")))
        tab = fits.BinTableHDU.from_columns(cols)
        for k, v in keys.items():
            tab.header[k] = v
        fits.HDUList([fits.PrimaryHDU(), tab]).writeto(fname, overwrite=True)
    elif fmt == "hdf5":
        import h5py
        # Group "2" mirrors the FITS HDU-2 binary table: one dataset per column,
        # keywords as group attributes (matches MoCafe's HDF5 io facade).
        with h5py.File(fname, "w") as f:
            g = f.create_group("CLUMPS")
            g.create_dataset("X", data=xyz[:, 0].astype("f4"))
            g.create_dataset("Y", data=xyz[:, 1].astype("f4"))
            g.create_dataset("Z", data=xyz[:, 2].astype("f4"))
            if write_radius:
                g.create_dataset("R_CLUMP", data=radius.astype("f4"))
            if write_rhokap:
                g.create_dataset("RHOKAP", data=rhokap.astype("f4"))
            for k, v in keys.items():
                g.attrs[k] = v
    else:
        raise SystemExit(f"ERROR: unknown format {fmt}")

    print(f"make_clumps: wrote {n} clumps to {fname} "
          f"(f_vol={f_vol:.4g}, f_cov={f_cov:.4g}, rhokap={rhokap[0]:.4g}, "
          f"R_clump={radius.max():.4g})")


# --------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description="Generate a MoCafe clump file.")
    ap.add_argument("--out", required=True, help="output file (.fits/.fits.gz/.h5)")
    ap.add_argument("--rmax", type=float, required=True, help="sphere radius [code units]")
    ap.add_argument("--clump-radius", type=float, required=True, help="clump radius [code units]")
    ap.add_argument("--rmin", type=float, default=0.0, help="inner cavity radius")
    ap.add_argument("--single", action="store_true",
                    help="one clump centered at the origin (validation reference)")
    ap.add_argument("--n-clumps", type=int, default=None)
    ap.add_argument("--f-vol", type=float, default=None)
    ap.add_argument("--f-cov", type=float, default=None)
    ap.add_argument("--tau0", type=float, default=None,
                    help="dust optical depth of one clump (center->surface)")
    ap.add_argument("--ndust", type=float, default=None, help="clump dust density")
    ap.add_argument("--nh", type=float, default=None, help="clump H density [cm^-3]")
    ap.add_argument("--cext-dust", type=float, default=1.6059e-21,
                    help="dust extinction per H [cm^2] (matches par%%cext_dust)")
    ap.add_argument("--dist-cm", type=float, default=1.0, help="cm per code unit")
    ap.add_argument("--distunit", default="", help="distance unit string")
    ap.add_argument("--no-fully-inside", action="store_true",
                    help="allow clumps to protrude past the sphere boundary")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--format", choices=["fits", "hdf5", "auto"], default="auto")
    args = ap.parse_args(argv)

    fmt = args.format
    if fmt == "auto":
        fmt = "hdf5" if args.out.endswith((".h5", ".hdf5")) else "fits"

    kap = opacity_from_inputs(args)
    rng = np.random.default_rng(args.seed)

    if args.single:
        xyz = np.zeros((1, 3))
    else:
        n = n_clumps_from_inputs(args)
        if n < 1:
            raise SystemExit("ERROR: derived N_clumps < 1")
        xyz = place_rsa(n, args.rmax, args.clump_radius, args.rmin,
                        not args.no_fully_inside, rng)

    radius = np.full(xyz.shape[0], args.clump_radius)
    rhokap = np.full(xyz.shape[0], kap)
    write_clumps(args.out, xyz, radius, rhokap, args, fmt)


if __name__ == "__main__":
    main()
