#!/usr/bin/env python3
"""Convert an Illustris/IllustrisTNG snapshot or cutout to a generic AMR dust
file for MoCafe (dust-only).

TNG uses an unstructured Voronoi mesh; MoCafe needs an octree.  This converter
(1) reads PartType0 gas cells (local HDF5 cutout or TNG.org API download),
(2) converts comoving+h code units to physical kpc and nH [cm^-3],
(3) builds an octree via AMR_grid.AMRGrid (uniform level, or adaptive on the
density gradient, optionally matched to the local Voronoi cell size),
(4) assigns each leaf nH / metallicity / xHI by nearest-neighbor lookup
(exact for a Voronoi tessellation), optionally precomputes ndust (Laursen+09),
and (5) writes the generic AMR file MoCafe reads with par%grid_type='amr'.

DUST-ONLY: gas velocity and temperature are NOT read or written (only
x,y,z,level,nH,metallicity [,xHI,ndust]).  See AMR_CONVERTERS_PLAN.md.

Usage
-----
    python convert_illustris_to_generic.py cutout.hdf5 -o galaxy.h5 \
        --grid-type amr --level-max 7 --boxsize 60 --emit-ndust
    python convert_illustris_to_generic.py --api-key KEY --simulation TNG50-1 \
        --snap 99 --subhalo-id 0 -o galaxy.h5
"""
import argparse
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from AMR_grid import AMRGrid

# physical constants (cgs)
MSUN_CGS   = 1.989e33
KPC_CM     = 3.0856775814913673e21
MASS_H_CGS = 1.6726e-24
X_H        = 0.76                       # primordial hydrogen mass fraction
ILLUSTRIS_MASS_CGS = 1.0e10 * MSUN_CGS  # 10^10 Msun / h


def laursen09_ndust(nH, xHI, Z, Z_ref=0.0134, f_ion=0.01):
    nHI, nHII = nH * xHI, nH * (1.0 - xHI)
    return (Z / max(Z_ref, 1e-30)) * (nHI + f_ion * nHII)


# --------------------------------------------------------------------------
def load_illustris_hdf5(filepath, center_kpc=None, boxsize_kpc=None):
    """Return dict of per-cell physical-unit arrays (dust-only fields)."""
    import h5py
    with h5py.File(filepath, "r") as f:
        hdr = dict(f["Header"].attrs)
        h = float(hdr.get("HubbleParam", 0.6774))
        z = float(hdr.get("Redshift", 0.0))
        a = float(hdr.get("Time", 1.0 / (1.0 + z)))
        pt0 = f["PartType0"]
        coords  = pt0["Coordinates"][:]                 # ckpc/h
        density = pt0["Density"][:]                      # (1e10 Msun/h)/(ckpc/h)^3
        Z   = pt0["GFM_Metallicity"][:]          if "GFM_Metallicity" in pt0 else None
        xHI = pt0["NeutralHydrogenAbundance"][:] if "NeutralHydrogenAbundance" in pt0 else None
        sfr = pt0["StarFormationRate"][:]        if "StarFormationRate" in pt0 else None
        vol = pt0["Volume"][:] if "Volume" in pt0 else None
        if vol is None and "Masses" in pt0:
            vol = pt0["Masses"][:] / np.maximum(density, 1e-30)   # (ckpc/h)^3

    # positions ckpc/h -> physical kpc
    xyz = coords * (a / h)
    if center_kpc is not None:
        xyz = xyz - np.asarray(center_kpc)[None, :]
        if boxsize_kpc is not None:
            half = 0.5 * boxsize_kpc
            xyz = np.where(xyz > half, xyz - boxsize_kpc, xyz)
            xyz = np.where(xyz < -half, xyz + boxsize_kpc, xyz)
    # density -> nH [cm^-3]
    rho_cgs = density * (ILLUSTRIS_MASS_CGS * h**2) / (KPC_CM * a)**3
    nH = rho_cgs * X_H / MASS_H_CGS
    # effective Voronoi cell radius [physical kpc]
    r_eff = None
    if vol is not None:
        r_eff = (3.0 * (vol * (a/h)**3) / (4.0*np.pi))**(1.0/3.0)
    return dict(x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2], nH=nH,
                metallicity=Z, xHI=xHI, sfr=sfr, r_eff=r_eff, h=h, a=a, redshift=z)


def download_tng_cutout(simulation, snap_num, subhalo_id, api_key, cache_dir="."):
    import requests
    base = "https://www.tng-project.org/api"
    headers = {"api-key": api_key}
    meta = requests.get(f"{base}/{simulation}/snapshots/{snap_num}/subhalos/{subhalo_id}",
                        headers=headers); meta.raise_for_status()
    cutout_url = meta.json()["cutouts"]["subhalo"]
    params = {"gas": "Coordinates,Density,Masses,GFM_Metallicity,NeutralHydrogenAbundance,StarFormationRate"}
    r = requests.get(cutout_url, headers=headers, params=params); r.raise_for_status()
    out = Path(cache_dir) / f"cutout_{simulation}_snap{snap_num}_sub{subhalo_id}.hdf5"
    out.write_bytes(r.content)
    print(f"downloaded {out} ({out.stat().st_size/1e6:.1f} MB)")
    return str(out)


# --------------------------------------------------------------------------
class VoronoiNN:
    """Nearest-neighbor lookup over the Voronoi cell centers (cKDTree)."""
    def __init__(self, data):
        from scipy.spatial import cKDTree
        self.pts = np.column_stack([data["x"], data["y"], data["z"]])
        self.tree = cKDTree(self.pts)
        self.nH = data["nH"]
        self.Z = data["metallicity"]
        self.xHI = data["xHI"]
        self.r_eff = data["r_eff"]

    def _idx(self, x, y, z):
        q = np.column_stack([np.atleast_1d(x), np.atleast_1d(y), np.atleast_1d(z)])
        return self.tree.query(q)[1]

    def nH_fn(self, x, y, z):   return self.nH[self._idx(x, y, z)]
    def Z_fn(self, x, y, z):    return self.Z[self._idx(x, y, z)]
    def xHI_fn(self, x, y, z):  return self.xHI[self._idx(x, y, z)]
    def size_fn(self, x, y, z): return self.r_eff[self._idx(x, y, z)]


# --------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description="Illustris/TNG -> generic AMR (dust-only).")
    ap.add_argument("input_file", nargs="?", default=None, help="local cutout/snapshot HDF5")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--grid-type", choices=["amr", "uniform"], default="amr")
    ap.add_argument("--level-min", type=int, default=3)
    ap.add_argument("--level-max", type=int, default=7)
    ap.add_argument("--dens-threshold", type=float, default=0.3)
    ap.add_argument("--match-resolution", action="store_true",
                    help="refine until the leaf size matches the local Voronoi cell")
    ap.add_argument("--resolution-factor", type=float, default=1.0)
    ap.add_argument("--center", type=float, nargs=3, default=None,
                    help="recenter on this physical-kpc point (default: density peak)")
    ap.add_argument("--boxsize", type=float, default=None,
                    help="octree box side [kpc] (default: enclose all gas)")
    ap.add_argument("--metallicity", type=float, default=None,
                    help="uniform Z fallback when GFM_Metallicity absent")
    ap.add_argument("--emit-ndust", action="store_true",
                    help="precompute ndust = laursen09(nH,xHI,Z) and write it")
    ap.add_argument("--Z-ref", type=float, default=0.0134)
    ap.add_argument("--f-ion", type=float, default=0.01)
    ap.add_argument("--sfr-mask", choices=["keep", "drop"], default="keep",
                    help="keep or drop star-forming cells (no temperature used)")
    ap.add_argument("--api-key", default=None)
    ap.add_argument("--simulation", default="TNG50-1")
    ap.add_argument("--snap", type=int, default=99)
    ap.add_argument("--subhalo-id", type=int, default=0)
    args = ap.parse_args(argv)

    if args.input_file is None:
        if not args.api_key:
            ap.error("provide a local input_file or --api-key for API download")
        args.input_file = download_tng_cutout(args.simulation, args.snap,
                                               args.subhalo_id, args.api_key)

    data = load_illustris_hdf5(args.input_file)
    n0 = data["x"].size

    # metallicity fallback
    if data["metallicity"] is None:
        if args.metallicity is None:
            print("note: no GFM_Metallicity; metallicity column omitted "
                  "(pass --metallicity Z for a uniform value).")
        else:
            data["metallicity"] = np.full(n0, args.metallicity)
    if data["xHI"] is None:
        print("note: no NeutralHydrogenAbundance; assuming xHI=1 (fully neutral).")
        data["xHI"] = np.ones(n0)

    # optional SFR mask (density-based, no temperature)
    if args.sfr_mask == "drop" and data["sfr"] is not None:
        keep = data["sfr"] <= 0.0
        for k in ("x", "y", "z", "nH", "metallicity", "xHI", "r_eff", "sfr"):
            if data[k] is not None:
                data[k] = data[k][keep]
        print(f"sfr-mask=drop: kept {data['x'].size}/{n0} non-star-forming cells")

    # center + box
    if args.center is not None:
        ctr = np.array(args.center, float)
    else:
        ctr = np.array([data["x"][np.argmax(data["nH"])],
                        data["y"][np.argmax(data["nH"])],
                        data["z"][np.argmax(data["nH"])]])
    for i, k in enumerate(("x", "y", "z")):
        data[k] = data[k] - ctr[i]
    if args.boxsize is not None:
        box = args.boxsize
        sel = ((np.abs(data["x"]) <= 0.5*box) & (np.abs(data["y"]) <= 0.5*box) &
               (np.abs(data["z"]) <= 0.5*box))
        for k in ("x", "y", "z", "nH", "metallicity", "xHI", "r_eff", "sfr"):
            if data[k] is not None:
                data[k] = data[k][sel]
        print(f"boxsize={box} kpc: {data['x'].size} cells inside the box")
    else:
        ext = max(np.abs(data["x"]).max(), np.abs(data["y"]).max(),
                  np.abs(data["z"]).max())
        box = 2.05 * ext
    print(f"center={ctr} kpc, boxlen={box:.3f} kpc, {data['x'].size} gas cells")

    interp = VoronoiNN(data)
    grid = AMRGrid(boxlen=box, origin=(-0.5*box, -0.5*box, -0.5*box))

    if args.grid_type == "uniform":
        grid.refine_uniform(args.level_max)
    else:
        grid.refine_by_density(interp.nH_fn, threshold=args.dens_threshold,
                               level_max=args.level_max, level_min=args.level_min)
        if args.match_resolution and data["r_eff"] is not None:
            grid.refine_to_resolution(interp.size_fn, level_max=args.level_max,
                                      factor=args.resolution_factor,
                                      level_min=args.level_min)

    grid.set_density(interp.nH_fn)
    if data["metallicity"] is not None:
        grid.set_metallicity(interp.Z_fn)
    grid.set_neutral_fraction(interp.xHI_fn)
    if args.emit_ndust and data["metallicity"] is not None:
        def ndust_fn(x, y, z):
            return laursen09_ndust(interp.nH_fn(x, y, z), interp.xHI_fn(x, y, z),
                                   interp.Z_fn(x, y, z), args.Z_ref, args.f_ion)
        grid.set_dust_density(ndust_fn)

    print(grid.info())
    grid.write(args.output)
    print("Hint: MoCafe namelist -> par%grid_type='amr', par%amr_file='"
          + args.output + "',")
    print("      par%dust_model='laursen09' (metallicity+xHI) | 'from_file' "
          "(--emit-ndust) | 'global_dgr' (nH only); par%distance_unit='kpc'.")


if __name__ == "__main__":
    main()
