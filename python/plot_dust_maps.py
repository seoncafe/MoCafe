#!/usr/bin/env python3
"""Plot the dust temperature and power maps from a MoCafe ``_dustsed`` (Lucy)
or ``_bwdust`` (Bjorkman & Wood) output.

On a Cartesian grid it shows a central slice; on an AMR octree it shows a
scatter of the leaf centers colored by value.

Usage:
    python plot_dust_maps.py <run>            # stem, or a _dustsed / _bwdust file
    python plot_dust_maps.py <run> --axis z --slice mid -o maps.pdf
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import mocafe_io as mio

AX = {"x": 0, "y": 1, "z": 2}


def resolve(run):
    if os.path.isfile(run):
        return run
    for suf in ("_dustsed", "_bwdust"):
        f = mio.find_mocafe_file(run, suf)
        if f is not None:
            return f
    raise SystemExit(f"no _dustsed / _bwdust output found for {run!r}")


def slice2d(cube, axis, idx):
    return np.take(cube, idx, axis=AX[axis])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run", help="run stem or a _dustsed / _bwdust file")
    ap.add_argument("-o", "--out", default=None)
    ap.add_argument("--axis", default="z", choices=list(AX), help="slice normal (Cartesian)")
    ap.add_argument("--slice", default="mid", help="'mid' or an integer index (Cartesian)")
    args = ap.parse_args()

    path = resolve(args.run)
    mf = mio.load_mocafe(path)
    m = mf.dust_maps()
    if m is None:
        raise SystemExit(f"{path} has no Tdust map")

    T = m["tdust"]
    P = m["ldust"] if m["ldust"] is not None else m["labs"]
    plabel = "Ldust" if m["ldust"] is not None else "Labs"

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4))

    if T.ndim == 3:  # Cartesian
        n = T.shape[AX[args.axis]]
        idx = n // 2 if args.slice == "mid" else int(args.slice)
        for ax, dat, name, cmap, logscale in (
            (axes[0], slice2d(T, args.axis, idx), r"$T_{\rm dust}\ [{\rm K}]$", "inferno", False),
            (axes[1], slice2d(P, args.axis, idx), plabel + r"$\ [{\rm erg\,s^{-1}}]$", "viridis", True),
        ):
            d = dat.T
            if logscale:
                pos = d[d > 0]
                vmin = pos.min() if pos.size else 1e-30
                im = ax.imshow(np.where(d > 0, d, np.nan), origin="lower",
                               cmap=cmap, norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=d.max()))
            else:
                im = ax.imshow(d, origin="lower", cmap=cmap)
            ax.set_title(name)
            ax.set_xlabel("cell"); ax.set_ylabel("cell")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        fig.suptitle(f"{os.path.basename(path)}  ({args.axis}-slice {idx})")
    else:  # AMR: leaf scatter
        xyz = m["leafxyz"]
        if xyz is None:
            raise SystemExit("AMR maps need a LeafXYZ table")
        xyz = np.asarray(xyz)
        for ax, val, name, cmap, logscale in (
            (axes[0], T, r"$T_{\rm dust}\ [{\rm K}]$", "inferno", False),
            (axes[1], P, plabel + r"$\ [{\rm erg\,s^{-1}}]$", "viridis", True),
        ):
            c = val
            norm = matplotlib.colors.LogNorm(vmin=c[c > 0].min(), vmax=c.max()) if logscale else None
            sc = ax.scatter(xyz[:, 0], xyz[:, 1], c=c, s=6, cmap=cmap, norm=norm)
            ax.set_title(name); ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_aspect("equal")
            fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        fig.suptitle(f"{os.path.basename(path)}  (AMR leaves)")

    out = args.out or (mio._strip_mocafe_ext(os.path.basename(path)) + "_maps.pdf")
    fig.tight_layout()
    fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
