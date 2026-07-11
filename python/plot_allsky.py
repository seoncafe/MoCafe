#!/usr/bin/env python3
"""Plot a MoCafe ``_allsky`` HEALPix map (the interior-observer all-sky output)
as a Mollweide projection with healpy.

Shows the surface brightness at a chosen wavelength, or the wavelength-integrated
(bolometric) map with ``--bol``.

Usage:
    python plot_allsky.py <run>               # stem or an _allsky file
    python plot_allsky.py <run> --lam 100     # nearest bin to 100 um
    python plot_allsky.py <run> --bol -o sky.pdf
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import healpy as hp

import mocafe_io as mio


def resolve(run):
    if os.path.isfile(run):
        return run
    f = mio.find_mocafe_file(run, "_allsky")
    if f is None:
        raise SystemExit(f"no _allsky output found for {run!r}")
    return f


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run", help="run stem or an _allsky file")
    ap.add_argument("-o", "--out", default=None)
    ap.add_argument("--lam", type=float, default=None, help="wavelength [um] (nearest bin)")
    ap.add_argument("--bol", action="store_true", help="wavelength-integrated map")
    args = ap.parse_args()

    path = resolve(args.run)
    a = mio.load_mocafe(path).allsky()
    if a is None:
        raise SystemExit(f"{path} has no AllSky map")

    sky = a["sky"]                       # (nlam, npix), RING
    lam = a["wavelength"]
    if sky.ndim == 1:                    # single-wavelength file
        m, title = sky, "all-sky"
    elif args.bol or lam is None:
        m = sky.sum(axis=0)
        title = r"$\int I_\lambda\,d\lambda$"
    else:
        target = args.lam if args.lam is not None else lam[len(lam) // 2]
        il = int(np.argmin(np.abs(lam - target)))
        m = sky[il]
        title = rf"$\lambda = {lam[il]:.3g}\ \mu$m"

    m = np.asarray(m, dtype=float)
    pos = m[m > 0]
    vmin = pos.min() if pos.size else None

    fig = plt.figure(figsize=(8.0, 5.0))
    hp.mollview(m, fig=fig.number, title=f"{os.path.basename(path)}   {title}",
                unit="surface brightness", norm="log", min=vmin, cmap="inferno")
    hp.graticule()

    out = args.out or (mio._strip_mocafe_ext(os.path.basename(path)) + "_allsky.pdf")
    fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
