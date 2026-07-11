"""Plot the radiation field in each cell from a MoCafe ``_jlam`` output.

Left: a central slice of the wavelength-integrated mean intensity J_bol.
Right: the J_lambda spectrum in the brightest cell.

Usage:
    python plot_jlam.py <run>                  # stem or a _jlam file
    python plot_jlam.py <run> --axis z -o jlam.pdf
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
    f = mio.find_mocafe_file(run, "_jlam")
    if f is None:
        raise SystemExit(f"no _jlam output found for {run!r}")
    return f


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run", help="run stem or a _jlam file")
    ap.add_argument("-o", "--out", default=None)
    ap.add_argument("--axis", default="z", choices=list(AX))
    args = ap.parse_args()

    path = resolve(args.run)
    j = mio.load_mocafe(path).jlambda()
    if j is None:
        raise SystemExit(f"{path} has no J_lambda")

    J = j["J"]              # Cartesian: (nx,ny,nz,nlam); AMR: (nleaf,nlam)
    jbol = j["jbol"]
    lam = j["wavelength"]

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4))

    if jbol.ndim == 3:      # Cartesian
        n = jbol.shape[AX[args.axis]]
        sl = np.take(jbol, n // 2, axis=AX[args.axis]).T
        pos = sl[sl > 0]
        norm = matplotlib.colors.LogNorm(vmin=pos.min(), vmax=sl.max()) if pos.size else None
        im = axes[0].imshow(np.where(sl > 0, sl, np.nan), origin="lower", cmap="magma", norm=norm)
        axes[0].set_title(r"$J_{\rm bol}$" + f" ({args.axis}-slice)")
        axes[0].set_xlabel("cell"); axes[0].set_ylabel("cell")
        fig.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04)
        pk = np.unravel_index(np.argmax(jbol), jbol.shape)
        spec = J[pk[0], pk[1], pk[2], :]
        cell_lbl = f"peak cell {pk}"
    else:                   # AMR
        ipk = int(np.argmax(jbol))
        spec = J[ipk, :]
        axes[0].hist(np.log10(jbol[jbol > 0]), bins=40, color="#555555")
        axes[0].set_xlabel(r"$\log_{10} J_{\rm bol}$"); axes[0].set_ylabel("leaves")
        axes[0].set_title(r"$J_{\rm bol}$ distribution")
        cell_lbl = f"peak leaf {ipk}"

    axes[1].plot(lam, lam * spec, color="#1f77b4", lw=1.6)
    axes[1].set_xscale("log"); axes[1].set_yscale("log")
    axes[1].set_xlabel(r"$\lambda\ [\mu\mathrm{m}]$")
    axes[1].set_ylabel(r"$\lambda J_\lambda$")
    axes[1].set_title(cell_lbl)
    axes[1].grid(True, which="both", ls=":", lw=0.4, alpha=0.5)

    fig.suptitle(os.path.basename(path))
    out = args.out or (mio._strip_mocafe_ext(os.path.basename(path)) + "_jlam.pdf")
    fig.tight_layout()
    fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
