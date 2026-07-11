#!/usr/bin/env python3
"""Plot the panchromatic SED of a MoCafe external-galaxy run.

Overlays the intrinsic stellar input spectrum and the dust thermal emission,
both as lambda*L_lambda [erg/s].  The stellar input is reconstructed from the
source parameters in the run's namelist (a sum of blackbodies, or a tabulated
source_spectrum); the dust emission is ``SED_intrinsic`` from the ``_dustsed``
file.  By energy conservation the dust integral equals the absorbed starlight,
so the plot shows the UV/optical energy reprocessed into the infrared (PAH
bands + FIR peak).

Usage:
    python plot_galaxy_sed.py <run>            # run stem (e.g. galaxy_sed)
    python plot_galaxy_sed.py <run> -o galaxy_sed.pdf
"""
import argparse
import glob
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import mocafe_io as mio

H, C, K = 6.62607015e-27, 2.99792458e10, 1.380649e-16   # cgs


def _nml_values(text, key):
    """Return the comma-separated values of `par%<key>` from namelist text."""
    m = re.search(rf"par%{key}\s*=\s*([^\n!]+)", text, re.IGNORECASE)
    if not m:
        return []
    out = []
    for tok in m.group(1).split(","):
        tok = tok.strip().strip("'\"")
        if tok:
            out.append(tok)
    return out


def parse_sources(in_path):
    """Return a list of stellar components from a .in file.

    Each component is ``('file', path, L)`` for a tabulated source spectrum or
    ``('bb', T, L)`` for a blackbody.  Spectrum-file paths are resolved relative
    to the namelist location.  L is the component luminosity (``src_lum``, or the
    single ``luminosity``).
    """
    if not in_path or not os.path.isfile(in_path):
        return []
    txt = open(in_path).read()
    base = os.path.dirname(os.path.abspath(in_path))
    ls = [float(x) for x in _nml_values(txt, "src_lum")] \
        or [float(x) for x in _nml_values(txt, "luminosity")]

    def _lum(i):
        return ls[i] if i < len(ls) else (ls[-1] if ls else 1.0)

    specs = _nml_values(txt, "src_spectrum") or _nml_values(txt, "source_spectrum")
    if specs:
        return [("file", s if os.path.isabs(s) else os.path.join(base, s), _lum(i))
                for i, s in enumerate(specs)]
    ts = _nml_values(txt, "src_tstar") or _nml_values(txt, "tstar")
    return [("bb", float(t), _lum(i)) for i, t in enumerate(ts) if float(t) > 0]


def blackbody_llam(lam_um, T):
    lam_cm = lam_um * 1e-4
    x = H * C / (lam_cm * K * T)
    with np.errstate(over="ignore", divide="ignore"):
        B = 1.0 / (lam_cm ** 5) / np.expm1(x)
    return np.where(np.isfinite(B), B, 0.0)


def _component_shape(lam_um, kind, data):
    """L_lambda shape of one component (blackbody or tabulated file)."""
    if kind == "file":
        tbl = np.loadtxt(data, comments="#")
        return np.interp(lam_um, tbl[:, 0], tbl[:, 1], left=0.0, right=0.0)
    return blackbody_llam(lam_um, data)


def stellar_input(lam_um, comps):
    """Sum of components (blackbody or file), each normalized to its L."""
    Llam = np.zeros_like(lam_um)
    for kind, data, L in comps:
        shape = _component_shape(lam_um, kind, data)
        norm = np.trapz(shape, lam_um)
        if norm > 0:
            Llam += L * shape / norm
    return Llam                              # erg/s per um


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run", help="run stem (e.g. galaxy_sed)")
    ap.add_argument("-o", "--out", default=None)
    ap.add_argument("--in", dest="infile", default=None,
                    help="namelist for the source spectrum (default: <stem>.in)")
    args = ap.parse_args()

    stem = mio._strip_mocafe_ext(args.run)
    dpath = mio.find_mocafe_file(stem, "_dustsed") or (stem + "_dustsed.h5")
    if not os.path.isfile(dpath):
        raise SystemExit(f"no _dustsed file for {stem!r}")
    md = mio.load_mocafe(dpath)
    d = md.dust_sed()
    lam, dlam, intr = d["wavelength"], d["dwavelength"], d["intrinsic"]
    sec = md.section("SED_emergent")
    L = float(sec.attr("TOT_LUM")) if sec and sec.attr("TOT_LUM") else None
    Labs = float(sec.attr("L_ABS")) if sec and sec.attr("L_ABS") else float(intr.sum())

    lLl_dust = lam * np.where(dlam > 0, intr / dlam, 0.0)

    fig, ax = plt.subplots(figsize=(7.4, 4.8))

    in_path = args.infile or (stem + ".in")
    comps = parse_sources(in_path)
    if comps:
        lLl_star = lam * stellar_input(lam, comps)
        ax.plot(lam, lLl_star, color="#1f77b4", lw=2.0, label="stellar input")

    ax.plot(lam, lLl_dust, color="#d62728", lw=2.0, label="dust emission")
    ax.fill_between(lam, lLl_dust, lLl_dust.max() * 1e-6, color="#d62728", alpha=0.08)

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.09, 2.0e3)
    ymax = ax.get_ylim()[1]
    ax.set_ylim(ymax * 1e-4, ymax)
    ax.set_xlabel(r"$\lambda\ [\mu\mathrm{m}]$")
    ax.set_ylabel(r"$\lambda L_\lambda\ [\mathrm{erg\,s^{-1}}]$")
    if L:
        ax.text(0.03, 0.06, rf"$L_{{\rm dust}}/L_\star = {Labs/L:.2f}$",
                transform=ax.transAxes, fontsize=10)
    ax.legend(frameon=False, loc="upper right")
    ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)

    out = args.out or (os.path.basename(stem) + "_sed.pdf")
    fig.tight_layout()
    fig.savefig(out)
    print("wrote", out, f"(L_dust/L_star = {Labs/L:.3f})" if L else "")


if __name__ == "__main__":
    main()
