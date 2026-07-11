#!/usr/bin/env python3
"""Generate stellar / interstellar source spectra for MoCafe SED runs.

MoCafe reads a source spectrum as a 2-column table (lambda[um], L_lambda) via
par%source_spectrum or par%src_spectrum(i); only the *shape* is used (it is
normalized to a per-bin sampling PDF, and the absolute scale is set by
par%luminosity / par%src_lum).  This script writes such tables for:

  spec_young_fsps.dat / spec_old_fsps.dat   -- FSPS SSP, 10 Myr / 10 Gyr
  spec_young_sb99.dat / spec_old_sb99.dat   -- Starburst99 SSP, 10 Myr / 10 Gyr
  isrf_mathis.dat                           -- Mathis (1983) MW ISRF (starlight)

Sources:
  FSPS         : import fsps (python-fsps).
  Starburst99  : ~/RT_Codes/starburst99 SSP spectrum output.
  Mathis ISRF  : the piecewise field of Mathis, Mezger & Panagia (1983),
                 mirroring ~/MoCafe/Grain/SEDust (radfield.f90 :: J_Mathis).

Run from anywhere; writes into ../data relative to this file unless --outdir.
"""
from __future__ import annotations
import argparse
import os
import numpy as np

CLIGHT_UM_HZ = 2.99792458e14   # c in um/s (for L_nu -> L_lambda)
HC_OVER_K_UM = 1.43877687750393e4   # h c / k  [um K]


def planck_shape(lam_um, T):
    """B_lambda shape (arbitrary normalization), lambda in um."""
    lam = np.asarray(lam_um, dtype=np.float64)
    x = HC_OVER_K_UM / (lam * T)
    b = np.zeros_like(lam)
    ok = x < 700.0
    b[ok] = 1.0 / (lam[ok] ** 5 * (np.exp(x[ok]) - 1.0))
    return b


def j_mathis_starlight(lam_um):
    """Mathis (1983) ISRF starlight component (no CMB), U = 1.

    Piecewise in lambda[um]; returns J_lambda in the paper's units.  The CMB
    term is intentionally omitted -- this is a stellar *source* spectrum, and
    the dust far-infrared emission is computed by MoCafe itself.
    """
    lam = np.asarray(lam_um, dtype=np.float64)
    J = np.zeros_like(lam)
    m2 = (lam >= 0.0912) & (lam < 0.110)
    m3 = (lam >= 0.110) & (lam < 0.134)
    m4 = (lam >= 0.134) & (lam < 0.250)
    m5 = lam >= 0.250
    J[m2] = 3069.0 * lam[m2] ** 3.4172
    J[m3] = 1.627
    J[m4] = 0.0566 * lam[m4] ** (-1.6678)
    J[m5] = (1.0e-14 * planck_shape(lam[m5], 7500.0)
             + 1.0e-13 * planck_shape(lam[m5], 4000.0)
             + 4.0e-13 * planck_shape(lam[m5], 3000.0))
    return J


def write_table(path, lam_um, Llam, header):
    """Write a 2-column (lambda[um], L_lambda) table, ascending in lambda,
    normalized to unit peak (shape only)."""
    lam_um = np.asarray(lam_um, dtype=np.float64)
    Llam = np.asarray(Llam, dtype=np.float64)
    order = np.argsort(lam_um)
    lam_um, Llam = lam_um[order], Llam[order]
    good = np.isfinite(Llam) & (Llam >= 0.0)
    lam_um, Llam = lam_um[good], Llam[good]
    peak = Llam.max()
    if peak > 0:
        Llam = Llam / peak
    with open(path, "w") as f:
        f.write("# " + header + "\n")
        f.write("# lambda[um]   L_lambda[arb, shape only; normalized to unit peak]\n")
        for x, y in zip(lam_um, Llam):
            f.write(f"{x:14.6e} {y:14.6e}\n")
    print(f"  wrote {path}  ({len(lam_um)} rows, {lam_um.min():.4g}-{lam_um.max():.4g} um)")


# ---------------------------------------------------------------- FSPS
def make_fsps(outdir):
    import fsps
    print("FSPS: building SSP (this loads isochrones, ~1 min the first time)...")
    sp = fsps.StellarPopulation(zcontinuous=1, sfh=0, logzsol=0.0, imf_type=2)
    for age_gyr, tag, label in [(0.01, "young", "10 Myr"), (10.0, "old", "10 Gyr")]:
        wave_aa, spec_nu = sp.get_spectrum(tage=age_gyr)   # AA, Lsun/Hz per Msun
        lam_um = wave_aa / 1.0e4
        # L_nu -> L_lambda: L_lambda = L_nu * c / lambda^2
        Llam = spec_nu * CLIGHT_UM_HZ / lam_um ** 2
        write_table(os.path.join(outdir, f"spec_{tag}_fsps.dat"), lam_um, Llam,
                    f"FSPS SSP (MIST/MILES, solar Z, Kroupa IMF), {label} -- {tag} population")


# ---------------------------------------------------------------- Starburst99
def load_sb99(path):
    """Return (times[yr], wave[AA], logtot) arrays from an SB99 .spectrum1."""
    t, w, lt = [], [], []
    with open(path) as fh:
        for ln in fh:
            p = ln.split()
            if len(p) >= 4:
                try:
                    t.append(float(p[0])); w.append(float(p[1])); lt.append(float(p[2]))
                except ValueError:
                    continue
    return np.array(t), np.array(w), np.array(lt)


def make_sb99(outdir, sb99_file):
    if not os.path.isfile(sb99_file):
        print(f"SB99: file not found ({sb99_file}); skipping.")
        return
    print(f"SB99: reading {sb99_file}")
    t, w, logtot = load_sb99(sb99_file)
    times = np.unique(t)
    for age_yr, tag, label in [(1.0e7, "young", "10 Myr"), (1.0e10, "old", "10 Gyr")]:
        it = times[np.argmin(np.abs(times - age_yr))]
        sel = t == it
        lam_um = w[sel] / 1.0e4
        Llam = 10.0 ** logtot[sel]              # erg/s/AA (shape)
        write_table(os.path.join(outdir, f"spec_{tag}_sb99.dat"), lam_um, Llam,
                    f"Starburst99 instantaneous SSP (Padova, Z=0.02), t={it:.3e} yr (~{label}) -- {tag} population")


# ---------------------------------------------------------------- Mathis ISRF
def make_mathis(outdir):
    print("Mathis (1983) MW ISRF starlight")
    lam_um = np.logspace(np.log10(0.0912), np.log10(2000.0), 400)
    Llam = j_mathis_starlight(lam_um)
    write_table(os.path.join(outdir, "isrf_mathis.dat"), lam_um, Llam,
                "Mathis, Mezger & Panagia (1983) MW ISRF starlight (no CMB), U=1")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=os.path.join(here, "..", "data"))
    ap.add_argument("--sb99-file",
                    default=os.path.expanduser(
                        "~/RT_Codes/starburst99/SSP_instant_Padova_STD/Z0.02.spectrum1"))
    ap.add_argument("--skip-fsps", action="store_true")
    ap.add_argument("--skip-sb99", action="store_true")
    args = ap.parse_args()
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    print(f"output directory: {outdir}\n")
    make_mathis(outdir)
    if not args.skip_sb99:
        make_sb99(outdir, args.sb99_file)
    if not args.skip_fsps:
        make_fsps(outdir)
    print("\ndone.")


if __name__ == "__main__":
    main()
