"""Mueller-matrix model functions for fitting dust scattering tables.

Direct Python port of the IDL routines `func_s11.pro`, `func_s12.pro`,
`func_s33.pro`, `func_s34.pro`.  Each model takes
    cos(theta)  : independent variable (1-D array)
    *p          : variable number of scalar parameters
and returns the model value at each cos(theta).  The signature matches
scipy.optimize.curve_fit (positional scalar parameters after `x`).

Also provides `read_mueller_file(fname)`, a Python replacement for the
IDL `readcol` calls that pull the header metadata and the
cos(theta)/S11/S12/S33/S34 columns out of a `mueller_*.dat` file.
"""
from __future__ import annotations
import numpy as np


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def read_mueller_file(fname):
    """Read a tabulated Mueller-matrix file.

    File layout (matches the data/mueller_*.dat files):
        line 1  : column labels for metadata (text, ignored)
        line 2  : wavelength[um], Cext[cm^2/H], albedo, <cos>, nang
        line 3  : column labels for the table (text, ignored)
        line 4+ : cos(theta), S11, S12, S33, S34

    Returns
    -------
    meta : dict with keys 'wavl', 'cext', 'albedo', 'hgg', 'nang'
    cost, s11, s12, s33, s34 : 1-D ndarrays
    """
    with open(fname) as f:
        lines = f.readlines()
    meta_vals = lines[1].split()
    meta = dict(
        wavl=float(meta_vals[0]),
        cext=float(meta_vals[1]),
        albedo=float(meta_vals[2]),
        hgg=float(meta_vals[3]),
        nang=int(meta_vals[4]),
    )
    data = np.loadtxt(fname, skiprows=3)
    cost = data[:, 0]
    s11 = data[:, 1]
    s12 = data[:, 2]
    s33 = data[:, 3]
    s34 = data[:, 4]
    return meta, cost, s11, s12, s33, s34


# ---------------------------------------------------------------------------
# Model functions
# ---------------------------------------------------------------------------
def func_s11(x, *p):
    """Henyey-Greenstein S11 model with 1, 2, 3, 4, or 5 parameters.

    np = 1 : single-HG, fixed unit scale         p = (g,)
    np = 2 : single-HG with free amplitude       p = (scale, g)
    np = 3 : double-HG, second scale = 1-first   p = (scale, g, g2)
    np = 4 : double-HG, both scales free         p = (scale, g, scale2, g2)
    np = 5 : double-HG + isotropic offset        p = (scale, g, scale2, g2, c)
    """
    np_ = len(p)
    if np_ == 1:
        scale, g = 1.0, p[0]
    else:
        scale, g = abs(p[0]), p[1]
    y = scale * 0.5 * (1.0 - g**2) / (1.0 + g**2 - 2.0 * g * x)**1.5

    if np_ == 3:
        scale2, g2 = 1.0 - scale, p[2]
    elif np_ >= 4:
        scale2, g2 = abs(p[2]), p[3]

    if np_ > 2:
        y = y + scale2 * 0.5 * (1.0 - g2**2) / (1.0 + g2**2 - 2.0 * g2 * x)**1.5

    if np_ == 5:
        y = y + p[4] / 0.5

    return y


def func_s12(x, *p):
    """S12/S11 model: a simple (1 - x^2) shape with one amplitude."""
    scale = p[0]
    return scale * (1.0 - x**2)


def func_s33(x, *p):
    """S33/S11 model: 2x / (1 + p0 * sqrt(|x|))."""
    return 2.0 * x / (1.0 + p[0] * np.sqrt(np.abs(x)))


def func_s34(x, *p):
    """S34/S11 model used for Lyman-alpha-like polarization curves.

    Parameters: pcir, s, se, th0, ex.
    The (theta-th0)**ex term is evaluated as |theta-th0|**ex so that
    fractional exponents explored by the optimizer don't produce NaN
    from negative bases (matches the way IDL's mpfit happened to walk
    this corner of parameter space).
    """
    pcir, s, se, th0, ex = p
    theta = np.degrees(np.arccos(np.clip(x, -1.0, 1.0)))
    arg = theta + s * 3.13 * theta * np.exp(
        -se * np.abs(theta - th0)**ex / 180.0)
    cc = np.cos(np.radians(arg))
    return pcir * (1.0 - cc**2) / (1.0 + cc**2)
