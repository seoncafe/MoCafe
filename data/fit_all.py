#!/usr/bin/env python3
"""Combined Mueller-matrix fit for all four components, ported from fit_all.pro.

Fits S11 with a 3-parameter double-HG and then S12/S11, S33/S11, S34/S11
with their model functions, plotting six panels (normalized fit and
unnormalized fit*S11 for each of S12, S33, S34).

Usage:
    python fit_all.py [filename]      # default: mueller_matrix.dat
"""
from __future__ import annotations
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from mueller_funcs import (read_mueller_file,
                           func_s11, func_s12, func_s33, func_s34)


def main(fname='mueller_matrix.dat'):
    _meta, cost, s11, s12, s33, s34 = read_mueller_file(fname)

    fig, axes = plt.subplots(3, 2, figsize=(10, 11))
    xr_full = (-1.0, 1.0)

    # --- S11 (used to multiply the normalized fits below) ---
    p11, _ = curve_fit(func_s11, cost, s11, p0=[0.7, 0.8, 0.4])
    yfit_s11 = func_s11(cost, *p11)
    print('S11   parameter', p11)

    # --- S12/S11 ---
    p12, _ = curve_fit(func_s12, cost, s12 / s11, p0=[-0.4])
    yfit12 = func_s12(cost, *p12)
    print('S12   parameter', p12)
    axes[0, 0].plot(cost, s12 / s11, 'k.', ms=2)
    axes[0, 0].plot(cost, yfit12, 'r-')
    axes[0, 0].set(xlim=xr_full, xlabel=r'$\cos\theta$', ylabel='S12 / S11')
    axes[0, 1].plot(cost, s12, 'k.', ms=2)
    axes[0, 1].plot(cost, yfit12 * yfit_s11, 'r-')
    axes[0, 1].set(xlim=xr_full, ylim=(-0.5, 0.0),
                   xlabel=r'$\cos\theta$', ylabel='S12')

    # --- S33/S11 ---
    p33, _ = curve_fit(func_s33, cost, s33 / s11, p0=[1.0])
    yfit33 = func_s33(cost, *p33)
    print('S33   parameter', p33)
    s33_ref = 2.0 * cost / (1.0 + np.sqrt(np.abs(cost)))
    axes[1, 0].plot(cost, s33 / s11, 'k.', ms=2)
    axes[1, 0].plot(cost, yfit33, 'r-')
    axes[1, 0].plot(cost, s33_ref, 'g--')
    axes[1, 0].set(xlim=xr_full, xlabel=r'$\cos\theta$', ylabel='S33 / S11')
    axes[1, 1].plot(cost, s33, 'k.', ms=2)
    axes[1, 1].plot(cost, yfit33 * yfit_s11, 'r-')
    axes[1, 1].set(xlim=xr_full, xlabel=r'$\cos\theta$', ylabel='S33')

    # --- S34/S11 (full 5-parameter form).  The shape is heavily
    # Lyman-alpha-tuned; on continuum-dust tables (mueller_V.dat etc.)
    # it may not converge.  Plot data even if the fit fails.
    axes[2, 0].plot(cost, s34 / s11, 'k.', ms=2)
    axes[2, 0].set(xlim=xr_full, xlabel=r'$\cos\theta$', ylabel='S34 / S11')
    axes[2, 1].plot(cost, s34, 'k.', ms=2)
    axes[2, 1].set(xlim=xr_full, xlabel=r'$\cos\theta$', ylabel='S34')
    try:
        p34, _ = curve_fit(func_s34, cost, s34 / s11,
                           p0=[0.39, 1.0, 3.0, 0.0, 1.0], maxfev=20000)
        yfit34 = func_s34(cost, *p34)
        print('S34   parameter', p34)
        axes[2, 0].plot(cost, yfit34, 'r-')
        axes[2, 1].plot(cost, yfit34 * yfit_s11, 'r-')
    except RuntimeError as e:
        print('S34   parameter  FIT FAILED:', e)

    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'mueller_matrix.dat')
