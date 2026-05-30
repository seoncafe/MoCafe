#!/usr/bin/env python3
"""S12/S33/S34 Mueller-matrix fits, ported from fit_s12.pro.

Each component is fit against its S11-normalized form, then the unnormalized
fit (yfit*S11) is overplotted on the raw component.  Six panels in a 2x3
grid: (S12/S11, S12), (S33/S11, S33), (S34/S11, S34).

Usage:
    python fit_s12.py [filename]      # default: mueller_matrix.dat
"""
from __future__ import annotations
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from mueller_funcs import read_mueller_file, func_s12, func_s33, func_s34


def main(fname='mueller_matrix.dat'):
    _meta, cost, s11, s12, s33, s34 = read_mueller_file(fname)

    fig, axes = plt.subplots(3, 2, figsize=(10, 11))
    xr_full = (-1.0, 1.0)
    xr_zoom = (0.5, 1.0)

    # --- S12 fit ---
    p12, _ = curve_fit(func_s12, cost, s12 / s11, p0=[-0.4])
    yfit12 = func_s12(cost, *p12)
    print('S12   parameter', p12)
    axes[0, 0].plot(cost, s12 / s11, 'k.', ms=2)
    axes[0, 0].plot(cost, yfit12, 'r-')
    axes[0, 0].set(xlim=xr_full, xlabel=r'$\cos\theta$', ylabel='S12 / S11')
    axes[0, 1].plot(cost, s12, 'k.', ms=2)
    axes[0, 1].plot(cost, yfit12 * s11, 'r-')
    axes[0, 1].set(xlim=xr_zoom, xlabel=r'$\cos\theta$', ylabel='S12')

    # --- S33 fit ---
    p33, _ = curve_fit(func_s33, cost, s33 / s11, p0=[1.0])
    yfit33 = func_s33(cost, *p33)
    print('S33   parameter', p33)
    s33_ref = 2.0 * cost / (1.0 + np.sqrt(np.abs(cost)))
    axes[1, 0].plot(cost, s33 / s11, 'k.', ms=2)
    axes[1, 0].plot(cost, yfit33, 'r-')
    axes[1, 0].plot(cost, s33_ref, 'g--')
    axes[1, 0].set(xlim=xr_full, xlabel=r'$\cos\theta$', ylabel='S33 / S11')
    axes[1, 1].plot(cost, s33, 'k.', ms=2)
    axes[1, 1].plot(cost, yfit33 * s11, 'r-')
    axes[1, 1].set(xlim=xr_zoom, xlabel=r'$\cos\theta$', ylabel='S33')

    # --- S34 fit (2-parameter abbreviated form from fit_s12.pro) ---
    # func_s34 takes 5 params; pad with reasonable defaults that don't
    # affect the fit much (se=1, th0=0, ex=1 mimic the simpler exponential).
    def s34_2par(x, p0, p1):
        return func_s34(x, p0, p1, 1.0, 0.0, 1.0)

    p34, _ = curve_fit(s34_2par, cost, s34 / s11, p0=[0.39, 1.0])
    yfit34 = s34_2par(cost, *p34)
    print('S34   parameter', p34)
    axes[2, 0].plot(cost, s34 / s11, 'k.', ms=2)
    axes[2, 0].plot(cost, yfit34, 'r-')
    axes[2, 0].set(xlim=xr_full, xlabel=r'$\cos\theta$', ylabel='S34 / S11')
    axes[2, 1].plot(cost, s34, 'k.', ms=2)
    axes[2, 1].plot(cost, yfit34 * s11, 'r-')
    axes[2, 1].set(xlim=xr_zoom, xlabel=r'$\cos\theta$', ylabel='S34')

    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'mueller_matrix.dat')
