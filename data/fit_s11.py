#!/usr/bin/env python3
"""S11-focused Mueller-matrix fit, ported from fit_s11.pro.

Fits the S11 phase function with a 1-parameter Henyey-Greenstein and then
a 3-parameter double-HG, plotting the comparison against the canonical
single-HG curve computed from the file's <cos> metadata.

Usage:
    python fit_s11.py [filename]      # default: mueller_matrix.dat
"""
from __future__ import annotations
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from mueller_funcs import read_mueller_file, func_s11


def main(fname='mueller_matrix.dat'):
    meta, cost, s11, _s12, _s33, _s34 = read_mueller_file(fname)
    hgg = meta['hgg']
    hgfun = 0.5 * (1.0 - hgg**2) / (1.0 + hgg**2 - 2.0 * hgg * cost)**1.5

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    # one-parameter HG fit (g only).
    p1, _ = curve_fit(func_s11, cost, s11, p0=[hgg])
    yfit1 = func_s11(cost, *p1)
    print('one   parameter', p1)

    # three-parameter double-HG fit (scale, g, g2).
    p3, _ = curve_fit(func_s11, cost, s11, p0=[0.7, 0.8, 0.4])
    yfit3 = func_s11(cost, *p3)
    print('three parameters', p3)

    for ax, ylog, xr in [(axes[0, 0], False, (0.5, 1.0)),
                         (axes[0, 1], True,  (-1.0, 1.0)),
                         (axes[1, 0], False, (0.5, 1.0)),
                         (axes[1, 1], True,  (-1.0, 1.0))]:
        ax.plot(cost, s11, 'k.', ms=2, label='data')
        if ax in (axes[0, 0], axes[0, 1]):
            ax.plot(cost, yfit1, 'r-', label='1-par HG fit')
        else:
            ax.plot(cost, yfit3, 'r-', label='3-par double-HG fit')
        ax.plot(cost, hgfun, 'g:', label='HG(<cos>)')
        ax.set_xlim(xr)
        if ylog:
            ax.set_yscale('log')
        ax.set_xlabel(r'$\cos\theta$')
        ax.set_ylabel('S11')
        ax.legend(fontsize=8, loc='best')

    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'mueller_matrix.dat')
