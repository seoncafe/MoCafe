#!/usr/bin/env python3
"""Plot a MoCafe HEALPix all-sky map (<base>_allsky.h5) with healpy.
Usage: python plot_healpix.py mw_allsky_allsky.h5 [lambda_um ...]"""
import sys, h5py, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import healpy as hp

fn = sys.argv[1] if len(sys.argv) > 1 else 'mw_allsky_allsky.h5'
lams = [float(x) for x in sys.argv[2:]] or [0.55, 100.0]
f = h5py.File(fn, 'r')
sky = f['AllSky/data'][:]                 # (nlambda, npix)
w   = f['Wavelength/data'][:]
npix = int(f['AllSky'].attrs['NPIX'])
nl = len(w)
ax = [i for i, s in enumerate(sky.shape) if s == nl][0]
# MoCafe uses 1-based RING; healpy uses 0-based RING with the same ordering,
# so array index i (= pixel i+1) is exactly healpy pixel i -> no offset.
fig = plt.figure(figsize=(7, 3*len(lams)))
for n, lam in enumerate(lams):
    i = int(np.argmin(np.abs(w - lam)))
    m = np.take(sky, i, axis=ax).astype(float).flatten()
    mm = np.where(m > 0, np.log10(m), hp.UNSEEN)
    hp.mollview(mm, fig=fig.number, sub=(len(lams), 1, n+1), nest=False,
                title=r'$\lambda=%.3g\,\mu$m' % lam, cmap='inferno',
                unit=r'$\log_{10} I$')
    hp.graticule()
out = fn.replace('.h5', '.pdf').replace('_allsky', '_healpix')
fig.savefig(out); print('wrote', out)
