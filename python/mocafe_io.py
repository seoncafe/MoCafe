#!/usr/bin/env python3
"""mocafe_io.py — format-agnostic MoCafe I/O.

Python counterpart of the Fortran ``iofile_mod`` facade introduced in
MoCafe v1.20.  Reads MoCafe output files in either FITS (.fits, .fits.gz)
or HDF5 (.h5, .hdf5) form and exposes a uniform structure, and converts
between formats by re-emitting the same logical content under the target
schema.

This module shares its core layout-mapping with LaRT's ``lart_io.py`` so
files round-trip cleanly between the two formats, but provides MoCafe-
specific affordances for the image-only output schema (no Spectrum
BinTable, no 3-D cubes), the ``_obs`` / ``_obs_tau`` / ``_stokes`` naming
pattern, and multi-observer runs (which write ``_001``, ``_002`` etc.).

Layout mapping (matches what MoCafe writes through ``iofile_mod``):

  FITS                                  HDF5
  -----------------------------------------------------------------
  Empty Primary HDU                     root group, attrs only
  Image Primary HDU (EXTNAME=X)         /X/data dataset + group attrs
  Image extension HDU (EXTNAME=X)       /X/data dataset + group attrs
  HDU header keyword K = v              attribute @K on the group

Public API
----------
``load_mocafe(path)``
    Read either a .fits[.gz] or .h5/.hdf5 file and return a
    :class:`MoCafeFile`.  Convenience accessors expose the standard
    image sections (Scattered, Direct, Direct0, Stokes_I/Q/U/V,
    TAU_dust) directly.

``find_mocafe_outputs(stem)``
    Locate the related files of a MoCafe run (the main peeloff,
    sightline tau, Stokes split, and any observer-suffixed files).

``convert(src, dst)``
    Convert between formats by extension.  ``foo_obs.fits.gz`` →
    ``foo_obs.h5`` and vice versa.

Command line
------------
    python mocafe_io.py info    <file>
    python mocafe_io.py peek    <file>     # image-stat summary
    python mocafe_io.py convert <src> <dst>
"""
from __future__ import annotations

import glob
import os
import sys
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional

import numpy as np


# ---------------------------------------------------------------------------
# Format detection
# ---------------------------------------------------------------------------

_HDF5_EXTS = ('.h5', '.hdf5')
_FITS_EXTS = ('.fits', '.fits.gz', '.fit')

# Order matters: HDF5 first because MoCafe v1.20 defaults to HDF5.
_MOCAFE_EXTS = ('.h5', '.hdf5', '.fits.gz', '.fits')

# The suffixes MoCafe attaches to par%base_name when writing files.
# Empty string covers the case where the user gave an explicit par%out_file.
_MOCAFE_SUFFIXES = ('_obs', '_obs_tau', '_stokes', '')


def _strip_mocafe_ext(name: str) -> str:
    """Return ``name`` with any recognised MoCafe extension stripped."""
    lower = name.lower()
    for ext in _MOCAFE_EXTS:
        if lower.endswith(ext):
            return name[: -len(ext)]
    return name


def detect_format(path: str) -> str:
    """Return ``'fits'`` or ``'hdf5'`` based on the filename extension.

    Falls back to a magic-byte sniff (``\x89HDF`` / ``SIMPLE``) if the
    extension is unrecognised; raises ``ValueError`` if both fail.
    """
    lower = path.lower()
    if lower.endswith(_HDF5_EXTS):
        return 'hdf5'
    for ext in _FITS_EXTS:
        if lower.endswith(ext):
            return 'fits'
    try:
        with open(path, 'rb') as fh:
            head = fh.read(8)
        if head[:4] == b'\x89HDF':
            return 'hdf5'
        if head[:6] == b'SIMPLE':
            return 'fits'
    except OSError:
        pass
    raise ValueError(f'Cannot determine MoCafe format for {path!r}')


def find_mocafe_file(stem: str, suffix: str = '') -> Optional[str]:
    """Return the first existing MoCafe file for ``<stem><suffix><ext>``.

    ``stem`` may be either a bare base name or a path that already
    carries an extension; the extension is stripped if present.
    ``suffix`` is an optional insertion (``'_obs'``, ``'_obs_tau'``,
    ``'_stokes'``, or observer ``'_001'``, …).  Returns ``None`` if
    no candidate exists.
    """
    base = _strip_mocafe_ext(stem)
    for ext in _MOCAFE_EXTS:
        candidate = base + suffix + ext
        if os.path.exists(candidate):
            return candidate
    return None


def find_mocafe_outputs(stem: str) -> Dict[str, List[str]]:
    """Locate the related files of a MoCafe run.

    Returns a dict with keys
        ``'obs'``     : the main peeloff output (Scattered + Direct + Direct0)
        ``'tau'``     : the sightline-tau image, if present
        ``'stokes'``  : the Stokes-I/Q/U/V split, if present
    Each value is a list — for a single-observer run it has at most one
    entry; for multi-observer runs it has one entry per observer
    (``_001``, ``_002``, …).
    """
    base = _strip_mocafe_ext(stem)
    found: Dict[str, List[str]] = {'obs': [], 'tau': [], 'stokes': []}

    def _push(key: str, suffix: str) -> None:
        # Single-observer file first.
        single = find_mocafe_file(base, suffix)
        if single is not None:
            found[key].append(single)
        # Files for each observer, _NNN form.
        for ext in _MOCAFE_EXTS:
            for path in sorted(glob.glob(f'{base}{suffix}_[0-9][0-9][0-9]{ext}')):
                if path not in found[key]:
                    found[key].append(path)

    _push('obs',    '_obs')
    _push('tau',    '_obs_tau')
    _push('stokes', '_stokes')
    return found


# ---------------------------------------------------------------------------
# In-memory representation
# ---------------------------------------------------------------------------

@dataclass
class Section:
    """One logical block — an image HDU or a BinTable HDU.

    ``name`` is the EXTNAME (or auto-generated ``section_NNN``).  For
    image sections ``data`` is the numpy array; for tables ``data`` is
    ``None`` and the columns live in ``columns``.  MoCafe writes only
    images today, but the table path is kept for forward compatibility
    and FITS round-trip with table-bearing files imported from elsewhere.
    """
    name: str
    data: Optional[np.ndarray] = None
    columns: Dict[str, np.ndarray] = field(default_factory=dict)
    attrs: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_table(self) -> bool:
        return self.data is None and len(self.columns) > 0

    def attr(self, name: str, default: Any = None) -> Any:
        """Case-insensitive attribute lookup."""
        for k, v in self.attrs.items():
            if k.lower() == name.lower():
                return v
        return default

    def col(self, name: str, default: Any = None) -> Any:
        """Case-insensitive column lookup."""
        for k, v in self.columns.items():
            if k.lower() == name.lower():
                return v
        return default


@dataclass
class MoCafeFile:
    path: str
    fmt: str                              # 'fits' or 'hdf5'
    sections: List[Section] = field(default_factory=list)
    root_attrs: Dict[str, Any] = field(default_factory=dict)

    # ---- generic section access ------------------------------------------------
    def section(self, name: str) -> Optional[Section]:
        """Case-insensitive section lookup; returns ``None`` if absent."""
        target = name.lower()
        for s in self.sections:
            if s.name.lower() == target:
                return s
        return None

    def __getitem__(self, key: str) -> Section:
        s = self.section(key)
        if s is None:
            raise KeyError(key)
        return s

    def __contains__(self, key: str) -> bool:
        return self.section(key) is not None

    # ---- MoCafe-specific convenience accessors --------------------------------
    @property
    def scattered(self) -> Optional[np.ndarray]:
        s = self.section('Scattered')
        return None if s is None else s.data

    @property
    def direct(self) -> Optional[np.ndarray]:
        s = self.section('Direct')
        return None if s is None else s.data

    @property
    def direct0(self) -> Optional[np.ndarray]:
        s = self.section('Direct0')
        return None if s is None else s.data

    @property
    def total(self) -> Optional[np.ndarray]:
        """Scattered + Direct, or ``None`` if either is missing.

        Axis layouts (numpy, reversed from the Fortran (nx,ny,...) order):
          plain        Scattered (ny,nx)            + Direct (ny,nx)
          (a,g)-scan   Scattered (ng,na,ny,nx)      + Direct (ny,nx)      [2-D]
          tau-scan     Scattered (nt,ng,na,ny,nx)   + Direct (nt,ny,nx)   [3-D]
        The lower-rank direct image is broadcast over the extra leading axes.
        """
        sc, di = self.scattered, self.direct
        if sc is None or di is None:
            return None
        if sc.ndim == 5 and di.ndim == 3:
            # (nt,ng,na,ny,nx) + (nt,ny,nx) -> broadcast direct over (ng,na)
            return sc + di[:, None, None, :, :]
        if sc.ndim == 4 and di.ndim == 2:
            return sc + di[None, None, :, :]
        return sc + di

    @property
    def is_ag_scan(self) -> bool:
        """True if this file holds an (albedo, asymmetry-factor) scan (4-D or 5-D)."""
        s = self.section('Scattered')
        return s is not None and s.data is not None and s.data.ndim in (4, 5)

    @property
    def is_tau_scan(self) -> bool:
        """True if this file holds a polychromatic (tau) scan: a 5-D (a,g,tau)
        Scattered cube, a 3-D tau-only Scattered cube (na = ng = 1), or a 3-D
        Direct cube."""
        s = self.section('Scattered')
        if s is not None and s.data is not None:
            if s.data.ndim == 5:
                return True
            if s.data.ndim == 3 and (s.attr('AG_NT') is not None or
                                     str(s.attr('CTYPE3', '')).upper() == 'TAU'):
                return True
        d = self.section('Direct')
        return d is not None and d.data is not None and d.data.ndim == 3

    def _ag_axes(self):
        """Return (albedo_axis, hgg_axis) as 1-D arrays, or (None, None).

        Works for both the 4-D (a,g) and 5-D (a,g,tau) Scattered cubes.
        Prefers the explicit AVALnnn / GVALnnn keyword lists (authoritative,
        recover non-uniform grids); falls back to the CRVAL3/CDELT3 (albedo)
        and CRVAL4/CDELT4 (hgg) WCS keywords.
        """
        s = self.section('Scattered')
        if s is None or s.data is None or s.data.ndim not in (4, 5):
            return None, None
        # numpy axes: 4-D (ng,na,ny,nx) -> na=shape[1],ng=shape[0];
        #             5-D (nt,ng,na,ny,nx) -> na=shape[2],ng=shape[1].
        if s.data.ndim == 5:
            ng, na = s.data.shape[1], s.data.shape[2]
        else:
            ng, na = s.data.shape[0], s.data.shape[1]
        na = int(s.attr('AG_NA', na))
        ng = int(s.attr('AG_NG', ng))
        a = [s.attr(f'AVAL{i:03d}') for i in range(1, na + 1)]
        g = [s.attr(f'GVAL{i:03d}') for i in range(1, ng + 1)]
        if all(v is not None for v in a) and all(v is not None for v in g):
            return np.asarray(a, float), np.asarray(g, float)
        a0, da = s.attr('CRVAL3'), s.attr('CDELT3')
        g0, dg = s.attr('CRVAL4'), s.attr('CDELT4')
        if None not in (a0, da, g0, dg):
            return (float(a0) + float(da) * np.arange(na),
                    float(g0) + float(dg) * np.arange(ng))
        return None, None

    def _tau_axis(self):
        """Return the 1-D tau (target taumax) axis of a tau scan, or None.

        The explicit TVALnnn list (authoritative) is written by
        write_tau_axis_keys on whichever HDU carries the tau axis -- the 5-D or
        3-D Scattered cube, or the 3-D Direct cube -- so read it from the first
        HDU that provides it.
        """
        for name in ('Scattered', 'Direct'):
            s = self.section(name)
            if s is None or s.data is None:
                continue
            nt = s.attr('AG_NT')
            if nt is None:
                continue
            nt = int(nt)
            t = [s.attr(f'TVAL{i:03d}') for i in range(1, nt + 1)]
            if all(v is not None for v in t):
                return np.asarray(t, float)
        return None

    @property
    def albedo_axis(self) -> Optional[np.ndarray]:
        """1-D array of albedo values along the albedo axis of a scan (else None)."""
        return self._ag_axes()[0]

    @property
    def hgg_axis(self) -> Optional[np.ndarray]:
        """1-D array of asymmetry (g) values along the g axis of a scan (else None)."""
        return self._ag_axes()[1]

    @property
    def tau_axis(self) -> Optional[np.ndarray]:
        """1-D array of tau (target taumax) values of a 5-D scan (else None)."""
        return self._tau_axis()

    @property
    def taumax_ref(self) -> Optional[float]:
        """The simulated reference taumax (s = 1 slice) of a tau scan, else None."""
        s = self.section('Scattered')
        if s is None:
            return None
        v = s.attr('AG_TAU0')
        return None if v is None else float(v)

    def scattered_at(self, a: Optional[float] = None,
                     g: Optional[float] = None,
                     tau: Optional[float] = None) -> Optional[np.ndarray]:
        """Return the 2-D Scattered image for the nearest (albedo, g, tau) on the grid.

        For a non-scan file the plain 2-D Scattered image is returned and the
        arguments are ignored.  Omitted axes default to the first grid point.
        Handles the 4-D (a,g) and 5-D (a,g,tau) cubes.
        """
        s = self.section('Scattered')
        if s is None or s.data is None:
            return None
        if s.data.ndim == 2:
            return s.data
        if s.data.ndim == 3:
            # tau-only cube (na = ng = 1): (nt, ny, nx)
            tax = self._tau_axis()
            it = 0 if (tau is None or tax is None) else int(np.argmin(np.abs(tax - tau)))
            return s.data[it, :, :]
        aax, gax = self._ag_axes()
        ia = 0 if (a is None or aax is None) else int(np.argmin(np.abs(aax - a)))
        ig = 0 if (g is None or gax is None) else int(np.argmin(np.abs(gax - g)))
        if s.data.ndim == 4:
            return s.data[ig, ia, :, :]
        tax = self._tau_axis()
        it = 0 if (tau is None or tax is None) else int(np.argmin(np.abs(tax - tau)))
        return s.data[it, ig, ia, :, :]

    def direct_at(self, tau: Optional[float] = None) -> Optional[np.ndarray]:
        """Return the 2-D Direct image for the nearest tau on the grid.

        For a non-tau-scan file the plain 2-D Direct image is returned and
        ``tau`` is ignored.  ``tau`` defaults to the first grid point.
        """
        di = self.direct
        if di is None:
            return None
        if di.ndim != 3:
            return di
        tax = self._tau_axis()
        it = 0 if (tau is None or tax is None) else int(np.argmin(np.abs(tax - tau)))
        return di[it, :, :]

    def stokes(self, component: str) -> Optional[np.ndarray]:
        """Return the named Stokes image (component in 'IQUV')."""
        c = component.upper()
        if c not in 'IQUV':
            raise ValueError(f'Stokes component must be I/Q/U/V, got {component!r}')
        s = self.section(f'Stokes_{c}')
        return None if s is None else s.data

    @property
    def tau_dust(self) -> Optional[np.ndarray]:
        s = self.section('TAU_dust')
        return None if s is None else s.data

    # ---- header inspection helpers --------------------------------------------
    def header(self, name: Optional[str] = None) -> Dict[str, Any]:
        """Attributes of one section (default: the first image section).

        Equivalent to reading an HDU header in FITS or a group's attrs in
        HDF5.  Use ``name`` to pick a non-default section.
        """
        if name is None:
            for s in self.sections:
                if s.data is not None:
                    return dict(s.attrs)
            return dict(self.root_attrs)
        s = self.section(name)
        return {} if s is None else dict(s.attrs)

    def info(self) -> str:
        lines = [f'MoCafeFile: {self.path}  ({self.fmt})']
        if self.root_attrs:
            lines.append('  root attrs: ' + ', '.join(sorted(self.root_attrs)))
        for s in self.sections:
            if s.is_table:
                cols = ', '.join(s.columns.keys())
                lines.append(f'  /{s.name} [table] cols=[{cols}]  attrs={len(s.attrs)}')
            elif s.data is not None:
                lines.append(f'  /{s.name} [image] shape={s.data.shape} '
                             f'dtype={s.data.dtype}  attrs={len(s.attrs)}')
            else:
                lines.append(f'  /{s.name} [empty]  attrs={len(s.attrs)}')
        return '\n'.join(lines)

    def peek(self) -> str:
        """One-line image stats for every image section.

        Useful for quick sanity checks: ``min``, ``max``, ``sum``, and the
        fraction of non-zero pixels.  Image sections only — tables are
        skipped.
        """
        lines = [f'MoCafeFile: {self.path}  ({self.fmt})']
        if self.is_ag_scan:
            aax, gax = self._ag_axes()
            if aax is not None:
                lines.append(f'  scan: albedo={np.array2string(aax, precision=2)}')
            if gax is not None:
                lines.append(f'        hgg   ={np.array2string(gax, precision=2)}')
        if self.is_tau_scan:
            tax = self._tau_axis()
            if tax is not None:
                lines.append(f'        tau   ={np.array2string(tax, precision=3)}'
                             f'   (ref taumax={self.taumax_ref})')
        for s in self.sections:
            if s.data is None:
                continue
            a = s.data
            nz_count = int(np.count_nonzero(a))
            lines.append(
                f'  /{s.name:<12}  shape={str(a.shape):<18}  '
                f'min={a.min():+.3e}  max={a.max():+.3e}  '
                f'sum={a.sum():+.3e}  nonzero={nz_count}/{a.size}'
            )
        return '\n'.join(lines)


# ---------------------------------------------------------------------------
# FITS reader
# ---------------------------------------------------------------------------

def _load_fits(path: str) -> MoCafeFile:
    from astropy.io import fits

    out = MoCafeFile(path=path, fmt='fits')
    with fits.open(path) as hdul:
        for i, hdu in enumerate(hdul):
            attrs = _fits_header_to_attrs(hdu.header)
            extname = attrs.pop('EXTNAME', None) or hdu.name or f'section_{i+1:03d}'

            if isinstance(hdu, fits.BinTableHDU):
                cols = {cn: np.array(hdu.data[cn]) for cn in hdu.columns.names}
                out.sections.append(Section(name=extname, columns=cols, attrs=attrs))
            elif hdu.data is None or hdu.data.size == 0:
                # An empty Primary HDU — fold its attrs into root and skip.
                if i == 0:
                    out.root_attrs.update(attrs)
                    continue
                out.sections.append(Section(name=extname, attrs=attrs))
            else:
                out.sections.append(
                    Section(name=extname, data=np.array(hdu.data), attrs=attrs))
    return out


def _fits_header_to_attrs(header) -> Dict[str, Any]:
    """Extract user-relevant FITS header keywords as a plain dict.

    Drops structural keywords (BITPIX, NAXIS*, EXTEND, ...) that the HDF5
    side reconstructs from dataset shape/dtype itself.
    """
    skip_exact = {
        'SIMPLE', 'BITPIX', 'NAXIS', 'EXTEND', 'BSCALE', 'BZERO',
        'XTENSION', 'PCOUNT', 'GCOUNT', 'TFIELDS', 'CHECKSUM', 'DATASUM',
    }
    skip_prefix = ('NAXIS', 'TFORM', 'TTYPE', 'TUNIT', 'TDIM',
                   'COMMENT', 'HISTORY')
    out: Dict[str, Any] = {}
    for card in header.cards:
        k = card.keyword
        if not k or k in skip_exact:
            continue
        if any(k.startswith(p) for p in skip_prefix):
            continue
        v = card.value
        if v is None:
            continue
        out[k] = v
    return out


# ---------------------------------------------------------------------------
# HDF5 reader
# ---------------------------------------------------------------------------

def _load_hdf5(path: str) -> MoCafeFile:
    import h5py

    out = MoCafeFile(path=path, fmt='hdf5')
    with h5py.File(path, 'r') as f:
        for k, v in f.attrs.items():
            out.root_attrs[k] = _h5_attr_to_py(v)

        for nm in _h5_children_in_order(f):
            obj = f[nm]
            if not isinstance(obj, h5py.Group):
                out.sections.append(Section(name=nm, data=np.array(obj), attrs={}))
                continue
            attrs = {k: _h5_attr_to_py(v) for k, v in obj.attrs.items()}
            children = list(obj.keys())
            if children == ['data']:
                out.sections.append(
                    Section(name=nm, data=np.array(obj['data']), attrs=attrs))
            else:
                cols = {cn: np.array(obj[cn]) for cn in children
                        if isinstance(obj[cn], h5py.Dataset)}
                out.sections.append(Section(name=nm, columns=cols, attrs=attrs))
    return out


def _h5_children_in_order(group) -> List[str]:
    """Best-effort enumerate of child names in creation order."""
    try:
        gcpl = group.id.get_create_plist()
        if gcpl.get_link_creation_order() & 0x1:  # H5P_CRT_ORDER_TRACKED
            n = group.id.get_num_objs()
            out = []
            for i in range(n):
                name = group.id.get_objname_by_idx(i)
                out.append(name.decode() if isinstance(name, bytes) else name)
            return out
    except Exception:
        pass
    return list(group.keys())


def _h5_attr_to_py(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        if value.shape == (1,):
            return _h5_attr_to_py(value[0])
        return value.tolist()
    if isinstance(value, (bytes, np.bytes_)):
        try:
            return value.decode().rstrip('\x00').rstrip()
        except UnicodeDecodeError:
            return value
    if isinstance(value, np.generic):
        return value.item()
    return value


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

def _write_fits(mf: MoCafeFile, path: str) -> None:
    from astropy.io import fits

    hdul = fits.HDUList()
    primary_filled = False
    primary_header = _attrs_to_fits_header(mf.root_attrs)

    for sec in mf.sections:
        hdr = _attrs_to_fits_header(sec.attrs)
        if sec.name and 'EXTNAME' not in hdr:
            hdr['EXTNAME'] = sec.name

        if sec.is_table:
            cols = [fits.Column(name=cn, array=arr,
                                format=_numpy_to_tform(arr))
                    for cn, arr in sec.columns.items()]
            if not primary_filled:
                hdul.append(fits.PrimaryHDU(header=primary_header))
                primary_filled = True
            hdul.append(fits.BinTableHDU.from_columns(cols, header=hdr, name=sec.name))
        elif sec.data is not None:
            if not primary_filled:
                merged = fits.Header()
                merged.update(primary_header)
                merged.update(hdr)
                hdul.append(fits.PrimaryHDU(data=sec.data, header=merged))
                primary_filled = True
            else:
                hdul.append(fits.ImageHDU(data=sec.data, header=hdr, name=sec.name))

    if not primary_filled:
        hdul.append(fits.PrimaryHDU(header=primary_header))

    if os.path.exists(path):
        os.remove(path)
    hdul.writeto(path)


def _attrs_to_fits_header(attrs: Dict[str, Any]):
    from astropy.io import fits
    h = fits.Header()
    for k, v in attrs.items():
        if v is None:
            continue
        try:
            if isinstance(v, (bytes, np.bytes_)):
                v = v.decode()
            # FITS keyword max 8 chars; preserve mixed case to match what
            # MoCafe writes (XMAX, nphotons, src_geom, ...).
            h[str(k)[:8]] = v
        except Exception:
            pass
    return h


def _numpy_to_tform(arr: np.ndarray) -> str:
    dt = arr.dtype
    if dt == np.float64:
        return 'D'
    if dt == np.float32:
        return 'E'
    if dt == np.int32:
        return 'J'
    if dt == np.int64:
        return 'K'
    if dt == np.int16:
        return 'I'
    if dt.kind in ('U', 'S'):
        n = max(1, arr.dtype.itemsize if dt.kind == 'S' else dt.itemsize // 4)
        return f'{n}A'
    return 'D'


def _write_hdf5(mf: MoCafeFile, path: str) -> None:
    import h5py

    if os.path.exists(path):
        os.remove(path)
    with h5py.File(path, 'w', libver='latest', track_order=True) as f:
        for k, v in mf.root_attrs.items():
            _h5_set_attr(f, k, v)
        for sec in mf.sections:
            g = f.create_group(sec.name, track_order=True)
            for k, v in sec.attrs.items():
                _h5_set_attr(g, k, v)
            if sec.is_table:
                for cn, arr in sec.columns.items():
                    _h5_create_dataset(g, cn, arr)
            elif sec.data is not None:
                _h5_create_dataset(g, 'data', sec.data)


def _h5_create_dataset(parent, name: str, arr: np.ndarray) -> None:
    kwargs: Dict[str, Any] = {}
    if arr.size > 4096:
        if arr.ndim == 1:
            chunks = (min(arr.shape[0], 4096),)
        else:
            chunks = tuple(min(s, 64) for s in arr.shape)
        kwargs.update(dict(chunks=chunks, compression='gzip', compression_opts=4))
    parent.create_dataset(name, data=arr, **kwargs)


def _h5_set_attr(obj, name: str, value: Any) -> None:
    if isinstance(value, str):
        obj.attrs.create(name, value)
    else:
        obj.attrs[name] = value


# ---------------------------------------------------------------------------
# Public reader / converter
# ---------------------------------------------------------------------------

def load_mocafe(path: str) -> MoCafeFile:
    """Read a MoCafe output file in either FITS or HDF5 format."""
    fmt = detect_format(path)
    if fmt == 'fits':
        return _load_fits(path)
    return _load_hdf5(path)


def convert(src: str, dst: str) -> None:
    """Convert ``src`` to ``dst``, picking direction from file extensions."""
    src_fmt = detect_format(src)
    dst_fmt = detect_format(dst)
    mf = load_mocafe(src)
    if dst_fmt == 'fits':
        _write_fits(mf, dst)
    else:
        _write_hdf5(mf, dst)
    print(f'Converted {src} [{src_fmt}] -> {dst} [{dst_fmt}]  '
          f'({len(mf.sections)} sections)')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _usage() -> str:
    return ('Usage:\n'
            '  python mocafe_io.py info    <file>\n'
            '  python mocafe_io.py peek    <file>\n'
            '  python mocafe_io.py convert <src> <dst>\n')


def _main(argv: List[str]) -> int:
    if len(argv) < 2:
        print(_usage(), file=sys.stderr)
        return 2
    cmd = argv[1]
    if cmd == 'info' and len(argv) == 3:
        print(load_mocafe(argv[2]).info())
        return 0
    if cmd == 'peek' and len(argv) == 3:
        print(load_mocafe(argv[2]).peek())
        return 0
    if cmd == 'convert' and len(argv) == 4:
        convert(argv[2], argv[3])
        return 0
    print(_usage(), file=sys.stderr)
    return 2


if __name__ == '__main__':
    sys.exit(_main(sys.argv))
