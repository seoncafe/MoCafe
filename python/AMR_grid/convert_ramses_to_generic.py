#!/usr/bin/env python3
"""Convert a RAMSES output to a generic AMR dust file for MoCafe (dust-only).

Reads the native RAMSES octree (amr_* / hydro_* Fortran-unformatted files,
parsed directly --- no yt) and extracts, per leaf cell:

    x, y, z, level, nH   (mandatory)
    metallicity          (optional; from a hydro var index or a uniform value)
    xHI                  (optional; from a hydro var index)
    ndust                (optional; --emit-ndust, Laursen+09 from nH,xHI,Z)

and writes the generic AMR file MoCafe reads with par%grid_type='amr'.

DUST-ONLY: gas velocity and temperature are NOT read or written (the
velocity/energy/pressure hydro variables are skipped).  This is the dust-only
slim of LaRT's convert_ramses_to_generic.py; the Fortran-record reader and the
AMR/hydro traversal are ported faithfully from it.  See AMR_CONVERTERS_PLAN.md.

Usage
-----
    python convert_ramses_to_generic.py output_00042 -o sim.fits \
        --output-unit kpc --metallicity-index 6 --emit-ndust
"""
import argparse
import re
import struct
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from AMR_grid import write_leaves

# unit conversions (cm)
KPC_TO_CM  = 3.0856775814913673e21
PC_TO_CM   = 3.0856775814913673e18
AU_TO_CM   = 1.495978707e13
MASS_H_CGS = 1.6726e-24


def laursen09_ndust(nH, xHI, Z, Z_ref=0.0134, f_ion=0.01):
    nHI, nHII = nH * xHI, nH * (1.0 - xHI)
    return (Z / max(Z_ref, 1e-30)) * (nHI + f_ion * nHII)


@dataclass
class RamsesInfo:
    repository: Path; snapnum: int; ncpu: int; ndim: int
    nx: int; ny: int; nz: int; nlevelmax: int
    boxlen_code: float; unit_l: float; unit_d: float; unit_t: float; gamma: float

    @property
    def boxlen_cm(self):
        return self.boxlen_code * self.unit_l


class FortranRecordReader:
    """Sequential Fortran unformatted reader (endianness auto-detected)."""
    def __init__(self, filename, endian=None):
        self.filename = Path(filename)
        self.handle = self.filename.open("rb")
        self.endian = endian or self._detect_endian()

    def close(self):           self.handle.close()
    def __enter__(self):       return self
    def __exit__(self, *a):    self.close()

    def _detect_endian(self):
        pos = self.handle.tell()
        lead = self.handle.read(4)
        if len(lead) != 4:
            raise EOFError(f"cannot read first record marker from {self.filename}")
        payload = self.handle.read(8)
        self.handle.seek(pos)
        for endian in ("<", ">"):
            recl = struct.unpack(f"{endian}i", lead)[0]
            if recl <= 0 or recl > 1024:
                continue
            if len(payload) < recl + 4:
                continue
            if struct.unpack(f"{endian}i", payload[recl:recl+4])[0] == recl:
                return endian
        raise ValueError(f"cannot determine endianness of {self.filename}")

    def read_record(self):
        lead = self.handle.read(4)
        if not lead or len(lead) != 4:
            raise EOFError(f"truncated record marker in {self.filename}")
        recl = struct.unpack(f"{self.endian}i", lead)[0]
        payload = self.handle.read(recl)
        trail = self.handle.read(4)
        if len(payload) != recl or len(trail) != 4 or \
           struct.unpack(f"{self.endian}i", trail)[0] != recl:
            raise EOFError(f"truncated/mismatched record in {self.filename}")
        return payload

    def skip(self, nrec=1):
        for _ in range(nrec):
            self.read_record()

    def read_array(self, dtype):
        dtype = np.dtype(dtype).newbyteorder(self.endian)
        payload = self.read_record()
        if len(payload) % dtype.itemsize != 0:
            raise ValueError(f"record/dtype size mismatch in {self.filename}")
        return np.frombuffer(payload, dtype=dtype).copy()

    def read_scalar(self, dtype):
        a = self.read_array(dtype)
        return a[0]


def normalize_repository_and_snapnum(path, snapnum):
    p = Path(path).expanduser().resolve()
    m = re.fullmatch(r"output_(\d{5})", p.name)
    if m:
        inferred = int(m.group(1))
        snapnum = inferred if snapnum is None else snapnum
        p = p.parent
    if snapnum is None:
        raise ValueError("specify --snapnum, or pass a path ending in output_00042")
    return p, snapnum


def _f(repo, snap, name):
    return repo / f"output_{snap:05d}" / name
def info_filename(repo, snap):       return _f(repo, snap, f"info_{snap:05d}.txt")
def amr_filename(repo, snap, icpu):  return _f(repo, snap, f"amr_{snap:05d}.out{icpu:05d}")
def hydro_filename(repo, snap, icpu):return _f(repo, snap, f"hydro_{snap:05d}.out{icpu:05d}")
def descriptor_filename(repo, snap): return _f(repo, snap, "hydro_file_descriptor.txt")


def read_info(repo, snap):
    values = {}
    with info_filename(repo, snap).open("r") as h:
        for raw in h:
            if "=" not in raw:
                continue
            k, v = raw.split("=", 1)
            try:
                values[k.strip()] = float(v.strip().split()[0].replace("d", "e").replace("D", "e"))
            except ValueError:
                pass
    return values


def read_hydro_descriptor(repo, snap):
    fn = descriptor_filename(repo, snap)
    mapping = {}
    if not fn.exists():
        return mapping
    for raw in fn.read_text().splitlines():
        m = re.match(r"\s*variable\s*#\s*(\d+)\s*:\s*(.+?)\s*$", raw, re.IGNORECASE)
        if m:
            mapping[int(m.group(1))] = m.group(2).strip().lower()
    return mapping


def read_basic_amr_header(filename):
    with FortranRecordReader(filename) as r:
        ncpu = int(r.read_scalar(np.int32)); ndim = int(r.read_scalar(np.int32))
        nx, ny, nz = r.read_array(np.int32)
        nlevelmax = int(r.read_scalar(np.int32)); r.read_scalar(np.int32)
        nboundary = int(r.read_scalar(np.int32)); r.read_scalar(np.int32)
        boxlen = float(r.read_scalar(np.float64))
    return (ncpu, ndim, int(nx), int(ny), int(nz), nlevelmax, nboundary, boxlen)


def build_ramses_info(repo, snap):
    values = read_info(repo, snap)
    ncpu_h, ndim, nx, ny, nz, nlevelmax, nboundary, boxlen_code = \
        read_basic_amr_header(amr_filename(repo, snap, 1))
    if ndim != 3:
        raise ValueError(f"only 3D RAMSES is supported (ndim={ndim})")
    return RamsesInfo(repo, snap, int(values.get("ncpu", ncpu_h)), ndim, nx, ny, nz,
                      nlevelmax, float(values.get("boxlen", boxlen_code)),
                      float(values.get("unit_l", 1.0)), float(values.get("unit_d", 1.0)),
                      float(values.get("unit_t", 1.0)), float(values.get("gamma", 5.0/3.0)))


def read_amr_layout(reader, ncpu, nlevelmax, nboundary):
    reader.skip(13)
    ngridlevel = reader.read_array(np.int32).reshape((ncpu, nlevelmax), order="F")
    ngridfile = np.zeros((ncpu + nboundary, nlevelmax), dtype=np.int32)
    ngridfile[:ncpu, :] = ngridlevel
    reader.skip(1)
    if nboundary > 0:
        reader.skip(2)
        ngridbound = reader.read_array(np.int32).reshape((nboundary, nlevelmax), order="F")
        ngridfile[ncpu:, :] = ngridbound
    reader.skip(6)
    return ngridfile


def read_hydro_header(reader):
    reader.skip(1)
    nvar = int(reader.read_scalar(np.int32))
    reader.skip(3)
    gamma = float(reader.read_scalar(np.float64))
    return nvar, gamma


def oct_offsets(ndim, ilevel):
    twotondim = 2 ** ndim
    dx = 0.5 ** ilevel
    xc = np.zeros((twotondim, ndim))
    for ind in range(twotondim):
        if ndim >= 3:
            xc[ind, 2] = ind // 4
        if ndim >= 2:
            xc[ind, 1] = (ind // 2) % 2
        xc[ind, 0] = ind % 2
        xc[ind, :] = (xc[ind, :] - 0.5) * dx
    return xc


def _unit_factor(unit):
    return {"cm": 1.0, "kpc": KPC_TO_CM, "pc": PC_TO_CM, "au": AU_TO_CM}[unit.lower()]


def convert_ramses_snapshot(repo, snap, hydro_precision=8, density_index=None,
                            metallicity_index=None, ionization_index=None):
    """Extract dust-only leaf data from a RAMSES output.  Returns
    (info, dict of arrays: x_cm,y_cm,z_cm,level,nH[,metallicity,xHI])."""
    info = build_ramses_info(repo, snap)
    names = read_hydro_descriptor(repo, snap)
    dens_idx = density_index or next((i for i, n in names.items() if n == "density"), 1)
    metal_idx = metallicity_index or next((i for i, n in names.items() if n == "metallicity"), None)
    ion_idx = ionization_index

    xb = np.array([info.nx // 2, info.ny // 2, info.nz // 2], float)
    twotondim = 2 ** info.ndim
    out = {k: [] for k in ("x", "y", "z", "level", "nH", "metallicity", "xHI")}

    for icpu in range(1, info.ncpu + 1):
        with FortranRecordReader(amr_filename(repo, snap, icpu)) as ar, \
             FortranRecordReader(hydro_filename(repo, snap, icpu)) as hr:
            ncpu_f = int(ar.read_scalar(np.int32)); ndim = int(ar.read_scalar(np.int32))
            ar.read_array(np.int32)  # nx,ny,nz
            nlevelmax = int(ar.read_scalar(np.int32)); ar.skip(1)
            nboundary = int(ar.read_scalar(np.int32)); ar.skip(2)
            ngridfile = read_amr_layout(ar, ncpu_f, nlevelmax, nboundary)
            nvar, _ = read_hydro_header(hr)
            need = max([dens_idx] + [i for i in (metal_idx, ion_idx) if i])
            if need > nvar:
                raise ValueError(f"hydro nvar={nvar} < requested index {need}")

            for ilevel in range(1, nlevelmax + 1):
                ngrida = int(ngridfile[icpu - 1, ilevel - 1])
                xc = oct_offsets(info.ndim, ilevel)
                if ngrida > 0:
                    xg = np.empty((ngrida, info.ndim)); son = np.empty((ngrida, twotondim), np.int32)
                    var = np.empty((ngrida, twotondim, nvar))
                else:
                    xg = son = var = None
                for j in range(1, ncpu_f + nboundary + 1):
                    ngrid_dom = int(ngridfile[j - 1, ilevel - 1])
                    if ngrid_dom > 0:
                        ar.skip(3)
                        for idim in range(info.ndim):
                            if j == icpu:
                                xg[:, idim] = ar.read_array(np.float64)
                            else:
                                ar.skip(1)
                        ar.skip(1 + 2 * info.ndim)
                        for ind in range(twotondim):
                            if j == icpu:
                                son[:, ind] = ar.read_array(np.int32)
                            else:
                                ar.skip(1)
                        ar.skip(2 * twotondim)
                    hr.skip(2)
                    if ngrid_dom > 0:
                        for ind in range(twotondim):
                            for ivar in range(nvar):
                                if j == icpu:
                                    if hydro_precision == 4:
                                        var[:, ind, ivar] = hr.read_array(np.float32).astype(np.float64)
                                    else:
                                        var[:, ind, ivar] = hr.read_array(np.float64)
                                else:
                                    hr.skip(1)
                if ngrida <= 0:
                    continue
                for ind in range(twotondim):
                    mask = son[:, ind] == 0
                    if not np.any(mask):
                        continue
                    x_code = (xg[:, 0] + xc[ind, 0] - xb[0]) / info.nx + 0.5
                    y_code = (xg[:, 1] + xc[ind, 1] - xb[1]) / info.ny + 0.5
                    z_code = (xg[:, 2] + xc[ind, 2] - xb[2]) / info.nz + 0.5
                    nH = var[:, ind, dens_idx - 1] * info.unit_d / MASS_H_CGS
                    out["x"].append(x_code[mask] * info.boxlen_cm)
                    out["y"].append(y_code[mask] * info.boxlen_cm)
                    out["z"].append(z_code[mask] * info.boxlen_cm)
                    out["level"].append(np.full(np.count_nonzero(mask), ilevel, np.int32))
                    out["nH"].append(nH[mask])
                    if metal_idx:
                        out["metallicity"].append(var[:, ind, metal_idx - 1][mask])
                    if ion_idx:
                        out["xHI"].append(var[:, ind, ion_idx - 1][mask])

    res = {"x_cm": np.concatenate(out["x"]) if out["x"] else np.array([]),
           "y_cm": np.concatenate(out["y"]) if out["y"] else np.array([]),
           "z_cm": np.concatenate(out["z"]) if out["z"] else np.array([]),
           "level": np.concatenate(out["level"]) if out["level"] else np.array([], np.int32),
           "nH": np.concatenate(out["nH"]) if out["nH"] else np.array([])}
    if metal_idx:
        res["metallicity"] = np.concatenate(out["metallicity"]) if out["metallicity"] else None
    if ion_idx:
        res["xHI"] = np.concatenate(out["xHI"]) if out["xHI"] else None
    return info, res


def main(argv=None):
    ap = argparse.ArgumentParser(description="RAMSES -> generic AMR (dust-only).")
    ap.add_argument("repository", help="RAMSES dir or .../output_00042")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("-s", "--snapnum", type=int, default=None)
    ap.add_argument("--output-unit", default="kpc", choices=["cm", "kpc", "pc", "au"])
    ap.add_argument("--hydro-precision", type=int, default=8, choices=[4, 8])
    ap.add_argument("--density-index", type=int, default=None,
                    help="1-based hydro var index for density (else from descriptor)")
    ap.add_argument("--metallicity-index", type=int, default=None,
                    help="1-based hydro var index for the metal mass fraction Z")
    ap.add_argument("--metallicity", type=float, default=None,
                    help="uniform Z when no metal column is read")
    ap.add_argument("--ionization-index", type=int, default=None,
                    help="1-based hydro var index for the neutral fraction xHI")
    ap.add_argument("--emit-ndust", action="store_true",
                    help="precompute ndust = laursen09(nH,xHI,Z)")
    ap.add_argument("--Z-ref", type=float, default=0.0134)
    ap.add_argument("--f-ion", type=float, default=0.01)
    args = ap.parse_args(argv)

    repo, snap = normalize_repository_and_snapnum(args.repository, args.snapnum)
    info, d = convert_ramses_snapshot(repo, snap, args.hydro_precision,
                                      args.density_index, args.metallicity_index,
                                      args.ionization_index)
    fac = _unit_factor(args.output_unit)
    x = d["x_cm"] / fac;  y = d["y_cm"] / fac;  z = d["z_cm"] / fac
    boxlen = info.boxlen_cm / fac
    nH = d["nH"]
    Z = d.get("metallicity")
    if Z is None and args.metallicity is not None:
        Z = np.full(nH.size, args.metallicity)
    xHI = d.get("xHI")
    ndust = None
    if args.emit_ndust:
        if Z is None:
            ap.error("--emit-ndust needs --metallicity-index or --metallicity")
        xh = xHI if xHI is not None else np.ones(nH.size)
        ndust = laursen09_ndust(nH, xh, Z, args.Z_ref, args.f_ion)

    print(f"RAMSES output_{snap:05d}: {nH.size} leaves, boxlen={boxlen:.4g} {args.output_unit}, "
          f"levels {d['level'].min()}-{d['level'].max()}, nH=[{nH.min():.3e},{nH.max():.3e}]")
    # RAMSES native box is [0, boxlen]; MoCafe recenters via ORIGIN.
    write_leaves(args.output, x, y, z, d["level"], nH, boxlen,
                 origin=(0.0, 0.0, 0.0), metallicity=Z, xHI=xHI, ndust=ndust)
    print("Hint: MoCafe namelist -> par%grid_type='amr', par%amr_file='"
          + args.output + "', par%distance_unit='" + args.output_unit + "',")
    print("      par%dust_model='global_dgr' | 'laursen09' (needs metallicity+xHI) "
          "| 'from_file' (--emit-ndust).")


if __name__ == "__main__":
    main()
