#!/usr/bin/env python3
"""Generic AMR octree builder + writer for MoCafe (dust-only).

Dust-only slim of LaRT's AMR_grid.py: builds an octree over a cubic box and
writes the generic AMR file MoCafe reads via ``par%grid_type='amr'``
(read_generic_amr.f90).  Only the geometry and the dust-relevant leaf fields
are kept --- there is no velocity, temperature, electron/ion density, or
Lyman-alpha emissivity.

Leaf fields:  nH (mandatory) + optional metallicity, xHI, ndust.
File schema:  binary table columns  x, y, z, level, nH [, metallicity, xHI,
ndust]  with header keywords  BOXLEN, ORIGINX/Y/Z, NLEAF, NAXIS2.

Octant convention matches the Fortran octree (octree_mod.f90 / amr_build_tree):
child center = parent center + (2*bit-1)*child_half, bit = ix + 2*iy + 4*iz.

See AMR_CLUMPS_PLAN.md Part A and AMR_CONVERTERS_PLAN.md.
"""
import numpy as np

_OPTIONAL_COLUMNS = ("metallicity", "xHI", "ndust")


# --------------------------------------------------------------------------
class Cell:
    __slots__ = ("cx", "cy", "cz", "h", "level", "children", "props")

    def __init__(self, cx, cy, cz, h, level):
        self.cx, self.cy, self.cz, self.h, self.level = cx, cy, cz, h, level
        self.children = None
        self.props = {}          # 'nH', 'metallicity', 'xHI', 'ndust'

    @property
    def is_leaf(self):
        return self.children is None

    def split(self):
        hc = 0.5 * self.h
        ch = []
        for b in range(8):
            ix, iy, iz = b & 1, (b >> 1) & 1, (b >> 2) & 1
            ch.append(Cell(self.cx + (2*ix - 1)*hc,
                           self.cy + (2*iy - 1)*hc,
                           self.cz + (2*iz - 1)*hc, hc, self.level + 1))
        self.children = ch


# --------------------------------------------------------------------------
class AMRGrid:
    def __init__(self, boxlen=2.0, origin=None):
        self.boxlen = float(boxlen)
        if origin is None:
            origin = (-0.5*boxlen, -0.5*boxlen, -0.5*boxlen)
        self.origin = tuple(float(o) for o in origin)
        cx = self.origin[0] + 0.5*boxlen
        cy = self.origin[1] + 0.5*boxlen
        cz = self.origin[2] + 0.5*boxlen
        self.root = Cell(cx, cy, cz, 0.5*boxlen, 0)
        self.metadata = {}

    # ---- traversal ----
    def leaves(self):
        out, stack = [], [self.root]
        while stack:
            c = stack.pop()
            if c.is_leaf:
                out.append(c)
            else:
                stack.extend(c.children)
        return out

    def leaves_at(self, level):
        return [c for c in self.leaves() if c.level == level]

    # ---- refinement ----
    def refine_uniform(self, level):
        """Refine every cell until all leaves reach `level`."""
        def rec(c):
            if c.level >= level:
                return
            if c.is_leaf:
                c.split()
            for ch in c.children:
                rec(ch)
        rec(self.root)

    def refine(self, criterion, level_max):
        """Split any leaf for which criterion(cell) is True, up to level_max,
        iterating level by level."""
        for _ in range(level_max + 1):
            todo = [c for c in self.leaves() if c.level < level_max and criterion(c)]
            if not todo:
                break
            for c in todo:
                c.split()

    def refine_by_density(self, dens_fn, threshold=0.3, level_max=8,
                          level_min=3, nprobe=2, floor=0.0):
        """Adaptive refinement on the density gradient.  At each level the
        sub-cell centers of all current leaves are sampled in ONE batched
        dens_fn call; a leaf splits when (max-min)/(max+min+floor) >= threshold.
        dens_fn(x, y, z) must accept and return numpy arrays."""
        if level_min > 0:
            self.refine_uniform(level_min)
        for lvl in range(level_min, level_max):
            leaves = self.leaves_at(lvl)
            if not leaves:
                continue
            pts = _subcenters(leaves, nprobe)            # (nleaf, nprobe^3, 3)
            shp = pts.shape
            vals = np.asarray(dens_fn(pts[..., 0].ravel(),
                                      pts[..., 1].ravel(),
                                      pts[..., 2].ravel())).reshape(shp[0], shp[1])
            vmax = vals.max(axis=1)
            vmin = vals.min(axis=1)
            grad = (vmax - vmin) / (vmax + vmin + floor + 1e-300)
            for c, g in zip(leaves, grad):
                if g >= threshold:
                    c.split()

    # convenience alias used by the Illustris converter (velocity args ignored)
    def refine_by_physics(self, dens_fn=None, threshold=0.3, level_max=8,
                          level_min=3, nprobe=2, floor=0.0, **_ignored):
        self.refine_by_density(dens_fn, threshold, level_max, level_min, nprobe, floor)

    def refine_to_resolution(self, size_fn, level_max, factor=1.0, level_min=0):
        """Refine a leaf while its width (2*h) exceeds factor * size_fn(center)
        (e.g. the local Voronoi cell size), up to level_max."""
        if level_min > 0:
            self.refine_uniform(level_min)
        for lvl in range(max(level_min, 0), level_max):
            leaves = self.leaves_at(lvl)
            if not leaves:
                continue
            cx = np.array([c.cx for c in leaves])
            cy = np.array([c.cy for c in leaves])
            cz = np.array([c.cz for c in leaves])
            target = np.asarray(size_fn(cx, cy, cz))
            for c, t in zip(leaves, target):
                if t > 0 and 2.0 * c.h > factor * t:
                    c.split()

    # ---- field assignment (callable fn(x,y,z) -> array, or scalar) ----
    def _set_field(self, name, fn):
        lv = self.leaves()
        xs = np.array([c.cx for c in lv])
        ys = np.array([c.cy for c in lv])
        zs = np.array([c.cz for c in lv])
        if callable(fn):
            vals = np.asarray(fn(xs, ys, zs), dtype=float)
        else:
            vals = np.full(len(lv), float(fn))
        for c, v in zip(lv, vals):
            c.props[name] = float(v)

    def set_density(self, fn):          self._set_field("nH", fn)
    def set_metallicity(self, fn):      self._set_field("metallicity", fn)
    def set_neutral_fraction(self, fn): self._set_field("xHI", fn)
    def set_dust_density(self, fn):     self._set_field("ndust", fn)

    # ---- summary ----
    def info(self):
        lv = self.leaves()
        levels = np.array([c.level for c in lv])
        nH = np.array([c.props.get("nH", 0.0) for c in lv])
        s = [f"AMRGrid: {len(lv)} leaves, levels {levels.min()}-{levels.max()}, "
             f"boxlen={self.boxlen:g}"]
        if nH.size:
            pos = nH[nH > 0]
            if pos.size:
                s.append(f"  nH (nonzero): min={pos.min():.3e} max={pos.max():.3e} "
                         f"fill={pos.size/nH.size:.3f}")
        for lev in range(levels.min(), levels.max() + 1):
            s.append(f"  level {lev}: {(levels == lev).sum()} leaves")
        return "\n".join(s)

    # ---- I/O ----
    def _leaf_arrays(self):
        lv = self.leaves()
        n = len(lv)
        x = np.fromiter((c.cx for c in lv), float, n)
        y = np.fromiter((c.cy for c in lv), float, n)
        z = np.fromiter((c.cz for c in lv), float, n)
        level = np.fromiter((c.level for c in lv), np.int32, n)
        nH = np.fromiter((c.props.get("nH", 0.0) for c in lv), float, n)
        opt = {}
        for name in _OPTIONAL_COLUMNS:
            if any(name in c.props for c in lv):
                opt[name] = np.fromiter((c.props.get(name, 0.0) for c in lv), float, n)
        return x, y, z, level, nH, opt

    def _keywords(self):
        return dict(BOXLEN=self.boxlen, ORIGINX=self.origin[0],
                    ORIGINY=self.origin[1], ORIGINZ=self.origin[2])

    def write(self, filename):
        """Write the generic AMR file.  Format from the extension:
        .h5/.hdf5 -> HDF5, .fits/.fits.gz -> FITS, else text."""
        x, y, z, level, nH, opt = self._leaf_arrays()
        keys = self._keywords()
        keys["NLEAF"] = x.size
        if filename.endswith((".h5", ".hdf5")):
            self._write_hdf5(filename, x, y, z, level, nH, opt, keys)
        elif filename.endswith((".fits", ".fits.gz")):
            self._write_fits(filename, x, y, z, level, nH, opt, keys)
        else:
            self._write_text(filename, x, y, z, level, nH, opt, keys)
        print(f"AMRGrid.write: {x.size} leaves -> {filename} "
              f"(cols: x,y,z,level,nH" + "".join(',' + k for k in opt) + ")")

    def _write_fits(self, fn, x, y, z, level, nH, opt, keys):
        from astropy.io import fits
        cols = [fits.Column(name="x",     format="1D", array=x),
                fits.Column(name="y",     format="1D", array=y),
                fits.Column(name="z",     format="1D", array=z),
                fits.Column(name="level", format="1J", array=level),
                fits.Column(name="nH",    format="1D", array=nH)]
        for name, arr in opt.items():
            cols.append(fits.Column(name=name, format="1D", array=arr))
        tab = fits.BinTableHDU.from_columns(cols)
        for k, v in keys.items():
            tab.header[k] = v
        fits.HDUList([fits.PrimaryHDU(), tab]).writeto(fn, overwrite=True)

    def _write_hdf5(self, fn, x, y, z, level, nH, opt, keys):
        import h5py
        with h5py.File(fn, "w") as f:
            g = f.create_group("AMR_GRID")
            g.create_dataset("x", data=x)
            g.create_dataset("y", data=y)
            g.create_dataset("z", data=z)
            g.create_dataset("level", data=level.astype("i4"))
            g.create_dataset("nH", data=nH)
            for name, arr in opt.items():
                g.create_dataset(name, data=arr)
            for k, v in keys.items():
                g.attrs[k] = v
            g.attrs["NAXIS2"] = x.size

    def _write_text(self, fn, x, y, z, level, nH, opt, keys):
        names = ["x", "y", "z", "level", "nH"] + list(opt)
        with open(fn, "w") as f:
            f.write(f"{x.size} {self.boxlen!r}\n")
            f.write("# " + " ".join(names) + "\n")
            for i in range(x.size):
                row = [f"{x[i]:.8e}", f"{y[i]:.8e}", f"{z[i]:.8e}",
                       f"{level[i]:d}", f"{nH[i]:.8e}"]
                row += [f"{opt[k][i]:.8e}" for k in opt]
                f.write(" ".join(row) + "\n")

    @classmethod
    def read(cls, filename):
        """Read a generic AMR file back into an AMRGrid (for round-trip tests)."""
        if filename.endswith((".h5", ".hdf5")):
            import h5py
            with h5py.File(filename, "r") as f:
                g = f[list(f.keys())[0]]
                cols = {k: g[k][...] for k in g}
                keys = {k: g.attrs[k] for k in g.attrs}
        else:
            from astropy.io import fits
            with fits.open(filename) as h:
                d = h[1].data
                cols = {nm: np.array(d[nm]) for nm in d.columns.names}
                keys = {k: h[1].header[k] for k in h[1].header
                        if k in ("BOXLEN", "ORIGINX", "ORIGINY", "ORIGINZ")}
        boxlen = float(keys.get("BOXLEN", 2.0))
        origin = (float(keys.get("ORIGINX", -0.5*boxlen)),
                  float(keys.get("ORIGINY", -0.5*boxlen)),
                  float(keys.get("ORIGINZ", -0.5*boxlen)))
        grid = cls(boxlen=boxlen, origin=origin)
        nm = "nH" if "nH" in cols else ("dens" if "dens" in cols else "gasDen")
        x, y, z, lvl, nH = (cols["x"], cols["y"], cols["z"],
                            cols["level"].astype(int), cols[nm])
        for i in range(x.size):
            c = grid._insert(float(x[i]), float(y[i]), float(z[i]), int(lvl[i]))
            c.props["nH"] = float(nH[i])
            for k in _OPTIONAL_COLUMNS:
                if k in cols:
                    c.props[k] = float(cols[k][i])
        return grid

    def _insert(self, x, y, z, target_level):
        c = self.root
        for _ in range(target_level):
            if c.is_leaf:
                c.split()
            b = (1 if x >= c.cx else 0) + 2*(1 if y >= c.cy else 0) + 4*(1 if z >= c.cz else 0)
            c = c.children[b]
        return c


# --------------------------------------------------------------------------
def write_leaves(filename, x, y, z, level, nH, boxlen, origin,
                 metallicity=None, xHI=None, ndust=None):
    """Write a generic AMR file directly from flat leaf arrays (used by the
    RAMSES converter, which has native leaf positions + levels).  Same schema
    as AMRGrid.write: columns x,y,z,level,nH [,metallicity,xHI,ndust]; header
    keywords BOXLEN, ORIGINX/Y/Z, NLEAF, NAXIS2."""
    x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
    level = np.asarray(level, np.int32); nH = np.asarray(nH, float)
    opt = {}
    for name, arr in (("metallicity", metallicity), ("xHI", xHI), ("ndust", ndust)):
        if arr is not None:
            opt[name] = np.asarray(arr, float)
    keys = dict(BOXLEN=float(boxlen), ORIGINX=float(origin[0]),
                ORIGINY=float(origin[1]), ORIGINZ=float(origin[2]), NLEAF=x.size)
    if filename.endswith((".h5", ".hdf5")):
        import h5py
        with h5py.File(filename, "w") as f:
            g = f.create_group("AMR_GRID")
            g.create_dataset("x", data=x);  g.create_dataset("y", data=y)
            g.create_dataset("z", data=z);  g.create_dataset("level", data=level)
            g.create_dataset("nH", data=nH)
            for nm, arr in opt.items():
                g.create_dataset(nm, data=arr)
            for k, v in keys.items():
                g.attrs[k] = v
            g.attrs["NAXIS2"] = x.size
    elif filename.endswith((".fits", ".fits.gz")):
        from astropy.io import fits
        cols = [fits.Column(name="x", format="1D", array=x),
                fits.Column(name="y", format="1D", array=y),
                fits.Column(name="z", format="1D", array=z),
                fits.Column(name="level", format="1J", array=level),
                fits.Column(name="nH", format="1D", array=nH)]
        for nm, arr in opt.items():
            cols.append(fits.Column(name=nm, format="1D", array=arr))
        tab = fits.BinTableHDU.from_columns(cols)
        for k, v in keys.items():
            tab.header[k] = v
        fits.HDUList([fits.PrimaryHDU(), tab]).writeto(filename, overwrite=True)
    else:
        names = ["x", "y", "z", "level", "nH"] + list(opt)
        with open(filename, "w") as f:
            f.write(f"{x.size} {float(boxlen)!r}\n")
            f.write("# " + " ".join(names) + "\n")
            for i in range(x.size):
                row = [f"{x[i]:.8e}", f"{y[i]:.8e}", f"{z[i]:.8e}",
                       f"{level[i]:d}", f"{nH[i]:.8e}"]
                row += [f"{opt[k][i]:.8e}" for k in opt]
                f.write(" ".join(row) + "\n")
    print(f"write_leaves: {x.size} leaves -> {filename} "
          f"(cols: x,y,z,level,nH" + "".join(',' + k for k in opt) + ")")


def _subcenters(leaves, nprobe):
    """(nleaf, nprobe^3, 3) array of sub-cell-center sample points."""
    n = len(leaves)
    cx = np.array([c.cx for c in leaves])[:, None]
    cy = np.array([c.cy for c in leaves])[:, None]
    cz = np.array([c.cz for c in leaves])[:, None]
    h = np.array([c.h for c in leaves])[:, None]
    # offsets in [-h, h] at nprobe points per axis
    if nprobe == 1:
        frac = np.array([0.0])
    else:
        frac = (np.arange(nprobe) + 0.5) / nprobe * 2.0 - 1.0   # in (-1,1)
    ox, oy, oz = np.meshgrid(frac, frac, frac, indexing="ij")
    ox = ox.ravel()[None, :]; oy = oy.ravel()[None, :]; oz = oz.ravel()[None, :]
    px = cx + h * ox
    py = cy + h * oy
    pz = cz + h * oz
    return np.stack([px, py, pz], axis=-1)
