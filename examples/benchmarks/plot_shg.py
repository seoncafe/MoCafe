#!/usr/bin/env python3
# Overlay the SEDust dust emission on the Camps et al. (2015) SHG
# multi-code benchmark envelope for the Mathis-ISRF radiation fields.
import numpy as np, glob, os, re
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=os.path.expanduser('~/MoCafe/Grain/SHG_Benchmark/Results_FullSolution')
codes=['CRT','DIRTY','SKIRT','DustEM','TRADING','MCFOST','DARTRAY']
def loadfix(f):
    rows=[]
    for line in open(f):
        s=line.strip()
        if not s or s[0]=='#': continue
        v=[]
        for t in s.split():
            t2=re.sub(r'^([+-]?\d*\.?\d+)([+-]\d{2,3})$', r'\1e\2', t)
            try: v.append(float(t2))
            except: v.append(0.0)
        rows.append(v)
    return np.array(rows)
def bench(c,U):
    f=f'{B}/{c}_SHG_Mathis_U_{U}.dat'
    return (loadfix(f)[:,0], loadfix(f)[:,2:].sum(1)) if os.path.exists(f) else None
def area(l,e):
    # Scale of a spectrum, integrated the way cmp_shg.py normalizes: over
    # d ln(lambda), on the points that carry it.
    i=e>0
    return np.trapz(e[i],np.log(l[i]))
fig,ax=plt.subplots(figsize=(7,5))
Us=[('1e-02','1.E-02'),('1e+00','1.E+00'),('1e+02','1.E+02'),('1e+04','1.E+04')]
for Utag,key in Us:
    curves=[bench(c,Utag) for c in codes]; curves=[c for c in curves if c is not None]
    lamb=curves[0][0]
    of=glob.glob('sedust_zubko_Mathis_U_*.dat')
    of=[g for g in of if key in g.replace(' ','')][0]
    od=loadfix(of)
    # The published spectra are in the benchmark's own units, ~7 decades below
    # ours at every U, so overlaying them raw puts every band in one clump at
    # the bottom of the axes instead of on the curve it belongs to.  Normalize
    # each code the way cmp_shg.py does and put the envelope back on our own
    # scale, which is what makes the two comparable -- the benchmark constrains
    # the SHAPE; the level is set by the absorbed power, checked separately.
    A=area(od[:,0],od[:,1])
    stack=np.array([c[1]/area(lamb,c[1]) for c in curves])
    lo,hi=stack.min(0)*A,stack.max(0)*A
    ax.fill_between(lamb, lo, hi, alpha=0.3, color='0.6')
    ax.plot(od[:,0], od[:,1], lw=1.3,
        label=r'$U=10^{'+str(int(np.log10(float(Utag))))+r'}$')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(1,2000); ax.set_xlabel(r'$\lambda\ (\mu{\rm m})$')
ax.set_ylabel(r'$\lambda\,\varepsilon_\lambda$ (per H)')
ax.set_ylim(1e-29,1e-18)
ax.legend(title=r'SEDust (lines) vs SHG codes (bands)', frameon=False, fontsize=9)
fig.tight_layout(); fig.savefig('shg_benchmark.pdf'); print('wrote shg_benchmark.pdf')
