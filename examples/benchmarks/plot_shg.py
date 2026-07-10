#!/usr/bin/env python3
# Overlay the vendored-SEDust dust emission on the Camps et al. (2015) SHG
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
fig,ax=plt.subplots(figsize=(7,5))
Us=[('1e-02','1.E-02'),('1e+00','1.E+00'),('1e+02','1.E+02'),('1e+04','1.E+04')]
for Utag,key in Us:
    curves=[bench(c,Utag) for c in codes]; curves=[c for c in curves if c is not None]
    lamb=curves[0][0]; stack=np.array([c[1] for c in curves])
    lo,hi=stack.min(0),stack.max(0)
    ax.fill_between(lamb, lamb*lo, lamb*hi, alpha=0.3, color='0.6')
    of=glob.glob('sedust_zubko_Mathis_U_*.dat')
    of=[g for g in of if key in g.replace(' ','')][0]
    od=loadfix(of); ax.plot(od[:,0], od[:,0]*od[:,1], lw=1.3,
        label=r'$U=10^{'+str(int(np.log10(float(Utag))))+r'}$')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(1,2000); ax.set_xlabel(r'$\lambda\ (\mu{\rm m})$')
ax.set_ylabel(r'$\lambda\,\varepsilon_\lambda$ (per H, arb.)')
ax.set_ylim(1e-27,1e-16)
ax.legend(title=r'SEDust (lines) vs SHG codes (bands)', frameon=False, fontsize=9)
fig.tight_layout(); fig.savefig('shg_benchmark.pdf'); print('wrote shg_benchmark.pdf')
