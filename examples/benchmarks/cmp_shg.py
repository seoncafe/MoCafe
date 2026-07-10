import numpy as np, glob, os, re
B=os.path.expanduser('~/MoCafe/Grain/SHG_Benchmark/Results_FullSolution')
codes=['CRT','DIRTY','SKIRT','DustEM','TRADING','MCFOST','DARTRAY']
def loadfix(f):
    rows=[]
    for line in open(f):
        line=line.strip()
        if not line or line[0]=='#': continue
        vals=[]
        for t in line.split():
            t2=re.sub(r'^([+-]?\d*\.?\d+)([+-]\d{2,3})$', r'\1e\2', t)
            try: vals.append(float(t2))
            except: vals.append(0.0)
        rows.append(vals)
    return np.array(rows)
def bench(code,U):
    f=f'{B}/{code}_SHG_Mathis_U_{U}.dat'
    if not os.path.exists(f): return None
    d=loadfix(f); return d[:,0], d[:,2:].sum(axis=1)
our={}
for g in glob.glob('sedust_zubko_Mathis_U_*.dat'):
    d=loadfix(g); our[g.replace(' ','')]=d
def norm(l,e):
    i=e>0; A=np.trapz(e[i],np.log(l[i])); return e/A if A>0 else e
for Utag,key in [('1e-02','1.E-02'),('1e+00','1.E+00'),('1e+02','1.E+02'),('1e+04','1.E+04')]:
    of=[g for g in our if key in g]
    if not of: continue
    od=our[of[0]]; ol,oe=od[:,0],od[:,1]
    curves=[bench(c,Utag) for c in codes]; curves=[c for c in curves if c is not None]
    if not curves: continue
    lamb=curves[0][0]; oe_i=np.interp(lamb,ol,oe)
    stack=np.array([norm(lamb,c[1]) for c in curves])
    med=np.median(stack,0); lo=stack.min(0); hi=stack.max(0); on=norm(lamb,oe_i)
    mask=(lamb>3)&(lamb<1000)&(med>0.02*med.max())
    rel=np.abs(on[mask]-med[mask])/med[mask]
    within=((on[mask]>=lo[mask])&(on[mask]<=hi[mask])).mean()
    ipk=np.argmax((lamb*med)[mask]); ipo=np.argmax((lamb*on)[mask])
    print(f'U={Utag}: median rel.diff vs code-median={np.median(rel):.1%}, within multi-code envelope={within:.0%}, '
          f'peak ours {lamb[mask][ipo]:.0f}um / codes {lamb[mask][ipk]:.0f}um ({len(curves)} codes)')
