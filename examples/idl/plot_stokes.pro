pro plot_stokes,fname

if n_elements(fname) eq 0 then fname = 'test1'
fname =fname+'_stokes.fits.gz'
inten = mrdfits(fname,0,/silen)
q     = mrdfits(fname,1,/silen)
u     = mrdfits(fname,2,/silen)
v     = mrdfits(fname,3,/silen)
pol   = sqrt(q*q + u*u)

pol_perc = pol/inten
theta    = 0.5d0*atan(u,q)/!dtor
w_nan    = where(inten le 0d0, nw)
if nw gt 0 then begin
   pol_perc[w_nan] = !values.f_nan
   theta[w_nan]    = !values.f_nan
endif

;mwrfits,theta,'theta.fits.gz',/create

sz   = size(inten)
xr   = [-0.5,0.5]*sz[1]
yr   = [-0.5,0.5]*sz[1]

!p.multi=[0,3,2]
!p.charsize=2.0
!x.title = 'X'
!y.title = 'Y'
margin=0.05
loadct,39
w = where(inten ne 0d0 and q ne 0d0 and u ne 0d0 and v ne 0d0, compl=wc, ncompl=nwc)
if (nwc gt 0) then pol_perc[wc] = 0d0
imdisp,asinh(inten/median(inten[w])),/er,/ax,tit='I',margin=margin,xr=xr,yr=yr
imdisp,asinh(q/abs(median(q[w]))),/er,/ax,tit='Q',margin=margin,xr=xr,yr=yr
imdisp,asinh(u/abs(median(u[w]))),/er,/ax,tit='U',margin=margin,xr=xr,yr=yr
imdisp,asinh(v/abs(median(v[w]))),/er,/ax,tit='V',margin=margin,xr=xr,yr=yr

;imdisp,asinh(pol/median(pol)),/er,/ax,tit=textoidl('(Q^2+U^2)^{1/2}'),margin=margin,xr=xr,yr=yr
print, median(pol_perc[w])
imdisp,asinh(pol_perc/median(pol_perc[w])),/er,/ax,tit='Polarization (%)',margin=margin,xr=xr,yr=yr
;imdisp,theta,/er,/ax

dn   = 8L
nn   = max([sz[1],sz[2]])
n    = nn/dn
x    = dblarr(n,n)
y    = dblarr(n,n)
xvec = dblarr(n,n)
yvec = dblarr(n,n)
scal = 1d0/median(pol_perc[w])*(nn/(dn*2d0))

for i=0L,n-1 do begin
for j=0L,n-1 do begin
   i1 = i*dn + dn/2L-1
   j1 = j*dn + dn/2L-1
   x[i,*] = i1
   y[*,j] = j1
   xvec[i,j] = -pol_perc[i1,j1]*scal*sin(theta[i1,j1]*!dtor)
   yvec[i,j] =  pol_perc[i1,j1]*scal*cos(theta[i1,j1]*!dtor)
endfor
endfor

plotpol,x,y,xvec,yvec,xr=[0d0,sz[1]],yr=[0d0,sz[2]],/xst,/yst,$
        pos=imdisp_pos(inten,margin=margin)

end
