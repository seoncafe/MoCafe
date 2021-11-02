a1 = mrdfits('rec111_rot000_tau50_hgg00_a10_obs.fits.gz',0)
a2 = mrdfits('rec111_rot000_tau50_hgg00_a10_obs.fits.gz',1)
b1 = mrdfits('rec111_rot045_tau50_hgg00_a10_obs.fits.gz',0)
b2 = mrdfits('rec111_rot045_tau50_hgg00_a10_obs.fits.gz',1)
c1 = mrdfits('rec111_rot090_tau50_hgg00_a10_obs.fits.gz',0)
c2 = mrdfits('rec111_rot090_tau50_hgg00_a10_obs.fits.gz',1)
d1 = mrdfits('rec111_rot135_tau50_hgg00_a10_obs.fits.gz',0)
d2 = mrdfits('rec111_rot135_tau50_hgg00_a10_obs.fits.gz',1)
e1 = mrdfits('rec111_rot180_tau50_hgg00_a10_obs.fits.gz',0)
e2 = mrdfits('rec111_rot180_tau50_hgg00_a10_obs.fits.gz',1)

a = a1 + a2
b = b1 + b2
c = c1 + c2
d = d1 + d2
e = e1 + e2

a1 = asinh(a1/mean(a1))
b1 = asinh(b1/mean(b1))
c1 = asinh(c1/mean(c1))
d1 = asinh(d1/mean(d1))
e1 = asinh(e1/mean(e1))

a  = asinh(a/mean(a))
b  = asinh(b/mean(b))
c  = asinh(c/mean(c))
d  = asinh(d/mean(d))
e  = asinh(e/mean(e))

loadct,39
!p.multi=[0,4,2]
!p.charsize=4

xr = [-1,1]
yr = [-1,1]
margin = 0.04

imdisp,a1,/er,/ax,xr=xr,yr=yr,margin=margin,tit=textoidl('I_{scatt} (\gamma = 0^{\circ})')
imdisp,b1,/er,/ax,xr=xr,yr=yr,margin=margin,tit=textoidl('I_{scatt} (\gamma = -45^{\circ})')
imdisp,c1,/er,/ax,xr=xr,yr=yr,margin=margin,tit=textoidl('I_{scatt} (\gamma = -90^{\circ})')
imdisp,d1,/er,/ax,xr=xr,yr=yr,margin=margin,tit=textoidl('I_{scatt} (\gamma = -135^{\circ})')
;imdisp,e1,/er,/ax,xr=xr,yr=yr,margin=margin,tit=textoidl('I_{scatt} (\gamma = -180^{\circ})')

imdisp,a,/er,/ax,xr=xr,yr=yr,margin=margin,tit=textoidl('I_{tot} (\gamma = 0^{\circ})')
imdisp,b,/er,/ax,xr=xr,yr=yr,margin=margin,tit=textoidl('I_{tot} (\gamma = -45^{\circ})')
imdisp,c,/er,/ax,xr=xr,yr=yr,margin=margin,tit=textoidl('I_{tot} (\gamma = -90^{\circ})')
imdisp,d,/er,/ax,xr=xr,yr=yr,margin=margin,tit=textoidl('I_{tot} (\gamma = -135^{\circ})')
;imdisp,e,/er,/ax,xr=xr,yr=yr,margin=margin,tit=textoidl('I_{tot} (\gamma = -180^{\circ})')

end
