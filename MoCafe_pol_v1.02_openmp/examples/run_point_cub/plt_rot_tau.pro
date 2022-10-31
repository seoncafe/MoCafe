a = mrdfits('rec111_rot000_tau50_hgg00_a10_tau.fits.gz',0)
b = mrdfits('rec111_rot045_tau50_hgg00_a10_tau.fits.gz',0)
c = mrdfits('rec111_rot090_tau50_hgg00_a10_tau.fits.gz',0)
d = mrdfits('rec111_rot135_tau50_hgg00_a10_tau.fits.gz',0)
e = mrdfits('rec111_rot180_tau50_hgg00_a10_tau.fits.gz',0)

;a  = asinh(a/mean(a))
;b  = asinh(b/mean(b))
;c  = asinh(c/mean(c))
;d  = asinh(d/mean(d))
;e  = asinh(e/mean(e))

loadct,39
!p.multi=[0,3,2]
!p.charsize=4

xr = [-1,1]
yr = [-1,1]
margin = 0.04

imdisp,a,/er,/ax,xr=xr,yr=yr,margin=margin,tit=textoidl('I_{tot} (\gamma = 0^{\circ})')
imdisp,b,/er,/ax,xr=xr,yr=yr,margin=margin,tit=textoidl('I_{tot} (\gamma = -45^{\circ})')
imdisp,c,/er,/ax,xr=xr,yr=yr,margin=margin,tit=textoidl('I_{tot} (\gamma = -90^{\circ})')
imdisp,d,/er,/ax,xr=xr,yr=yr,margin=margin,tit=textoidl('I_{tot} (\gamma = -135^{\circ})')
imdisp,e,/er,/ax,xr=xr,yr=yr,margin=margin,tit=textoidl('I_{tot} (\gamma = -180^{\circ})')

end
