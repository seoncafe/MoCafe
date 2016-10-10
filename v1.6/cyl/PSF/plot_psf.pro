pro plot_psf,nuv=nuv

;---------
if keyword_set(nuv) then begin
  a  = mrdfits('PSFnuv_faint.fits.gz')
endif else begin
  a  = mrdfits('PSFfuv.fits.gz')
endelse
a = a/max(a)

delt = 4.16666700000000E-4
delt = delt * 60d0 * 60d0
n    = n_elements(a)
w    = array_indices(a,where(a eq max(a)))
idx0 = w[0]
idy0 = w[1]
idx  = reform((array_indices(a,lindgen(n)))[0,*])
idy  = reform((array_indices(a,lindgen(n)))[1,*])
help,idx,idy

ang  = sqrt((idx-idx0)^2 + (idy-idy0)^2)
ang_uniq = myuniq(ang)
nang = n_elements(ang_uniq)
psf  = dblarr(nang)

for i=0,nang-1 do begin
  w = where(ang eq ang_uniq[i],nw)
  if nw gt 0 then psf[i] = mean(a[w])
endfor
;---------
if keyword_set(nuv) then begin
  b  = mrdfits('PSF_GALEX_NUV.fits.gz')
endif else begin
  b  = mrdfits('PSF_GALEX_FUV.fits.gz')
endelse
b = b/max(b)

delt = 4.16666700000000E-4
delt = delt * 60d0 * 60d0
n    = n_elements(b)
w    = array_indices(b,where(b eq max(b)))
idx0 = w[0]
idy0 = w[1]
idx  = reform((array_indices(b,lindgen(n)))[0,*])
idy  = reform((array_indices(b,lindgen(n)))[1,*])
help,idx,idy

ang1  = sqrt((idx-idx0)^2 + (idy-idy0)^2)
ang1_uniq = myuniq(ang1)
nang1 = n_elements(ang1_uniq)
psf1  = dblarr(nang1)

for i=0,nang-1 do begin
  w = where(ang1 eq ang1_uniq[i],nw)
  if nw gt 0 then psf1[i] = mean(b[w])
endfor

;-------
!p.multi=[0,2,1]
plot,ang_uniq*delt,psf,yr=[1e-6,2.0],xr=[0,60],/ylog,/yst,$
     xtit='Angle (arcsec)',ytit='PSF'
oplot,ang1_uniq*delt,psf1,color=cgcolor('red')

dang    = ang_uniq[1:nang-1]-ang_uniq[0:nang-2]
psf_cum = total([0,psf[0:nang-2]*dang],/cum)
psf_cum = psf_cum/max(psf_cum)
radius  = findgen(61)
psf_r   = interpol(psf_cum,ang_uniq*delt,radius)
plot,ang_uniq*delt,psf_cum,yr=[0,1.1],xr=[0,60],/yst,$
;plot,radius,psf_r,yr=[0,1.1],xr=[0,60],/yst,$
     xtit='Angle (arcsec)',ytit='PSF'

dang    = ang1_uniq[1:nang1-1]-ang1_uniq[0:nang1-2]
psf_cum = total([0,psf1[0:nang1-2]*dang],/cum)
psf_cum = psf_cum/max(psf_cum)
radius  = findgen(61)
psf_r   = interpol(psf_cum,ang1_uniq*delt,radius)
oplot,ang1_uniq*delt,psf_cum,color=cgcolor('red')

end
