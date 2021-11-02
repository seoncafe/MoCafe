pro plot_peel, fname, all=all, old=old,ylog=ylog,ps=ps,median=median

if n_elements(fname) eq 0 then fname = 'model_001'
if keyword_set(all) then begin
   a0 = mrdfits(fname+'_peel.fits.gz',0,/sil)
   a1 = mrdfits(fname+'_peel.fits.gz',1,/sil)
   a  = a0 + a1
endif else begin
   a  = mrdfits(fname+'_peel.fits.gz',0,/sil)
endelse

if keyword_set(old) then begin
  a0  = mrdfits(fname+'_peel.fits.gz',0,/sil)
  a1  = mrdfits(fname+'_peel.fits.gz',1,/sil)
  a   = a0 + a1
  d   = mrdfits(fname+'_dens.fits.gz',0,hdr,/sil)
  dz  = fxpar(hdr,'delz')
  tau = total(d,3)*dz
endif else begin
  if keyword_set(all) then begin
     a0 = mrdfits(fname+'_peel.fits.gz',0,/sil)
     a1 = mrdfits(fname+'_peel.fits.gz',1,/sil)
     a  = a0 + a1
  endif else begin
     a  = mrdfits(fname+'_peel.fits.gz',0,/sil)
  endelse
  tau = mrdfits(fname+'_tau.fits.gz',0,hdr,/sil)
endelse

!p.multi=[0,2,2]
loadct,39
margin=0.06
img=lindgen(2,2)

mach = strmid(fname,1,3)*0.1d0
str  = strmid(fname,11,4)
if strlen(str) eq 4 then begin
   tauh = str * 0.001d0
   tau_str = string(tauh,format='(f5.3)')
endif else begin
   tauh = str * 0.01d0
   tau_str = string(tauh,format='(f4.2)')
endelse

if keyword_set(ps) then begin
   mydevice = !d.name
   set_thick,2.0
   !p.font = 0
   !p.charsize = 1.2
   psfile = fname + '.ps'
   set_plot,'ps'
   device,ysize=16,yoffset=2,/color,bit=16,/helvetica,file=psfile
endif

;imdisp,asinh(a/median(a)),/er,/ax
;imdisp,asinh(tau/median(tau)),/er,/ax
imdisp,a,  /er,/ax, tit='Intensity', margin=margin
imdisp,tau,/er,/ax, tit=textoidl('\tau_{dust}'), margin=margin

ex = long(alog10(max(tau)))
if (ex lt 1d0) then begin
   ex    = ex-1d0
   tau10 = tau/10d0^ex
   xtit = textoidl('\tau_{dust} \times10^{'+string(ex,format='(I2)')+'}')
endif else begin
   tau10 = tau
   xtit = textoidl('\tau_{dust}')
endelse

ex = long(alog10(max(a)))
if (ex lt 1d0) then begin
   ex   = ex-1d0
   a10  = a/10d0^ex
   ytit = textoidl('Intensity \times10^{'+string(ex,format='(I2)')+'}')
endif else begin
   a10  = a
   ytit = textoidl('Intensity')
endelse

xr=[0,max(tau10)]
yr=[0,max(a10)]

tit = textoidl('M = '+string(mach,format='(f3.1)')+', '+'\tau_h = '+tau_str)
if keyword_set(ylog) then yr = [1d-4,1]*max(a10)
plot,tau10,a10,psym=3, xr=xr,yr=yr,/xst, ylog=ylog,xtit=xtit,ytit=ytit,tit=tit,$
            pos=imdisp_pos(img,margin=margin)

xmin = 0.0d0
xmax = max(tau10)
nx = 50
dx = (xmax-xmin)/(nx-1d0)
x1 = dindgen(nx)*dx + xmin
x2 = x1 + dx
x  = (x1 + x2)/2.0d0

y   = dblarr(nx) + !values.f_nan
err = dblarr(nx) + !values.f_nan
for i=0,nx-1 do begin
   w = where(tau10 ge x1[i] and tau10 lt x2[i],nw)
   if nw gt 0 then begin
      if keyword_set(median) then y[i]   = median(a10[w]) $
      else y[i]   = mean(a10[w])
      err[i] = stddev(a10[w])
   endif
endfor

plot,x,y,xr=xr,yr=yr,/xst,ylog=ylog,xtit=xtit, ytit=ytit,$
         pos=imdisp_pos(img,margin=margin)
errplot,x,y-err,y+err

if keyword_set(ps) then begin
  device,/close
  set_plot,mydevice
endif

end
