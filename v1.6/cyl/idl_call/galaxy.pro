;--------
function colorbar_pos,pos,dx=dx,width=witdh
   if n_elements(dx)    ne 1 then dx    = (pos[2]-pos[0])*0.02d0
   if n_elements(width) ne 1 then width = (pos[2]-pos[0])*0.05d0
   pos   = [pos[2]+dx,pos[1],pos[2]+width,pos[3]]
   return,pos
end
function string_set,str,stringlength=stringlength
   if n_elements(stringlength) ne 1 then stringlength = 128
   return,strn(str,length=stringlength,padtype=0,padchar=' ')
end
;---------

function galaxy,plot=plot

no_photons    = 5d6
nprint        = 1000000L

hgg           = 0.41d0
albedo        = 0.4d0
luminosity    = 1d48
dust1         = 'exponential'
tauface1      = 0.8d0
dust_rscale1  = 6.6d0
dust_zscale1  = 0.25d0
dust_rmax1    = 18.0d0
;dust_zmax1    = 7.0d0
dust_zmax1    = 18.0d0
dust2         = 'exponential'
tauface2      = 0.0d0
;tauface2      = 0.2d0
dust_rscale2  = 6.6d0
dust_zscale2  = 1.50d0
dust_rmax2    = 18.0d0
;dust_zmax2    = 7.0d0
dust_zmax2    = 18.0d0
rmax          = 18.0d0
;zmax          = 7.0d0
zmax          = 18.0d0
nr            = 30L
np            = 1L
nz            = 61L
disk_name     = 'exponential'
disk_rscale   = 4.4d0
disk_zscale   = 0.5d0
disk_rmax     = 18.0d0
;disk_zmax     = 7.0d0
disk_zmax     = 18.0d0
bulge_name    = 'sersic'
sersic_index  = 2.5d0
Reff          = 2.0d0
axial_ratio   = 0.5d0
BulgeToDisk   = 0.33d0

inclination_angle = 89.0d0
position_angle    = 15.0d0
phase_angle       = 0.0d0
distance          = 9500.0d0
nxim          = 500L
nyim          = 250L
dxim          = 1.5d0/60d0/60d0
dyim          = dxim
; left and right sides of the output image are fold if tile_angle = position_angle = 0.0
; with half number of photons, we can obtain the same quality.
left_right_fold = 1

;--- 2016-03-20, now left-right side are doubled regardless of position_angle.
;if position_angle ne 0.0d0 then begin
;   left_right_fold = 0
;endif

; PSF funcion
; if you don't want to convolve with a PSF, then set psf_file = '', do not define psf_file, or undefine,psf_file
; Warning - 2015/11/21, in fims.kasi.re.kr machine, IDL version does not work with compressed fits file.
;           standalone version has no problem with a compressed fits PSF file.
;           I don't know why.
if getenv('HOST') eq 'fims.kasi.re.kr' then begin
   ;psf_file = '../PSF/PSF_GALEX_FUV.fits'
   psf_file = '../PSF/PSF_TNG_V.fits'
endif else begin
   ;psf_file = '../PSF/PSF_GALEX_FUV.fits.gz'
   psf_file = '../PSF/PSF_TNG_V.fits.gz'
endelse

;----------
; These output arrays should be decalred as single precision.
im_scatt      = fltarr(nxim,nyim)
im_direc      = fltarr(nxim,nyim)
im_scatt_sig  = fltarr(nxim,nyim)
im_direc_sig  = fltarr(nxim,nyim)

if n_elements(psf_file) eq 0 then psf_file = ''
; To match the string length to that defined in fortran and C programs.
dust1         = string_set(dust1)
dust2         = string_set(dust2)
disk_name     = string_set(disk_name)
bulge_name    = string_set(bulge_name)
psf_file       = string_set(psf_file)

; We do not have phi-anlge dependence in grid system.
pmax          = 0.0d0
galaxy_idl,double(no_photons),long(nprint),double(hgg),double(albedo),double(luminosity),$
           dust1,double(tauface1),double(dust_rscale1),double(dust_zscale1),double(dust_rmax1),double(dust_zmax1),$
           dust2,double(tauface2),double(dust_rscale2),double(dust_zscale2),double(dust_rmax2),double(dust_zmax2),$
           double(rmax),double(pmax),double(zmax),long(nr),long(np),long(nz),$
           disk_name,double(disk_rscale),double(disk_zscale),double(disk_rmax),double(disk_zmax),$
           bulge_name,double(sersic_index),double(Reff),double(axial_ratio),double(BulgeToDisk),$
           double(inclination_angle),double(position_angle),double(phase_angle),double(distance),$
           long(nxim),long(nyim),double(dxim),double(dyim),im_scatt,im_direc,im_scatt_sig,im_direc_sig,$
           psf_file,long(left_right_fold)

tot = im_scatt + im_direc
sig = sqrt(im_scatt_sig^2 + im_direc_sig^2)
out = {scatt:im_scatt,direc:im_direc,scatt_sig:im_scatt_sig,direc_sig:im_direc_sig,tot:tot,sig:sig}

if keyword_set(plot) then begin
   loadct,39
   if !d.name eq 'PS' then begin
      !p.font=0
      !p.charsize=0.8
      device,/color,bit=8,/helve,xsize=20,ysize=11
      ;neg = 1
      ;loadct,39
      ;loadct,5
      margin=0.07
   endif else begin
      !p.charsize=1.2
      margin=0.07
   endelse

   nx = (size(out.tot))[1]
   ny = (size(out.tot))[2]
   x  = dindgen(nx)*dxim-(nx-1)/2d0*dxim
   y  = dindgen(ny)*dyim-(ny-1)/2d0*dyim
   xr = minmax(x)
   yr = minmax(y)
   xtit = 'X (arcmin)'
   ytit = 'Y (arcmin)'

   !p.multi=[0,2,2]
   !x.title = xtit
   !y.title = ytit
   range = [0,max(out.tot)]
   imdisp,out.scatt,/er,/ax,tit='Scattered Light',negative=neg,margin=margin,out_pos=pos,xr=xr,yr=yr
   cgcolorbar,/vertical,/right,range=range,pos=colorbar_pos(pos),reverse=neg
   imdisp,out.direc,/er,/ax,tit='Direct Light',negative=neg,margin=margin,out_pos=pos,xr=xr,yr=yr
   cgcolorbar,/vertical,/right,range=range,pos=colorbar_pos(pos),reverse=neg
   imdisp,out.tot,/er,/ax,tit='Total Light',negative=neg,margin=margin,out_pos=pos,xr=xr,yr=yr
   cgcolorbar,/vertical,/right,range=range,pos=colorbar_pos(pos),reverse=neg
   snr = out.tot[*,*]
   w   = where(out.sig le 0.d0, nw, compl=wc)
   snr[wc] = out.tot[wc]/out.sig[wc]
   if nw ge 1 then snr[w] = 0.0d0
   range = [1.0, median(snr[wc]) + stddev(snr[wc],/nan)*5]
   imdisp,snr,/er,/ax,range=range,tit='Signal-to-Noise Ratio',negativ=neg,margin=margin,out_pos=pos,xr=xr,yr=yr
   cgcolorbar,/vertical,/right,range=range,pos=colorbar_pos(pos),reverse=neg
   if !d.name eq 'PS' then begin
      device,/close
      cleanplot,/sil
   endif
endif

return, out
end
