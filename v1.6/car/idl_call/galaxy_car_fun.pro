function string_set,str,stringlength=stringlength
   if n_elements(stringlength) ne 1 then stringlength = 100
   return,strn(str,length=stringlength,padtype=0,padchar=' ')
end
;---------

function galaxy_car_fun, x, params,$
         noPSF=noPSF,$
         luminosity=luminosity,iseed=iseed,$
         dust1=dust1,dust2=dust2,source_name=source_name,$
         nx=nx,ny=ny,nz=nz,$
         obs_distance=obs_distance,$
         nxim=nxim,nyim=nyim,dxim=dxim,dyim=dyim,$
         no_photons=no_photons,nprint=nprint

if n_elements(no_photons) ne 1 then no_photons = 1d7
if n_elements(luminosity) ne 1 then luminosity = 1d38
if n_elements(iseed)      ne 1 then iseed = 0L
if n_elements(dust1)      ne 1 then dust1 = 'exponential'
if n_elements(dust2)      ne 1 then dust2 = 'exponential'
if n_elements(nx)         ne 1 then nx = 80L
if n_elements(ny)         ne 1 then ny = 80L
if n_elements(nz)         ne 1 then nz = 500L
if n_elements(source_name)  ne 1 then source_name  = 'exponential'
if n_elements(obs_distance) ne 1 then obs_distance = 9500.0d0
if n_elements(nxim)   ne 1 then nxim = 500L
if n_elements(nyim)   ne 1 then nyim = 350L
if n_elements(dxim)   ne 1 then dxim = 1.5d0/60d0/60d0
if n_elements(dyim)   ne 1 then dyim = dxim
if n_elements(nprint) ne 1 then nprint = 0L

if n_params() ne 2 then begin
   tauface1      = 2.0d0
   tauface2      = 2.0d0
   dust_rscale1  = 8.0d0
   dust_zscale1  = 0.10d0
   dust_rscale2  = 8.0d0
   dust_zscale2  = 1.0d0
   dust_rmax     = 18.0d0
   dust_zmax     = 7.0d0
   source_rscale = 5.0d0
   source_zscale = 0.1d0
   source_rmax   = 16.0d0
   inclination   = 89.8d0
   hgg           = 0.41d0
   albedo        = 0.4d0
endif else begin
   tauface1      = params[0]
   tauface2      = params[1]
   dust_rscale1  = params[2]
   dust_zscale1  = params[3]
   dust_rscale2  = params[4]
   dust_zscale2  = params[5]
   dust_rmax     = params[6]
   dust_zmax     = params[7]
   source_rscale = params[8]
   source_zscale = params[9]
   source_rmax   = params[10]
   inclination   = params[11]
   hgg           = params[12]
   albedo        = params[13]
endelse

dust_rmax1  = dust_rmax
dust_zmax1  = dust_zmax
dust_rmax2  = dust_rmax
dust_zmax2  = dust_zmax
xmax        = dust_rmax
ymax        = dust_rmax
zmax        = dust_zmax
source_zmax = dust_zmax
obs_angles  = [0.0d0, 90d0-inclination, 0.0d0]

;im_scatt      = dblarr(nxim,nyim)
;im_direc      = dblarr(nxim,nyim)
im_scatt      = fltarr(nxim,nyim)
im_direc      = fltarr(nxim,nyim)

dust1         = string_set(dust1)
dust2         = string_set(dust2)
source_name   = string_set(source_name)

galaxy_car_idl,double(no_photons),long(nprint),double(hgg),double(albedo),double(luminosity),long(iseed),$
               dust1,double(tauface1),double(dust_rscale1),double(dust_zscale1),$
                                      double(dust_rmax1),double(dust_zmax1),$
               dust2,double(tauface2),double(dust_rscale2),double(dust_zscale2),$
                                      double(dust_rmax2),double(dust_zmax2),$
               double(xmax),double(ymax),double(zmax),long(nx),long(ny),long(nz),$
               source_name,double(source_rscale),double(source_zscale),$
                           double(source_rmax),double(source_zmax),$
               double(obs_angles),double(obs_distance),$
               long(nxim),long(nyim),double(dxim),double(dyim),im_scatt,im_direc

;galaxy_car_idl,no_photons,long(nprint),hgg,albedo,luminosity,long(iseed),$
;               dust1,tauface1,dust_rscale1,dust_zscale1,$
;                                      dust_rmax1,dust_zmax1,$
;               dust2,tauface2,dust_rscale2,dust_zscale2,$
;                                      dust_rmax2,dust_zmax2,$
;               xmax,ymax,zmax,long(nx),long(ny),long(nz),$
;               source_name,source_rscale,source_zscale,$
;                           source_rmax,source_zmax,$
;               obs_angles,obs_distance,$
;               long(nxim),long(nyim),dxim,dyim,im_scatt,im_direc

out_image = im_scatt + im_direc
if not(keyword_set(noPSF)) then begin
   kernel    = mrdfits('/Users/kiseon/Galaxies/GALEX_PSF/PSFfuv.fits.gz')
   out_image = convolve(out_image,kernel)
endif

intensity_scale = 5863.8510/max(total(out_image,1)/nyim)
out_image = out_image * intensity_scale

return, out_image

end
