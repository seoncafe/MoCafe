function string_set,str,stringlength=stringlength
   if n_elements(stringlength) ne 1 then stringlength = 100
   return,strn(str,length=stringlength,padtype=0,padchar=' ')
end
;---------

function galaxy_car

no_photons    = 1d7
nprint        = 1000000
hgg           = 0.41d0
albedo        = 0.4d0
luminosity    = 1d38
iseed         = 0L
dust1         = 'exponential'
tauface1      = 2.0d0
dust_rscale1  = 8.0d0
dust_zscale1  = 0.10d0
dust_rmax1    = 18.0d0
dust_zmax1    = 7.0d0
dust2         = 'exponential'
tauface2      = 2.0d0
dust_rscale2  = 8.0d0
dust_zscale2  = 1.0d0
dust_rmax2    = 18.0d0
dust_zmax2    = 7.0d0
xmax          = 18.0d0
ymax          = 18.0d0
zmax          = 7.0d0
nx            = 80L
ny            = 80L
nz            = 500L
source_name   = 'exponential'
source_rscale = 5.0d0
source_zscale = 0.1d0
source_rmax   = 16.0d0
source_zmax   = 7.0d0
obs_angles    = [0.0d0, 0.2d0, 0.0d0]
obs_distance  = 9500.0d0
nxim          = 500L
nyim          = 350L
dxim          = 1.5d0/60d0/60d0
dyim          = dxim
;im_scatt      = dblarr(nxim,nyim)
;im_direc      = dblarr(nxim,nyim)
im_scatt      = fltarr(nxim,nyim)
im_direc      = fltarr(nxim,nyim)

dust1         = string_set(dust1)
dust2         = string_set(dust2)
source_name   = string_set(source_name)

;galaxy_car_idl,double(no_photons),long(nprint),double(hgg),double(albedo),double(luminosity),long64(iseed),$
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

return, {scatt:im_scatt,direc:im_direc}
end
