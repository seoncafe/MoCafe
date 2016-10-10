@/Users/kiseon/Galaxies/NGC_891/get_ngc891.pro
pro fit1

fold = 1
a   = get_ngc891(fold=fold)
map = a.map[*,75:424]

;  0:  tauface1      = 2.0d0
;  1:  tauface2      = 2.0d0
;  2:  dust_rscale1  = 8.0d0
;  3:  dust_zscale1  = 0.15d0
;  4:  dust_rscale2  = 8.0d0
;  5:  dust_zscale2  = 2.0d0
;  6:  dust_rmax     = 20.0d0
;  7:  dust_zmax     = 7.0d0
;  8:  source_rscale = 6.0d0
;  9:  source_zscale = 0.1d0
; 10:  source_rmax   = 18.0d0
; 11:  inclination   = 89.8d0
; 12:  hgg           = 0.41d0
; 13:  albedo        = 0.4d0
p0 = [2.0d0, 2.0d0, 8.0d0, 0.15d0, 8.0d0, 2.0d0, 20.0d0, 7.0d0, 6.0d0, 0.1d0, 18.0d0, 89.8d0, 0.41d0, 0.4d0]
perror = 1.0
parinfo = replicate({value:0.d0,fixed:0,limited:[0,0],limits:[0.d,0],relstep:0.d0},14)
parinfo[*].value = p0
parinfo[*].fixed = 1
parinfo[3].fixed = 0
parinfo[3].limited = [1,1]
parinfo[3].limits  = [0.05d0, 0.25d0]
parinfo[3].relstep = 0.1d0

x = map
params = mpfitfun('galaxy_car_fun',x,map,err,p0,yfit=yfit,perror=perror,parinfo=parinfo)
print,'best-fit parameters'
print,params

!p.multi=[0,2,1]
loadct,39
imdisp,map,/er,/ax
imdisp,yfit,/er,/ax

end
