fwhm = 3.0
range = [0.0, 7300.0]
a = get_ngc891(/fold)
restore,'best_fit.sav'
a1 = filter_image(a.map,fwhm=fwhm)
b1 = filter_image(map_best,fwhm=fwhm)

!p.multi=[0,2,1]
imdisp,a1,/er,/ax,range=range
imdisp,b1,/er,/ax,range=range
