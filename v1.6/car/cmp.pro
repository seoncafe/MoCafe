; to compare with the result calculated with cylindrical code.
;
a1 = mrdfits('out.fits.gz',0)/1d-10
a2 = mrdfits('out.fits.gz',1)/1d-10
a  = a1+a2
b1 = mrdfits('../extraplanar_dust/Galaxy_cylsym/out.fits.gz',0)/1d-10
b2 = mrdfits('../extraplanar_dust/Galaxy_cylsym/out.fits.gz',1)/1d-10
b  = b1+b2

!p.multi=[0,3,2]
plot,a1,b1,psym=3,/iso
oplot,!x.crange,!y.crange,color=cgcolor('red')
plot,a,b,psym=3,/iso,xr=[0,20],yr=[0,20]
oplot,!x.crange,!y.crange,color=cgcolor('red')

loadct,39
imdisp,a1,/er,/ax,range=minmax([a1,b1])
imdisp,b1,/er,/ax,range=minmax([a1,b1])
imdisp,a, /er,/ax,range=minmax([a2,b2])
imdisp,b, /er,/ax,range=minmax([a2,b2])

end
