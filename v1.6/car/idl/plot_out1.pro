a1 = mrdfits('out.fits.gz',0,hdr)/1d-10
a2 = mrdfits('out.fits.gz',1,hdr)/1d-10
a  = a1+a2

dx = sqrt(fxpar(hdr,'CD1_1')^2+fxpar(hdr,'CD1_2')^2) * 60d0
dy = sqrt(fxpar(hdr,'CD2_2')^2+fxpar(hdr,'CD2_2')^2) * 60d0
nx = n_elements(a1[*,0])
ny = n_elements(a1[0,*])
x  = dindgen(nx)*dx-(nx-1)/2d0*dx
y  = dindgen(ny)*dy-(ny-1)/2d0*dy

!p.multi=[0,3,2]
loadct,39
if !d.name eq 'PS' then begin
   !p.font=0
   device,/color,bit=8,/helve
   negative = 1
   set_thick,2.0
   !p.charsize=1.0
endif else begin
   negative = 0
   !p.charsize = 2.2
   margin = 0.05
endelse

xr = minmax(x)
yr = minmax(y)

xtit = 'X (arcmin)'
ytit = 'Y (arcmin)'
imdisp,a1,/er,/ax,range=minmax(a1),xr=xr,yr=yr,xtit=xtit,ytit=ytit,margin=margin,negative=negative,tit='scattered'
imdisp,a2,/er,/ax,range=minmax(a2),xr=xr,yr=yr,xtit=xtit,ytit=ytit,margin=margin,negative=negative,tit='direct'
imdisp,a, /er,/ax,range=minmax(a),xr=xr,yr=yr,xtit=xtit,ytit=ytit,margin=margin,negative=negative,tit='total'

xtit = 'Y (arcmin)'
ytit = 'Profile'
plot,y,total(a1,1),/xst,psym=10,xtit=xtit,ytit=ytit,pos=imdisp_pos(a,margin=margin),tit='scattered'
plot,y,total(a2,1),/xst,psym=10,xtit=xtit,ytit=ytit,pos=imdisp_pos(a,margin=margin),tit='direct'
plot,y,total(a, 1),/xst,psym=10,xtit=xtit,ytit=ytit,pos=imdisp_pos(a,margin=margin),tit='total'

end
