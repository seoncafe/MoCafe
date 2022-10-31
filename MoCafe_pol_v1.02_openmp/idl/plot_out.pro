a = mrdfits('../out.fits.gz',0,hdr)
b = mrdfits('../out.fits.gz',1)
c = a+b

nx = n_elements(a[*,0])
ny = n_elements(a[0,*])

xyad,hdr,0.0,0.0,xmin,ymin
xyad,hdr,nx-1,ny-1,xmax,ymax
xr = [xmin,xmax]
yr = [ymin,ymax]

xx = findgen(nx,ny)
yy = findgen(nx,ny)
for i=0,nx-1 do begin
for j=0,ny-1 do begin
   xyad,hdr,i,j,xx1,yy1
   xx[i,j]=xx1
   yy[i,j]=yy1
endfor
endfor
rr = sqrt(xx*xx + yy*yy)

if !d.name eq 'PS' then begin
  set_thick,2
  !p.font = 0
  device,ysize=16,/color,bit=8,/helvetica,file='out.ps'
endif

!p.charsize=1.8
!p.multi=[0,2,2]
loadct,5
margin=0.05

imdisp,asinh(a/median(a)),/er,/ax,margin=margin,$
       xtit='X (degree)',ytit='Y (degree)',tit='Scattered Light',xr=xr,yr=yr
imdisp,asinh(c/median(c)),/er,/ax,margin=margin,$
       xtit='X (degree)',ytit='Y (degree)',tit='Total Light',xr=xr,yr=yr
plot,rr,a,psym=3,/ylog,xr=[0,xmax],pos=imdisp_pos(a,margin=margin),$
      xtit='Radius (degree)',ytit='Scattered Intensity'
plot,rr,c,psym=3,/ylog,xr=[0,xmax],pos=imdisp_pos(a,margin=margin),$
      xtit='Radius (degree)',ytit='Total Intensity'

if !d.name eq 'PS' then device,/close

end
