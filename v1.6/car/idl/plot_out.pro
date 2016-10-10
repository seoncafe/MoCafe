fname = 'out.fits.gz'
;fname = 'out1.fits.gz'
;fname = 'out2.fits.gz'
;fname = 'out3.fits.gz'
;fname = 'out4.fits.gz'
a=mrdfits(fname,0,hdr)
b=mrdfits(fname,1,hdr)
a_sig=mrdfits(fname,2,hdr)
b_sig=mrdfits(fname,3,hdr)
c = a+b
c_sig = sqrt(a_sig^2 + b_sig^2)
snr = div(c,c_sig)
c = c/1d-10
;------
h = hist1(c,min=0,max=max(c),bin=max(c)/100d0,xarr=xarr)
maxx  = interpol(reverse(xarr),reverse(h),max(h)*0.01d0)
h = hist1(c,min=0,max=maxx,bin=maxx/100d0,xarr=xarr)
range = minmax(xarr)

loadct,39
!p.multi=[0,2,1]
if !d.name eq 'PS' then begin
   !p.font=0
   device,/color,bit=8,/helve
   negative = 1
   set_thick,2.0
   !p.charsize=1.0
endif else begin
   negative = 0
   !p.charsize = 1.7
endelse

dx = sqrt(fxpar(hdr,'CD1_1')^2+fxpar(hdr,'CD1_2')^2) * 60d0
dy = sqrt(fxpar(hdr,'CD2_2')^2+fxpar(hdr,'CD2_2')^2) * 60d0
nx = n_elements(c[*,0])
ny = n_elements(c[0,*])
x  = dindgen(nx)*dx-(nx-1)/2d0*dx
y  = dindgen(ny)*dy-(ny-1)/2d0*dy
xr = minmax(x)
yr = minmax(y)
xtit = 'X (arcmin)'
ytit = 'Y (arcmin)'

margin = 0.07
imdisp,c,/er,/ax,xr=xr,yr=yr,range=range,margin=margin,out_pos=pos,negative=negative,tit='Intensity',xtit=xtit,ytit=ytit
width = (pos[2]-pos[0])*0.1d0
dx    = (pos[2]-pos[0])*0.02d0
pos=[pos[2]+dx,pos[1],pos[2]+width,pos[3]]
cgcolorbar,/vertical,/right,range=range,pos=pos,reverse=negative

;------
h_snr = hist1(snr,min=0,max=max(snr),bin=max(snr)/100d0,xarr=xarr)
maxx  = interpol(reverse(xarr),reverse(h_snr),max(h_snr)*0.01d0)
h_snr = hist1(snr,min=0,max=maxx,bin=maxx/100d0,xarr=xarr)

undefine,range
range = minmax(xarr)
;range = [0,12]
imdisp,snr,/er,/ax,xr=xr,yr=yr,range=range,margin=margin,out_pos=pos,negative=negative,tit='SNR',xtit=xtit,ytit=ytit
width = (pos[2]-pos[0])*0.1d0
dx    = (pos[2]-pos[0])*0.02d0
pos=[pos[2]+dx,pos[1],pos[2]+width,pos[3]]
cgcolorbar,/vertical,/right,range=range,pos=pos,reverse=negative

end
