fname = 'out.fits.gz'
;fname = 'out_nobulge.fits.gz'

a=mrdfits(fname,0,hdr)
b=mrdfits(fname,1,hdr)
a_sig=mrdfits(fname,2,hdr)
b_sig=mrdfits(fname,3,hdr)
c = a+b
c_sig = sqrt(a_sig^2 + b_sig^2)
snr = div(c,c_sig)
c   = c/1d-10
;------
h    = histogram(c,min=0,max=max(c,/nan),bin=max(c,/nan)/100d0,locations=xarr)
minx = 0.0d0
maxx = interpol(reverse(xarr),reverse(h),max(h)*0.0001d0)
h    = histogram(c,min=minx,max=maxx,bin=maxx/100d0,locations=xarr)
range = [min(c[where(c gt 0.0)]),maxx]

loadct,39
!p.multi=[0,2,1]
if !d.name eq 'PS' then begin
   ;loadct,5
   !p.font=0
   device,/color,bit=8,/helve
   ;negative = 1
   set_thick,2.0
   !p.charsize=0.8
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
scale = 10d0^long(alog10(max(range)))
cgcolorbar,/vertical,/right,range=range/scale,pos=colorbar_pos(pos),reverse=negative

;------
h_snr = histogram(snr,min=0,max=max(snr),bin=max(snr)/100d0,locations=xarr)
maxx  = interpol(reverse(xarr),reverse(h_snr),max(h_snr)*0.01d0)
h_snr = histogram(snr,min=0,max=maxx,bin=maxx/100d0,locations=xarr)

undefine,range
range = minmax(xarr)
;range = [0,12]
imdisp,snr,/er,/ax,xr=xr,yr=yr,range=range,margin=margin,out_pos=pos,negative=negative,tit='SNR',xtit=xtit,ytit=ytit
cgcolorbar,/vertical,/right,range=range,pos=colorbar_pos(pos),reverse=negative

end
