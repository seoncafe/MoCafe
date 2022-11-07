pro gen_fbm3d,mach=mach,noexcute=noexcute,bvalue=bvalue,$
              iseed=iseed,ngrid=ngrid,seq=seq
;
; version 2013-09-20
;   default is for natural mixing mode (bvalue = 0.4)
;
if n_elements(iseed)  eq 0 then iseed  = 0
if n_elements(bvalue) eq 0 then bvalue = 0.4
if n_elements(mach)   eq 0 then mach   = 1.0
if n_elements(ngrid)  eq 0 then ngrid  = 256
if n_elements(seq)    eq 0 then seq    = 0

nx   = ngrid
st_m = strn(long(mach*10),   length=3,padtype=1,padchar='0')
st_b = strn(long(bvalue*100),length=3,padtype=1,padchar='0')
st_s = strn(long(seq),       length=3,padtype=1,padchar='0')

comm    = 'M'+st_m+'b'+st_b+'_'+st_s
parfile = comm+'.in'
outfile = comm+'.fits.gz'

openw,1,parfile
printf,1,'&input'
printf,1,'  iseed    = '+strtrim(string(iseed,format='(i)'),2)+','
printf,1,'  bvalue   = '+string(bvalue,format='(f8.4)')+','
printf,1,'  mach     = '+string(mach,format='(f8.4)')+','
printf,1,'  nx       = '+string(nx,format='(i)')+','
printf,1,"  outfile  = '"+outfile+"',"
printf,1,"  out_mode = 1"
printf,1,'/'
close,1

print,'.....'
if not(keyword_set(noexcute)) then spawn,'../../fbm3d/make_fbm3d.x < '+parfile

end
