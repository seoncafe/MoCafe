bvalue = 0.4d0
mach = [0.5d0,1d0,2d0,3d0,4d0,6d0,8d0,10d0,12d0,16d0,20d0]
seq  = lindgen(10)+1

nmach = n_elements(mach)
nseq  = n_elements(seq)

clumping_factor = dblarr(nmach,nseq)

for i=0,nmach-1 do begin
for j=0,nseq-1 do begin
   st_m = strn(long(mach[i]*10),length=3,padtype=1,padchar='0')
   st_b = strn(long(bvalue*100),length=3,padtype=1,padchar='0')
   st_s = strn(long(seq[j]),    length=3,padtype=1,padchar='0')
   comm  = 'M'+st_m+'b'+st_b+'_'+st_s
   fname = comm+'.fits.gz'
   a = mrdfits(fname,0,/silen)
   clumping_factor[i,j] = mean(a^2)/mean(a)^2
   ;print,mach[i],clumping_factor[i,j]
endfor
   print,mach[i],mean(clumping_factor[i,*]),stddev(clumping_factor[i,*])
endfor

end
