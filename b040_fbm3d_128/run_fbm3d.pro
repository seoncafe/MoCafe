ngrid   = 128
;mach_in = [0.5,1.0,2.0,3.0,4.0]
;mach_in = [6.0,8.0,10.0]
;mach_in = [12.0]
;mach_in = [16.0,20.0]
mach_in = [16.0]
nmach   = n_elements(mach_in)
nseq    = 10
seq_in  = lindgen(nseq)+1
bvalue  = 0.4

for i=0,nmach-1 do for k=0,nseq-1 do $
    gen_fbm3d,mach=mach_in[i],seq=seq_in[k],bvalue=bvalue,ngrid=ngrid
