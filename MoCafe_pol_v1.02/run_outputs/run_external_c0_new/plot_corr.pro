pro plot_corr,nsample

if n_elements(nsample) eq 0 then nsample = 1
flist = file_search('M010_'+strn(nsample,length=3,padchar='0')+'K_t*_peel.fits.gz')
n     = n_elements(flist)

tau  = strmid(flist,11,4)*1d-3
flux = dblarr(n)
for i=0,n-1 do begin
   fname   = flist[i]
   a       = mrdfits(fname,0,/sil)
   flux[i] = total(a)
endfor

plot,tau,flux,psym=-4

end
