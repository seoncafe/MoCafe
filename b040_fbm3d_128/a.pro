a1 = mrdfits('M010b040_001.fits.gz',0,hdr1)
a2 = mrdfits('M020b040_001.fits.gz',0,hdr2)
a3 = mrdfits('M040b040_001.fits.gz',0,hdr3)
a1 = a1/mean(a1)
a2 = a2/mean(a2)
a3 = a3/mean(a3)
mwrfits,a1,'M010_001.fits.gz',hdr1
mwrfits,a2,'M020_001.fits.gz',hdr2
mwrfits,a3,'M040_001.fits.gz',hdr3
