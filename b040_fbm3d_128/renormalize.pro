a1 = mrdfits('M010b040_001.fits.gz',0,hdr1)
a2 = mrdfits('M020b040_001.fits.gz',0,hdr2)
a3 = mrdfits('M040b040_001.fits.gz',0,hdr3)
a4 = mrdfits('M060b040_001.fits.gz',0,hdr4)
a5 = mrdfits('M080b040_001.fits.gz',0,hdr5)
a6 = mrdfits('M100b040_001.fits.gz',0,hdr6)
a7 = mrdfits('M120b040_001.fits.gz',0,hdr7)
a8 = mrdfits('M160b040_001.fits.gz',0,hdr8)
a9 = mrdfits('M200b040_001.fits.gz',0,hdr9)
a1 = a1/mean(a1)
a2 = a2/mean(a2)
a3 = a3/mean(a3)
a4 = a4/mean(a4)
a5 = a5/mean(a5)
a6 = a6/mean(a6)
a7 = a7/mean(a7)
a8 = a8/mean(a8)
a9 = a9/mean(a9)
mwrfits,a1,'M010_001.fits',hdr1,/create
mwrfits,a2,'M020_001.fits',hdr2,/create
mwrfits,a3,'M040_001.fits',hdr3,/create
mwrfits,a4,'M060_001.fits',hdr4,/create
mwrfits,a5,'M080_001.fits',hdr5,/create
mwrfits,a6,'M100_001.fits',hdr6,/create
mwrfits,a7,'M120_001.fits',hdr7,/create
mwrfits,a8,'M160_001.fits',hdr8,/create
mwrfits,a9,'M200_001.fits',hdr9,/create
