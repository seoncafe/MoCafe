#!/bin/sh
#ds9 -width 420 -colorbar vertical -colorbar fontsize 16 -3d -3d az 30 -3d el 30 -3d border color white $1.fits.gz \
#    -zoom to fit -cmap bb -scale mode zmax -scale sqrt -print filename $1.ps -print
#ds9 -width 420 -colorbar vertical -colorbar fontsize 16 -plot title {test} -3d -3d az 30 -3d el 30 -3d border color white $1.fits.gz \
#    -zoom to fit -cmap bb -scale mode zmax -scale sqrt
ds9 -width 420 -colorbar vertical -colorbar fontsize 16 -3d -3d az 30 -3d el 30 -3d border color white $1.fits.gz \
    -zoom to fit -cmap bb -scale mode zmax -scale sqrt -regions load $2.reg
