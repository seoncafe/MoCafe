!-------------------------------------------------------------------
  subroutine write_output(fname,out_mode,iseed,bvalue,mach,slope_ln,slope_gauss,sigma_lnrho, &
                          nx,ny,nz,array,karray1,karray2)
  implicit none
  character(len=*) :: fname
  integer :: iseed,nx,ny,nz,out_mode
  real    :: array(nx,ny,nz),karray1(nx/2+1,ny,nz),karray2(nx/2+1,ny,nz)
  real    :: bvalue,mach,slope_ln,slope_gauss,sigma_lnrho

  integer :: status,unit,blocksize,bitpix,naxis,naxes(3),group
  integer :: fpixel,nelements
  logical :: simple,extend

  call unlink(trim(fname))

  status = 0
  unit   = 1
  blocksize = 1
  call ftinit(unit,trim(fname),blocksize,status)

!---
  simple = .true.
  bitpix = -32
  naxis  = 3
  naxes(1) = nx
  naxes(2) = ny
  naxes(3) = nz
  extend   = .true.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  group  = 1
  fpixel = 1
  nelements = naxes(1)*naxes(2)*naxes(3)
  call ftppre(unit,group,fpixel,nelements,array,status)
  call ftpkyj(unit,'ISEED',iseed,  'Random Number iseed',status)
  call ftpkye(unit,'B-VALUE', bvalue, -8,'b-value',status)
  call ftpkye(unit,'MACH', mach, -8,'mach',status)
  call ftpkye(unit,'SLOPE_LN', slope_ln, -8,'Power Spectral Index of Lognormal Field',status)
  call ftpkye(unit,'SLOPE_GA', slope_gauss, -8,'Power Spectral Index of Gaussian Field',status)
  call ftpkye(unit,'SIGMA_LN', sigma_lnrho,-8,'Standard Deviation in Ln scale',status)
!---
  if (out_mode == 2 .or. out_mode == 3) then
     naxes(1)  = nx/2+1
     nelements = naxes(1)*naxes(2)*naxes(3)
     call ftiimg(unit,bitpix,naxis,naxes,status)
     call ftppre(unit,group,fpixel,nelements,karray1,status)
     call ftpkye(unit,'B-VALUE', bvalue, -8,'b-value',status)
     call ftpkye(unit,'MACH', mach, -8,'mach',status)
     call ftpkye(unit,'SLOPE_LN', slope_ln, -8,'Power Spectral Index of Lognormal Field',status)
     call ftpkye(unit,'SLOPE_GA', slope_gauss, -8,'Power Spectral Index of Gaussian Field',status)
     call ftpkye(unit,'SIGMA_LN', sigma_lnrho,-8,'Standard Deviation in Ln scale',status)
  endif
!---
  if (out_mode == 3) then
     naxes(1)  = nx/2+1
     nelements = naxes(1)*naxes(2)*naxes(3)
     call ftiimg(unit,bitpix,naxis,naxes,status)
     call ftppre(unit,group,fpixel,nelements,karray2,status)
     call ftpkye(unit,'B-VALUE', bvalue, -8,'b-value',status)
     call ftpkye(unit,'MACH', mach, -8,'mach',status)
     call ftpkye(unit,'SLOPE_LN', slope_ln, -8,'Power Spectral Index of Lognormal Field',status)
     call ftpkye(unit,'SLOPE_GA', slope_gauss, -8,'Power Spectral Index of Gaussian Field',status)
     call ftpkye(unit,'SIGMA_LN', sigma_lnrho,-8,'Standard Deviation in Ln scale',status)
  endif
!---
  call ftclos(unit,status)

  end
