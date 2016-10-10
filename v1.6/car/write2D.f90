module write2D
contains
  !-----------------------------------------------
  subroutine write_2D(fname,array,bitpix)

  use define
  implicit none
  character(len=*),  intent(in) :: fname
  real(kind=wp),     intent(in) :: array(:,:)
  integer, optional, intent(in) :: bitpix

  integer :: status,unit,blocksize,bitpix_image,naxis,naxes(2),group,fpixel,nelements
  logical :: simple,extend
  real(kind=dp) :: cd1_1,cd2_1,cd1_2,cd2_2,crpix1,crpix2,crval1,crval2

  status    = 0
  blocksize = 1
  simple    = .true.
  naxis     = 2
  naxes(1)  = size(array,1)
  naxes(2)  = size(array,2)
  extend    = .true.
  group     = 1
  fpixel    = 1
  nelements = naxes(1)*naxes(2)
  if (.not. present(bitpix)) then
     bitpix_image = -32
  else
     bitpix_image = bitpix
     if (bitpix /= -32 .or. bitpix /= -64) bitpix_image = -64
  endif

  cd1_1   = par%dxim
  cd1_2   = 0.0d0
  cd2_1   = 0.0d0
  cd2_2   = par%dyim
  crpix1  = (par%nxim+1)/2.0d0
  crpix2  = (par%nyim+1)/2.0d0
  crval1  = 0.0d0
  crval2  = 0.0d0

  call unlink(trim(fname))
  call ftgiou(unit,status)
  call ftinit(unit,trim(fname),blocksize,status)
  call ftphpr(unit,simple,bitpix_image,naxis,naxes,0,1,extend,status)
  if (bitpix_image /= -32) then
     call ftpprd(unit,group,fpixel,nelements,real(array,8),status)
  else
     call ftppre(unit,group,fpixel,nelements,real(array,4),status)
  endif
  call ftpkys(unit,'EXTNAME','IMAGE','IMAGE',status)
  call ftpkys(unit,'CTYPE1' ,'RA--TAN','Coordinate Type',status)
  call ftpkys(unit,'CTYPE2' ,'DEC-TAN','Coordinate Type',status)

  call ftpkyd(unit,'CD1_1'  ,cd1_1  ,-8,'Degree / Pixel',status)
  call ftpkyd(unit,'CD2_1'  ,cd2_1  ,-8,'Degree / Pixel',status)
  call ftpkyd(unit,'CD1_2'  ,cd1_2  ,-8,'Degree / Pixel',status)
  call ftpkyd(unit,'CD2_2'  ,cd2_2  ,-8,'Degree / Pixel' ,status)
  call ftpkyd(unit,'CRPIX1' ,crpix1 ,-8,'Reference Pixel in X',status)
  call ftpkyd(unit,'CRPIX2' ,crpix2 ,-8,'Reference Pixel in Y',status)
  call ftpkyd(unit,'CRVAL1' ,crval1 ,-8,'R.A. (Degree)',status)
  call ftpkyd(unit,'CRVAL2' ,crval2 ,-8,'Dec  (Degree)',status)
  call ftclos(unit,status)
  call ftfiou(unit,status)

  end subroutine write_2D
end module write2D
