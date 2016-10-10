module fits_module
  use define, only : wp, dp, sp
contains
  subroutine read_2D(filename,image,status)
  implicit none
  character(len=*),           intent(in)  :: filename
  real(kind=wp), allocatable, intent(out) :: image(:,:)
  integer,                    intent(out) :: status

!--- local variables
  real(kind=sp), allocatable :: img32(:,:)
  real(kind=dp), allocatable :: img64(:,:)
  integer :: bitpix
  integer :: unit,group=1,blocksize,naxes(2),nelements,dim1
  integer :: nfound,readwrite
  logical :: simple,anyf
  integer :: nx,ny
  real(kind=wp) :: nullval = -999.9
  character(len=40) :: comm

!---------------------------------
  simple    = .true.
  readwrite = 0

  status = 0
  call ftgiou(unit,status)
  call ftopen(unit,trim(filename),readwrite,blocksize,status)
  if (status == 0) then
     call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
     nx = naxes(1)
     ny = naxes(2)
     nelements = naxes(1)*naxes(2)

     if (.not. allocated(image)) allocate(image(nx,ny))
     call ftgkyj(unit,'BITPIX',bitpix,comm,status)
     dim1 = naxes(1)

     if ((bitpix == -32 .and. wp == sp) .or. (bitpix == -64 .and. wp == dp)) then
        call ftg2de(unit,group,nullval,dim1,naxes(1),naxes(2),image,anyf,status)
     endif
     if (bitpix == -32 .and. wp == dp) then
        if (.not. allocated(img32)) allocate(img32(nx,ny))
        call ftg2de(unit,group,real(nullval,4),dim1,naxes(1),naxes(2),img32,anyf,status)
        image(:,:) = img32(:,:)
        if (allocated(img32)) deallocate(img32)
     endif
     if (bitpix == -64 .and. wp == sp) then
        if (.not. allocated(img64)) allocate(img64(nx,ny))
        call ftg2dd(unit,group,real(nullval,8),dim1,naxes(1),naxes(2),img64,anyf,status)
        image(:,:) = img64(:,:)
        if (allocated(img64)) deallocate(img64)
     endif
     call ftclos(unit,status)
  else
     write(*,*) 'Error in reading the data file :', trim(filename),' status = ',status
  endif
  call ftfiou(unit,status)

  end subroutine read_2D
!------------------
  !-- 2015-10-16, Kwangil Seon
  !-- Uncertainties for direct and scattered light are now calculated.
  subroutine write_output(filename,output,bitpix)
  use define
  implicit none
  character(len=*),  intent(in) :: filename
  type(output_type), intent(in) :: output
  integer, optional, intent(in) :: bitpix
!--------------------------------------------------------------
  integer :: equinox=2000
  real(kind=wp) :: cd1_1,cd2_1,cd1_2,cd2_2,crpix1,crpix2,crval1,crval2

  integer :: bitpix_image
  integer :: unit,status=0,blocksize=1,naxis=2,naxes(2),nelements
  logical :: simple,extend

!---------------------------------
  cd1_1   = output%dx
  cd1_2   = 0.0d0
  cd2_1   = 0.0d0
  cd2_2   = output%dy
  crpix1  = (output%nx+1)/2.0d0
  crpix2  = (output%ny+1)/2.0d0
  crval1  = 0.0d0
  crval2  = 0.0d0
!---------------------------------
  simple=.true.
  extend=.true.
  naxes(1)=output%nx
  naxes(2)=output%ny
  nelements=naxes(1)*naxes(2)
  call unlink(trim(filename))
  call ftgiou(unit,status)
  call ftinit(unit,trim(filename),blocksize,status)

  if (.not. present(bitpix)) then
     bitpix_image = -32
  else
     bitpix_image = bitpix
  endif
  call ftphpr(unit,simple,bitpix_image,naxis,naxes,0,1,extend,status)
  if (bitpix_image /= -32) then
     if (par%output_mode == 1) call ftpprd(unit,1,1,nelements,output%tot,status)
     if (par%output_mode == 2) call ftpprd(unit,1,1,nelements,output%scatt,status)
  else
     if (par%output_mode == 1) call ftppre(unit,1,1,nelements,real(output%tot,  4),status)
     if (par%output_mode == 2) call ftppre(unit,1,1,nelements,real(output%scatt,4),status)
  endif
  if (par%output_mode == 1) call ftpkys(unit,'EXTNAME','TOTAL','Total light',status)
  if (par%output_mode == 2) call ftpkys(unit,'EXTNAME','SCATTERED','Scattered light',status)

  call ftpkyj(unit,'EQUINOX',equinox,'Equinox of Ref. Coord.' ,status)
  call ftpkyd(unit,'CD1_1'  ,cd1_1  ,-8,'Degree / Pixel',status)
  call ftpkyd(unit,'CD2_1'  ,cd2_1  ,-8,'Degree / Pixel',status)
  call ftpkyd(unit,'CD1_2'  ,cd1_2  ,-8,'Degree / Pixel',status)
  call ftpkyd(unit,'CD2_2'  ,cd2_2  ,-8,'Degree / Pixel' ,status)
  call ftpkyd(unit,'CRPIX1' ,crpix1 ,-8,'Reference Pixel in X',status)
  call ftpkyd(unit,'CRPIX2' ,crpix2 ,-8,'Reference Pixel in Y',status)
  call ftpkyd(unit,'CRVAL1' ,crval1 ,-8,'R.A. (Degree)',status)
  call ftpkyd(unit,'CRVAL2' ,crval2 ,-8,'Dec  (Degree)',status)
  call ftpkys(unit,'CTYPE1' ,'RA--TAN','Coordinate Type',status)
  call ftpkys(unit,'CTYPE2' ,'DEC-TAN','Coordinate Type',status)
  call ftpkyk(unit,'NPHOTONS',par%nphotons,'Number of Input Photons',status)
  call ftpkyj(unit,'NTHREADS',par%nthreads,'Number of Threads',status)
  if (len_trim(par%psf_file) > 0) then
     call ftpkys(unit,'PSF' ,trim(par%psf_file),'Point Spread Function',status)
  else
     call ftpkys(unit,'PSF' ,'NONE','Point Spread Function',status)
  endif

!---------------------------------
  if (par%output_mode == 2) then
     simple=.true.
     extend=.true.
     naxes(1)=output%nx
     naxes(2)=output%ny
     nelements=naxes(1)*naxes(2)

     call ftiimg(unit,bitpix_image,naxis,naxes,status)
     if (bitpix_image /= -32) then
        call ftpprd(unit,1,1,nelements,output%direc,status)
     else
        call ftppre(unit,1,1,nelements,real(output%direc,4),status)
     endif
     call ftpkys(unit,'EXTNAME','DIRECT','Direct light',status)
     call ftpkyk(unit,'NPHOTONS',par%nphotons,'Number of Input Photons',status)

     call ftpkyj(unit,'EQUINOX',equinox,'Equinox of Ref. Coord.' ,status)
     call ftpkyd(unit,'CD1_1'  ,cd1_1  ,-8,'Degree / Pixel',status)
     call ftpkyd(unit,'CD2_1'  ,cd2_1  ,-8,'Degree / Pixel',status)
     call ftpkyd(unit,'CD1_2'  ,cd1_2  ,-8,'Degree / Pixel',status)
     call ftpkyd(unit,'CD2_2'  ,cd2_2  ,-8,'Degree / Pixel' ,status)
     call ftpkyd(unit,'CRPIX1' ,crpix1 ,-8,'Reference Pixel in X',status)
     call ftpkyd(unit,'CRPIX2' ,crpix2 ,-8,'Reference Pixel in Y',status)
     call ftpkyd(unit,'CRVAL1' ,crval1 ,-8,'R.A. (Degree)',status)
     call ftpkyd(unit,'CRVAL2' ,crval2 ,-8,'Dec  (Degree)',status)
     call ftpkys(unit,'CTYPE1' ,'RA--TAN','Coordinate Type',status)
     call ftpkys(unit,'CTYPE2' ,'DEC-TAN','Coordinate Type',status)
     call ftpkyk(unit,'NPHOTONS',par%nphotons,'Number of Input Photons',status)
     call ftpkyj(unit,'NTHREADS',par%nthreads,'Number of Threads',status)
     if (len_trim(par%psf_file) > 0) then
        call ftpkys(unit,'PSF' ,trim(par%psf_file),'Point Spread Function',status)
     else
        call ftpkys(unit,'PSF' ,'NONE','Point Spread Function',status)
     endif
  endif
!---------------------------------
  simple=.true.
  extend=.true.
  naxes(1)=output%nx
  naxes(2)=output%ny
  nelements=naxes(1)*naxes(2)

  call ftiimg(unit,bitpix_image,naxis,naxes,status)
  if (bitpix_image /= -32) then
     if (par%output_mode == 1) call ftpprd(unit,1,1,nelements,output%tot_sig,status)
     if (par%output_mode == 2) call ftpprd(unit,1,1,nelements,output%scatt_sig,status)
  else
     if (par%output_mode == 1) call ftppre(unit,1,1,nelements,real(output%tot_sig,4),status)
     if (par%output_mode == 2) call ftppre(unit,1,1,nelements,real(output%scatt_sig,4),status)
  endif
  if (par%output_mode ==1) call ftpkys(unit,'EXTNAME','TOT_SIG',  'Total light',status)
  if (par%output_mode ==2) call ftpkys(unit,'EXTNAME','SCATT_SIG','Scattered light',status)

  call ftpkyj(unit,'EQUINOX',equinox,'Equinox of Ref. Coord.' ,status)
  call ftpkyd(unit,'CD1_1'  ,cd1_1  ,-8,'Degree / Pixel',status)
  call ftpkyd(unit,'CD2_1'  ,cd2_1  ,-8,'Degree / Pixel',status)
  call ftpkyd(unit,'CD1_2'  ,cd1_2  ,-8,'Degree / Pixel',status)
  call ftpkyd(unit,'CD2_2'  ,cd2_2  ,-8,'Degree / Pixel' ,status)
  call ftpkyd(unit,'CRPIX1' ,crpix1 ,-8,'Reference Pixel in X',status)
  call ftpkyd(unit,'CRPIX2' ,crpix2 ,-8,'Reference Pixel in Y',status)
  call ftpkyd(unit,'CRVAL1' ,crval1 ,-8,'R.A. (Degree)',status)
  call ftpkyd(unit,'CRVAL2' ,crval2 ,-8,'Dec  (Degree)',status)
  call ftpkys(unit,'CTYPE1' ,'RA--TAN','Coordinate Type',status)
  call ftpkys(unit,'CTYPE2' ,'DEC-TAN','Coordinate Type',status)
!  call ftpkyj(unit,'NPHOTONS',par%nphotons,'Number of Input Photons',status)
  call ftpkyk(unit,'NPHOTONS',par%nphotons,'Number of Input Photons',status)
  call ftpkyj(unit,'NTHREADS',par%nthreads,'Number of Threads',status)
  if (len_trim(par%psf_file) > 0) then
     call ftpkys(unit,'PSF' ,trim(par%psf_file),'Point Spread Function',status)
  else
     call ftpkys(unit,'PSF' ,'NONE','Point Spread Function',status)
  endif
!---------------------------------
  if (par%output_mode == 2) then
     simple=.true.
     extend=.true.
     naxes(1)=output%nx
     naxes(2)=output%ny
     nelements=naxes(1)*naxes(2)

     call ftiimg(unit,bitpix_image,naxis,naxes,status)
     if (bitpix_image /= -32) then
        call ftpprd(unit,1,1,nelements,output%direc_sig,status)
     else
        call ftppre(unit,1,1,nelements,real(output%direc_sig,4),status)
     endif
     call ftpkys(unit,'EXTNAME','DIREC_SIG','Direct light',status)

     call ftpkyj(unit,'EQUINOX',equinox,'Equinox of Ref. Coord.' ,status)
     call ftpkyd(unit,'CD1_1'  ,cd1_1  ,-8,'Degree / Pixel',status)
     call ftpkyd(unit,'CD2_1'  ,cd2_1  ,-8,'Degree / Pixel',status)
     call ftpkyd(unit,'CD1_2'  ,cd1_2  ,-8,'Degree / Pixel',status)
     call ftpkyd(unit,'CD2_2'  ,cd2_2  ,-8,'Degree / Pixel' ,status)
     call ftpkyd(unit,'CRPIX1' ,crpix1 ,-8,'Reference Pixel in X',status)
     call ftpkyd(unit,'CRPIX2' ,crpix2 ,-8,'Reference Pixel in Y',status)
     call ftpkyd(unit,'CRVAL1' ,crval1 ,-8,'R.A. (Degree)',status)
     call ftpkyd(unit,'CRVAL2' ,crval2 ,-8,'Dec  (Degree)',status)
     call ftpkys(unit,'CTYPE1' ,'RA--TAN','Coordinate Type',status)
     call ftpkys(unit,'CTYPE2' ,'DEC-TAN','Coordinate Type',status)
!     call ftpkyj(unit,'NPHOTONS',par%nphotons,'Number of Input Photons',status)
     call ftpkyk(unit,'NPHOTONS',par%nphotons,'Number of Input Photons',status)
     call ftpkyj(unit,'NTHREADS',par%nthreads,'Number of Threads',status)
     if (len_trim(par%psf_file) > 0) then
        call ftpkys(unit,'PSF' ,trim(par%psf_file),'Point Spread Function',status)
     else
        call ftpkys(unit,'PSF' ,'NONE','Point Spread Function',status)
     endif
  endif
!---------------------------------

  call ftclos(unit,status)
  call ftfiou(unit,status)

  end subroutine write_output
end module fits_module
