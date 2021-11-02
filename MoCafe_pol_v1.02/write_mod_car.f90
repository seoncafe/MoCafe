module write_mod
  use define, only : wp
  use utility
  use fitsio_mod
contains
!------------------
!-- 2017-09-03, Kwang-il Seon
  subroutine write_output_car(filename,grid)
  use define
  implicit none
  character(len=*),  intent(in) :: filename
  type(grid_type),   intent(in) :: grid

!--------------------------------------------------------------
  character(len=128) :: filename1
  integer :: unit,status=0

  real(kind=wp) :: cd1_1, cd1_2, cd2_1, cd2_2
  real(kind=wp) :: crpix1, crpix2, crval1, crval2
  integer       :: equinox = 2000

  !--- header keyword for spectral image
  cd1_1  = par%dxim
  cd1_2  = 0.0_wp
  cd2_1  = 0.0_wp
  cd2_2  = par%dyim
  crpix1 = (par%nxim+1)/2.0_wp
  crpix2 = (par%nyim+1)/2.0_wp
  crval1 = 0.0_wp
  crval2 = 0.0_wp

  !--- open new fits file
  call fits_open_new(unit,trim(filename),status)

  !--- write Image for scattered light
  call fits_append_image(unit,observer%scatt,status,bitpix=par%out_bitpix)

  !--- write keywords
  call fits_put_keyword(unit,'ExeTime',  par%exetime,        'Excution Time (min)',status)
  call fits_put_keyword(unit,'Nproc',    mpar%nproc,         'No. of Threads',status)
  call fits_put_keyword(unit,'xyz_sym',  par%xyz_symmetry,   'xyz_symmetry',status)
  call fits_put_keyword(unit,'xy_per',   par%xy_periodic,    'xy_periodic',status)
  call fits_put_keyword(unit,'nphotons', par%no_photons,     'number of photons',status)
  call fits_put_keyword(unit,'taumax',   par%taumax,         'tau_max',status)
  call fits_put_keyword(unit,'Nsc_tot',  par%nscatt_tot,     'Total Nscatt/photon',status)
  call fits_put_keyword(unit,'nx',       grid%nx,            'No. of x cells',status)
  call fits_put_keyword(unit,'ny',       grid%ny,            'No. of y cells',status)
  call fits_put_keyword(unit,'nz',       grid%nz,            'No. of z cells',status)
  call fits_put_keyword(unit,'xmax',     par%xmax,           'xmax',status)
  call fits_put_keyword(unit,'ymax',     par%ymax,           'ymax',status)
  call fits_put_keyword(unit,'zmax',     par%zmax,           'zmax',status)

  !--- write keywords
  call fits_put_keyword(unit,'EQUINOX',equinox,   'Equinox of Ref. Coord.' ,status)
  call fits_put_keyword(unit,'CD1_1'  ,cd1_1  ,   'Degree / Pixel',status)
  call fits_put_keyword(unit,'CD2_1'  ,cd2_1  ,   'Degree / Pixel',status)
  call fits_put_keyword(unit,'CD1_2'  ,cd1_2  ,   'Degree / Pixel',status)
  call fits_put_keyword(unit,'CD2_2'  ,cd2_2  ,   'Degree / Pixel' ,status)
  call fits_put_keyword(unit,'CRPIX1' ,crpix1 ,   'Reference Pixel in X',status)
  call fits_put_keyword(unit,'CRPIX2' ,crpix2 ,   'Reference Pixel in Y',status)
  call fits_put_keyword(unit,'CRVAL1' ,crval1 ,   'R.A. (Degree)',status)
  call fits_put_keyword(unit,'CRVAL2' ,crval2 ,   'Dec  (Degree)',status)
  call fits_put_keyword(unit,'CTYPE1' ,'RA--TAN', 'Coordinate Type',status)
  call fits_put_keyword(unit,'CTYPE2' ,'DEC-TAN', 'Coordinate Type',status)
  call fits_put_keyword(unit,'DISTANCE', par%distance,     'Distance',status)
  call fits_put_keyword(unit,'DISTUNIT', par%distance_unit,'Distance Unit',status)
  call fits_put_keyword(unit,'EXTNAME','Scattered','J(x,y) (intensity)',status)

  !--- write Image for direct light
  call fits_append_image(unit,observer%direc,status,bitpix=par%out_bitpix)
  call fits_put_keyword(unit,'EXTNAME','Direct','J(x,y) (intensity)',status)

  !--- close fits file
  call fits_close(unit,status)

  if (par%use_stokes) then
    !--- open new fits file
    !filename1 = trim(get_base_name(filename))//'_stokes.fits.gz'
    filename1 = trim(par%base_name)//'_stokes.fits.gz'
    call fits_open_new(unit,trim(filename1),status)

    !--- write Images for Stokes parameters
    call fits_append_image(unit,observer%I,status,bitpix=par%out_bitpix)
    call fits_put_keyword(unit,'EXTNAME','Stokes_I','Stokes I image',status)

    !--- write keywords
    call fits_put_keyword(unit,'EQUINOX',equinox,   'Equinox of Ref. Coord.' ,status)
    call fits_put_keyword(unit,'CD1_1'  ,cd1_1  ,   'Degree / Pixel',status)
    call fits_put_keyword(unit,'CD2_1'  ,cd2_1  ,   'Degree / Pixel',status)
    call fits_put_keyword(unit,'CD1_2'  ,cd1_2  ,   'Degree / Pixel',status)
    call fits_put_keyword(unit,'CD2_2'  ,cd2_2  ,   'Degree / Pixel' ,status)
    call fits_put_keyword(unit,'CRPIX1' ,crpix1 ,   'Reference Pixel in X',status)
    call fits_put_keyword(unit,'CRPIX2' ,crpix2 ,   'Reference Pixel in Y',status)
    call fits_put_keyword(unit,'CRVAL1' ,crval1 ,   'R.A. (Degree)',status)
    call fits_put_keyword(unit,'CRVAL2' ,crval2 ,   'Dec  (Degree)',status)
    call fits_put_keyword(unit,'CTYPE1' ,'RA--TAN', 'Coordinate Type',status)
    call fits_put_keyword(unit,'CTYPE2' ,'DEC-TAN', 'Coordinate Type',status)
    call fits_put_keyword(unit,'DISTANCE', par%distance,     'Distance',status)
    call fits_put_keyword(unit,'DISTUNIT', par%distance_unit,'Distance Unit',status)

    call fits_append_image(unit,observer%Q,status,bitpix=par%out_bitpix)
    call fits_put_keyword(unit,'EXTNAME','Stokes_Q','Stokes Q image',status)

    call fits_append_image(unit,observer%U,status,bitpix=par%out_bitpix)
    call fits_put_keyword(unit,'EXTNAME','Stokes_U','Stokes U image',status)

    call fits_append_image(unit,observer%V,status,bitpix=par%out_bitpix)
    call fits_put_keyword(unit,'EXTNAME','Stokes_V','Stokes V image',status)

    !--- close fits file
    call fits_close(unit,status)
  endif

  if (par%sightline_tau) then
    !--- open new fits file
    filename1 = trim(par%base_name)//'_tau.fits.gz'
    call fits_open_new(unit,trim(filename1),status)

    !--- write Images for Stokes parameters
    call fits_append_image(unit,observer%tau,status,bitpix=par%out_bitpix)
    call fits_put_keyword(unit,'EXTNAME','TAU','optical depth',status)

    !--- write keywords
    call fits_put_keyword(unit,'EQUINOX',equinox,   'Equinox of Ref. Coord.' ,status)
    call fits_put_keyword(unit,'CD1_1'  ,cd1_1  ,   'Degree / Pixel',status)
    call fits_put_keyword(unit,'CD2_1'  ,cd2_1  ,   'Degree / Pixel',status)
    call fits_put_keyword(unit,'CD1_2'  ,cd1_2  ,   'Degree / Pixel',status)
    call fits_put_keyword(unit,'CD2_2'  ,cd2_2  ,   'Degree / Pixel' ,status)
    call fits_put_keyword(unit,'CRPIX1' ,crpix1 ,   'Reference Pixel in X',status)
    call fits_put_keyword(unit,'CRPIX2' ,crpix2 ,   'Reference Pixel in Y',status)
    call fits_put_keyword(unit,'CRVAL1' ,crval1 ,   'R.A. (Degree)',status)
    call fits_put_keyword(unit,'CRVAL2' ,crval2 ,   'Dec  (Degree)',status)
    call fits_put_keyword(unit,'CTYPE1' ,'RA--TAN', 'Coordinate Type',status)
    call fits_put_keyword(unit,'CTYPE2' ,'DEC-TAN', 'Coordinate Type',status)
    call fits_put_keyword(unit,'DISTANCE', par%distance,     'Distance',status)
    call fits_put_keyword(unit,'DISTUNIT', par%distance_unit,'Distance Unit',status)

    !--- close fits file
    call fits_close(unit,status)
  endif

  end subroutine write_output_car
  !=================================================
end module write_mod
