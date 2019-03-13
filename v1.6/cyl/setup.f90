module setup_mod
contains
  subroutine read_input
  use define
#ifdef MPI
  use mpi
#endif
  implicit none

! local_variables
  integer(kind=i8b)  :: nphotons,nprint
  character(len=128) :: distance_unit
  character(len=128) :: model_infile,out_file,psf_file

  real(kind=wp) :: no_photons
  real(kind=wp) :: luminosity
  real(kind=wp) :: inclination_angle, position_angle, phase_angle, distance

  integer :: nr,np,nz,nh,nxim,nyim
  real(kind=wp) :: dxim,dyim,xshift,yshift
  real(kind=wp) :: rmax,zmax,r_alpha,z_alpha

  integer :: i
  type(dust_type)   :: dust(MAX_dust_geometry)
  type(source_type) :: source

  real(kind=wp) :: hgg,albedo,hgg2
  real(kind=wp) :: tau_crit,dist_crit,wgt_min,p_survival
  logical       :: left_right_fold
  integer       :: output_mode

#ifdef MPI
  integer :: ierr
#endif
  integer :: myid
!--- Read in parameters from params.par using namelist command
  namelist /input_parameters/ no_photons,hgg,albedo,luminosity,distance_unit,&
                              out_file,nprint,tau_crit,wgt_min,psf_file,left_right_fold,output_mode
  namelist /geometry_parameters/ nr,np,nz,rmax,zmax,r_alpha,z_alpha,dust,source
  namelist /observer_parameters/ inclination_angle,position_angle,phase_angle,distance,nxim,nyim,dxim,dyim,xshift,yshift

! default values
!   if these variables are not included in input file, then the default values will be used.
  no_photons = 1e6
  nprint     = 1e7
  hgg        = 0.5
  albedo     = 0.5
  luminosity = 1.0
  tau_crit   = 2.0_wp
  dist_crit  = 1.0_wp
  wgt_min    = 1e-6_wp
  p_survival = 0.5_wp
  inclination_angle = 0.0_wp
  position_angle    = 0.0_wp
  phase_angle       = 0.0_wp
  distance          = 9500.0_wp
  nr   = 30
  np   = 1
  nz   = 61
  nxim = 100
  nyim = 100
  rmax = 18.0
  zmax = 18.0
  r_alpha = -999.0_wp
  z_alpha = -999.0_wp
  dxim = 0.00041666667
  dyim = 0.00041666667
  xshift = 0.0_wp
  yshift = 0.0_wp
  distance_unit = 'kpc'
  out_file  = 'out.fits.gz'
  psf_file  = ''
  left_right_fold = .true.
  output_mode     = 1

  dust(:)%name    = ''
  source%diskname = 'exponential'

  myid = 0
#ifdef MPI
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
#endif

! read model input
! Note: command_argument_count are get_command_argument are intrinsic functions of fortran 2003.
!       use iargc and getarg for old fortran.
  if (command_argument_count() >= 1) then
     call get_command_argument(1,model_infile)
  else
     model_infile = 'model.in'
  endif
  open(10,file=trim(model_infile),status='old')
  read(10,input_parameters)
  read(10,geometry_parameters)
  read(10,observer_parameters)
  close(10)

!--- dust and source geometries
  par%ndust   = count(len_trim(dust%name) > 0 .and. dust%tau_faceon > 0.0_wp)
  if (myid==0) then
     write(6,'(a)') '----------------------------------------------'
     if (par%ndust == 0) then
        write(*,'(a)') 'At least one dust geometry should be defined.'
        stop
     endif
     write(*,'(i3,a)') par%ndust,' dust geometry(ies) is(are) defined.'
  endif
  if (trim(source%diskname)=='exponential') then
     if (source%rmax > rmax) then
        if (myid==0) write(*,*) ' source_rmax should be <= system rmax. Setting "source%rmax = rmax."'
        source%rmax = rmax
     endif
     if (source%zmax > zmax) then
        if (myid==0) write(*,*) ' source_zmax should be <= system zmax. Setting "source%zmax = zmax"'
        source%zmax = zmax
     endif
  endif
  if (trim(source%bulgename)=='sersic') then
     if (source%sersic_index < 0.5_wp) then
        if (myid==0) write(*,*) ' sersic_index of bulge should be >= 0.5. Setting "sersic_index = 0.5"'
        source%sersic_index = 0.5_wp
     endif
     if (source%axial_ratio <= 0.0_wp) then
        if (myid==0) write(*,*) ' axial_ratio of bulge should be > 0.0. Setting "axial_ratio = 0.1"'
        source%axial_ratio = 0.1_wp
     endif
  endif
  par%source = source
  par%dust   = dust

  np = 1
  if (nprint >= no_photons) nprint = no_photons/10
  if (no_photons < 10)      nprint = 1

  par%nphotons = no_photons
  par%nprint   = nprint
  par%hgg      = hgg
  par%albedo   = albedo
  par%luminosity   = luminosity
  par%tau_crit     = tau_crit
  par%dist_crit    = dist_crit
  par%wgt_min      = wgt_min
  par%p_survival   = p_survival
  par%inclination_angle = inclination_angle
  par%position_angle    = position_angle
  par%phase_angle       = phase_angle
  par%distance          = distance
  par%distance_unit     = distance_unit
  par%nr   = nr
  par%np   = np
  par%nz   = nz
  par%rmax = rmax
  par%zmax = zmax
  par%r_alpha = r_alpha
  par%z_alpha = z_alpha
  par%nxim = nxim
  par%nyim = nyim
  par%dxim = dxim
  par%dyim = dyim
  par%xshift   = xshift
  par%yshift   = yshift
  par%out_file = out_file
  par%psf_file = psf_file
  par%left_right_fold = left_right_fold
  par%output_mode     = output_mode

  hgg2     = par%hgg * par%hgg
  par%gg1  = 1.0_wp - hgg2
  par%gg2  = 1.0_wp + hgg2
  par%gg3  = 2.0_wp * par%hgg

  if (myid==0) write(6,'(a,i12)') 'Total photons: ',par%nphotons

  end subroutine read_input
  !-------------------
  !-------------------
  subroutine setup(photon,grid,observer)
  use define
  use grid_mod
  implicit none

  type(photon_type),   intent(out) :: photon
  type(grid_type),     intent(out) :: grid
  type(observer_type), intent(out) :: observer

! local variables
  real(kind=wp) :: cosa,sina,cosb,sinb,cosg,sing
  real(kind=wp) :: steradian_pix

! read input parameters
  call read_input
  call grid_create(grid)

!--- setup observer and image plane
! alpha = phase angle
! beta  = inclination angle
! gamma = - posisiont angle (note a negative sign).
!           Position angle is measured counter-clock wise from X axis of the detector plane.
  cosa = cos(par%phase_angle*deg2rad)
  sina = sin(par%phase_angle*deg2rad)
  cosb = cos(par%inclination_angle*deg2rad)
  sinb = sin(par%inclination_angle*deg2rad)
  cosg = cos(-par%position_angle*deg2rad)
  sing = sin(-par%position_angle*deg2rad)
  observer%rmatrix(1,1) =  cosa*cosg - sina*cosb*sing
  observer%rmatrix(2,1) = -cosa*sing - sina*cosb*cosg
  observer%rmatrix(3,1) =  sina*sinb
  observer%rmatrix(1,2) =  sina*cosg + cosa*cosb*sing
  observer%rmatrix(2,2) = -sina*sing + cosa*cosb*cosg
  observer%rmatrix(3,2) = -cosa*sinb
  observer%rmatrix(1,3) =  sinb*sing
  observer%rmatrix(2,3) =  sinb*cosg
  observer%rmatrix(3,3) =  cosb

! observer's coordinates in grid system.
  observer%x =  par%distance*sina*sinb
  observer%y = -par%distance*cosa*sinb
  observer%z =  par%distance*cosb

!--- setup distance scale
  select case(trim(par%distance_unit))
     case ('kpc')
        distance2cm = kpc2cm
     case ('au')
        distance2cm = au2cm
     case default
        distance2cm = pc2cm
  end select

!--- setup sources luminosity scale
  steradian_pix = (par%dxim*par%dyim)*(deg2rad**2)
  !photon%lscale = (par%luminosity/dble(par%nphotons))/(distance2cm**2)/steradian_pix
  photon%lscale = (par%luminosity/(dble(par%nphotons)/dble(par%nthreads))/(distance2cm**2)/steradian_pix
  if (par%left_right_fold .and. par%phase_angle == 0.0_wp) photon%lscale = photon%lscale/2.0_wp

  return
  end subroutine setup
!---------------------------------------------------------------------------------
end module setup_mod
