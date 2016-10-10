module setup_idl_mod
contains
 subroutine setup_par(no_photons,nprint4,hgg,albedo,luminosity,&
                      dust1,tauface1,dust_rscale1,dust_zscale1,dust_rmax1,dust_zmax1,&
                      dust2,tauface2,dust_rscale2,dust_zscale2,dust_rmax2,dust_zmax2,&
                      rmax,pmax,zmax,nr,np,nz,&
                      disk_name,disk_rscale,disk_zscale,disk_rmax,disk_zmax,&
                      bulge_name,sersic_index,Reff,axial_ratio,BulgeToDisk,&
                      inclination_angle,position_angle,phase_angle,distance,&
                      nxim,nyim,dxim,dyim,&
                      psf_file,left_right_fold)
  use define
  implicit none
  ! input/output variables
  real(kind=wp),      intent(in) :: no_photons,hgg,albedo,luminosity
  integer,            intent(in) :: nprint4
  character(len=128), intent(in) :: dust1,dust2,disk_name,bulge_name
  real(kind=wp),      intent(in) :: tauface1,dust_rscale1,dust_zscale1,dust_rmax1,dust_zmax1,&
                                    tauface2,dust_rscale2,dust_zscale2,dust_rmax2,dust_zmax2,&
                                    disk_rscale,disk_zscale,disk_rmax,disk_zmax,&
                                    sersic_index,Reff,axial_ratio,BulgeToDisk
  real(kind=wp),      intent(in) :: rmax,pmax,zmax
  integer,            intent(in) :: nr,np,nz,nxim,nyim
  real(kind=wp),      intent(in) :: inclination_angle,position_angle,phase_angle,distance,dxim,dyim
  character(len=128), intent(in) :: psf_file
  integer,            intent(in) :: left_right_fold

  ! local variables
  character(len=128) :: distance_unit
  type(source_type)  :: source
  type(dust_type)    :: dust(MAX_dust_geometry)
  real(kind=wp)      :: hgg2

  distance_unit = 'kpc'

  par%nphotons = no_photons
  par%nprint   = nprint4
  par%hgg      = hgg
  par%albedo   = albedo
  par%luminosity = luminosity
  !par%tau_crit = tau_crit
  !par%dist_crit = dist_crit
  !par%wgt_min   = wgt_min
  !par%p_survival = p_survival
  par%inclination_angle = inclination_angle
  par%position_angle    = position_angle
  par%phase_angle       = phase_angle
  par%distance          = distance
  par%nr = nr
  !par%np = np
  par%nz = nz
  par%rmax = rmax
  par%zmax = zmax
  par%nxim = nxim
  par%nyim = nyim
  par%dxim = dxim
  par%dyim = dyim
  par%distance_unit = distance_unit
  !par%out_file = out_file
  par%psf_file = psf_file
  !if (left_right_fold == 1 .and. position_angle == 0.0_wp) then
  !   par%left_right_fold = .true.
  !else
  !   par%left_right_fold = .false.
  !endif
  ! Note output_mode = 2 in IDL-interface routine.
  par%output_mode = 2

  par%np = 1
  if (par%nprint >= no_photons) par%nprint = no_photons/10
  if (no_photons < 10)  par%nprint = 1

  dust(1)%name = trim(dust1)
  dust(2)%name = trim(dust2)
  dust(1)%tau_faceon = tauface1
  dust(2)%tau_faceon = tauface2
  dust(1)%rscale = dust_rscale1
  dust(2)%rscale = dust_rscale2
  dust(1)%zscale = dust_zscale1
  dust(2)%zscale = dust_zscale2
  dust(1)%rmax = dust_rmax1
  dust(2)%rmax = dust_rmax2
  dust(1)%zmax = dust_zmax1
  dust(2)%zmax = dust_zmax2

  source%diskname = trim(disk_name)
  source%rscale   = disk_rscale
  source%zscale   = disk_zscale
  source%rmax     = disk_rmax
  source%zmax     = disk_zmax
  source%bulgename    = trim(bulge_name)
  source%sersic_index = sersic_index
  source%Reff         = Reff
  source%axial_ratio  = axial_ratio
  source%BulgeToDisk  = BulgeToDisk

  par%source  = source
  par%dust    = dust

  hgg2     = par%hgg * par%hgg
  par%gg1  = 1.0_wp - hgg2
  par%gg2  = 1.0_wp + hgg2
  par%gg3  = 2.0_wp * par%hgg

  if (trim(source%diskname)=='exponential') then
     if (source%rmax > rmax) then
        write(*,*) ' source_rmax should be <= system rmax. Setting "source%rmax = rmax."'
        source%rmax = rmax
     endif
     if (source%zmax > zmax) then
        write(*,*) ' source_zmax should be <= system zmax. Setting "source%zmax = zmax"'
        source%zmax = zmax
     endif
  endif
  if (trim(source%bulgename)=='sersic') then
     if (source%sersic_index < 0.5_wp) then
        write(*,*) ' sersic_index of bulge component should be >= 0.5. Setting "sersic_index = 0.5"'
        source%sersic_index = 0.5_wp
     endif
     if (source%axial_ratio <= 0.0_wp) then
        write(*,*) ' axial_ratio of bulge component should be > 0.0. Setting "axial_ratio = 0.1"'
        source%axial_ratio = 0.1_wp
     endif
  endif

  end subroutine setup_par
!------------------------
  !-------------------
  subroutine setup_idl(photon,grid,observer)
  use define
  use grid_mod
  implicit none

  type(photon_type),   intent(out) :: photon
  type(grid_type),     intent(out) :: grid
  type(observer_type), intent(out) :: observer

! local variables
  integer :: i,j,k,nzcen
  real(kind=wp) :: x,y,z,r,taueq,taupole
  real(kind=wp) :: cosa,sina,cosb,sinb,cosg,sing
  integer :: id
  real(kind=wp) :: steradian_pix

  write(6,'(a)') '----------------------------------------------'
  write(6,'(a,i12)') 'Total photons: ',par%nphotons
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
  photon%lscale = (par%luminosity/par%nphotons)/(distance2cm**2)/steradian_pix
  if (par%left_right_fold .and. par%phase_angle == 0.0_wp) photon%lscale = photon%lscale/2.0_wp

  return
  end subroutine setup_idl
end module setup_idl_mod
