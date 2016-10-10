  subroutine setup_idl(no_photons,hgg_in,albedo_in,luminosity,&
                       dust1,tauface1,dust_rscale1,dust_zscale1,dust_rmax1,dust_zmax1,&
                       dust2,tauface2,dust_rscale2,dust_zscale2,dust_rmax2,dust_zmax2,&
                       xmax,ymax,zmax,nx,ny,nz,&
                       source_name,source_rscale,source_zscale,source_rmax,source_zmax,&
                       obs_angles,obs_distance,nxim,nyim,dxim,dyim,&
                       photon,grid,source,observer,output,nphotons)

  use define
  use random
  implicit none
! input/output variables
  real(kind=wp), intent(in) :: no_photons,hgg_in,albedo_in,luminosity
  character(len=100), intent(in) :: dust1,dust2,source_name
  real(kind=wp), intent(in) :: tauface1,dust_rscale1,dust_zscale1,dust_rmax1,dust_zmax1,&
                               tauface2,dust_rscale2,dust_zscale2,dust_rmax2,dust_zmax2,&
                               source_rscale,source_zscale,source_rmax,source_zmax
  real(kind=wp), intent(in) :: xmax,ymax,zmax
  integer, intent(in) :: nx,ny,nz,nxim,nyim
  real(kind=wp), intent(in) :: obs_angles(3),obs_distance,dxim,dyim

  type(photon_type),   intent(out) :: photon
  type(grid_type),     intent(out) :: grid
  type(observer_type), intent(out) :: observer
  type(output_type),   intent(out) :: output
  integer(kind=i8b),   intent(out) :: nphotons
  type(source_type),   intent(out) :: source

  integer, parameter :: ndust = 2
  type(dust_type) :: dust(ndust)

  integer :: i,j,k,nxcen,nycen,nzcen
  real(kind=wp) :: x,y,z,r,taueq,taupole

! local variables
  real(kind=wp) :: cosx,sinx,cosy,siny,cosz,sinz
  character(len=100) :: distance_unit
  real(kind=wp) :: opac
  integer :: id

  hgg    = hgg_in
  albedo = albedo_in

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

  source%name = trim(source_name)
  source%rscale = source_rscale
  source%zscale = source_zscale
  source%rmax = source_rmax
  source%zmax = source_zmax

  nphotons    = no_photons

  distance_unit = 'kpc'
  nscatt_max    = 20
  wgt_min       = 1.0e-6_wp

!--- setup dust property
  write(*,'(a,2f7.4)') 'Dust Parameters(a, g): ',albedo,hgg

  hgg2 = hgg * hgg
  gg1  = 1.0_wp - hgg2
  gg2  = 1.0_wp + hgg2
  gg3  = 2.0_wp * hgg

!--- setup distance scale
  select case(trim(distance_unit))
     case ('kpc')
        distance2cm = kpc2cm
     case ('au')
        distance2cm = au2cm
     case default
        distance2cm = pc2cm
  end select

!--- Set up grid faces.
  write(6,'(a)') '----------------------------------------------'
  write(6,'(a)') 'Setting up density grid and radiation sources.'
  grid%nx = nx
  grid%ny = ny
  grid%nz = nz
  grid%xmax = xmax
  grid%ymax = ymax
  grid%zmax = zmax
  grid%xmin = -xmax
  grid%ymin = -ymax
  grid%zmin = -zmax
  grid%xrange = grid%xmax - grid%xmin
  grid%yrange = grid%ymax - grid%ymin
  grid%zrange = grid%zmax - grid%zmin
  grid%dx     = grid%xrange/grid%nx
  grid%dy     = grid%yrange/grid%ny
  grid%dz     = grid%zrange/grid%nz
  nullify(grid%xface)
  nullify(grid%yface)
  nullify(grid%zface)
  if (.not. associated(grid%xface)) allocate(grid%xface(grid%nx+1))
  if (.not. associated(grid%yface)) allocate(grid%yface(grid%ny+1))
  if (.not. associated(grid%zface)) allocate(grid%zface(grid%nz+1))
  grid%xface(:) = (/ ((i-1)*grid%dx + grid%xmin, i=1,grid%nx+1) /)
  grid%yface(:) = (/ ((i-1)*grid%dy + grid%ymin, i=1,grid%ny+1) /)
  grid%zface(:) = (/ ((i-1)*grid%dz + grid%zmin, i=1,grid%nz+1) /)

!--- setup observer and image plane
  cosz = cos(obs_angles(1)*deg2rad)
  sinz = sin(obs_angles(1)*deg2rad)
  cosy = cos(obs_angles(2)*deg2rad)
  siny = sin(obs_angles(2)*deg2rad)
  cosx = cos(obs_angles(3)*deg2rad)
  sinx = sin(obs_angles(3)*deg2rad)
  observer%tmatrix(1,1) = cosy*cosz
  observer%tmatrix(2,1) = sinx*siny*cosz-cosx*sinz
  observer%tmatrix(3,1) = cosx*siny*cosz+sinx*sinz
  observer%tmatrix(1,2) = cosy*sinz
  observer%tmatrix(2,2) = sinx*siny*sinz+cosx*cosz
  observer%tmatrix(3,2) = cosx*siny*sinz-sinx*cosz
  observer%tmatrix(1,3) = -siny
  observer%tmatrix(2,3) = sinx*cosy
  observer%tmatrix(3,3) = cosx*cosy
!! Transpose the matrix
  observer%tmatrix = transpose(observer%tmatrix)

! observer's coordinates
  observer%x =  obs_distance*cosy*cosz
  observer%y =  obs_distance*cosy*sinz
  observer%z = -obs_distance*siny

  output%nx = nxim
  output%ny = nyim
  output%dx = dxim
  output%dy = dyim
!  output%dx = atan(grid%xmax/obs_distance)/(output%nx/2.0_wp) * rad2deg
!  output%dy = atan(grid%xmax/obs_distance)/(output%ny/2.0_wp) * rad2deg
  output%steradian_pix = (output%dx*output%dy)*(deg2rad**2)

!--- setup sources luminosity scale
  photon%lscale = (luminosity/nphotons)/(distance2cm**2)/output%steradian_pix

!--- setup dust density
  nullify(grid%opacity)
  if (.not. associated(grid%opacity)) allocate(grid%opacity(grid%nx,grid%ny,grid%nz))
  grid%opacity(:,:,:) = 0.0_wp

  !$OMP PARALLELDO &
  !$OMP private(i,j,k,x,y,z,r,opac) shared(grid,dust) &
  !$OMP schedule(static,1)
  do i=1,grid%nx
     x = (grid%xface(i)+grid%xface(i+1))/2.0_wp
     do j=1,grid%ny
        y = (grid%yface(j)+grid%yface(j+1))/2.0_wp
        r = sqrt(x*x + y*y)
        do k=1,grid%nz
           z = (grid%zface(k)+grid%zface(k+1))/2.0_wp
           do id=1,ndust
             call dust_density(dust(id),r,z,opac)
             grid%opacity(i,j,k) = grid%opacity(i,j,k) + opac
           enddo
        enddo
     enddo
  enddo
  !$OMP END PARALLELDO

!--- Calculate equatorial and polar optical depths
  nxcen = (grid%nx+1)/2
  nycen = (grid%ny+1)/2
  nzcen = (grid%nz+1)/2
  taueq   = sum(grid%opacity(:,nycen,nzcen))*2.0_wp*grid%xmax/grid%ny
  taupole = sum(grid%opacity(nxcen,nycen,:))*2.0_wp*grid%zmax/grid%nz
  write(6,'(a,f7.3)') 'tau_equat  : ',taueq
  write(6,'(a,f7.3)') 'tau_pole   : ',taupole

  return
  end
!---------------------------------------------------------------------------------
  subroutine setup_pointers(output)
  use define
  implicit none

  type(output_type), intent(out) :: output

! The pointer should be allocated after the statement of
!    private(output)
!--- setup output array
  nullify(output%scatt)
  nullify(output%direc)
  if (.not. associated(output%scatt)) allocate(output%scatt(output%nx,output%ny))
  if (.not. associated(output%direc)) allocate(output%direc(output%nx,output%ny))
  output%scatt(:,:) = 0.0_wp
  output%direc(:,:) = 0.0_wp

  return
  end
