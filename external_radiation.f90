module external_radiation
  use define
  use random
  use mathlib
  use raytrace, only : raytrace_to_edge_car

  !--- external radiation distribution vs. cos(theta), theta = angle from z-axis.
  type ISRF_type
     integer :: ncost
     logical :: initialized = .false.
     real(kind=wp), allocatable :: cost(:)
     real(kind=wp), allocatable :: PDF(:)
  end type ISRF_type
  type(ISRF_type) :: ISRF

  !--- procedure
  procedure(), pointer :: peeling_direct_external_cyl => peeling_direct_external_cyl1
  !--- external_illumination_sph1 is faster than external_illumination_sph2
  procedure(), pointer :: external_illumination_sph   => external_illumination_sph1
  !procedure(), pointer :: external_illumination_sph   => external_illumination_sph2
  ! peeling_direct_external_sph1 is not approriate for external_illumination_sph2.
  !procedure(), pointer :: peeling_direct_external_sph => peeling_direct_external_sph1
  procedure(), pointer :: peeling_direct_external_sph => peeling_direct_external_sph2
  procedure(), pointer :: peeling_direct_external_rec => peeling_direct_external_rec1
contains
!--------------------------------------------------
!-- 2025.03.16/2024.05.03/2023.10.14
!--------------------------------------------------
subroutine setup_radiation_angular_PDF(fname)
  use define
  use memory_mod
  use mpi
  implicit none
  character(len=128) :: fname
  integer :: unit, info, ncost
  real(kind=wp) :: cost1, PDF1, PDF_sum

  if (.not.ISRF%initialized) then
     if (.not. allocated(ISRF%cost)) allocate(ISRF%cost(0))
     if (.not. allocated(ISRF%PDF))  allocate(ISRF%PDF(0))
     open(newunit=unit,file=trim(fname),status='old')
     do while(.true.)
        read(unit,*,iostat=info) cost1, PDF1
        if (info /= 0) exit
        ISRF%cost = [ISRF%cost, cost1]
        ISRF%PDF  = [ISRF%PDF,  PDF1]
     enddo
     close(unit)
     !-- Note that the sum of PDF should be 2.
     ncost   = size(ISRF%cost)
     PDF_sum = 0.5_wp * sum((ISRF%PDF(2:ncost)+ISRF%PDF(1:ncost-1)) * abs(ISRF%cost(2:ncost)-ISRF%cost(1:ncost-1)))
     ISRF%PDF(:)      = ISRF%PDF(:)/PDF_sum * 2.0_wp
     ISRF%initialized = .true.
  endif
end subroutine setup_radiation_angular_PDF
!--------------------------------------------------
function rand_radiation_angle_PP() result(cost)
  use define
  implicit none
  real(wp)      :: cost
  logical, save :: rand_radiation_PP_init = .false.
  real(wp), allocatable, save :: cosx(:), Px(:), Palias(:)
  integer,  allocatable, save :: alias(:)
  !$OMP THREADPRIVATE(rand_radiation_PP_init,cosx,Px,Palias,alias)
  integer :: nx, i, unit, info, ncost
  real(kind=wp) :: cost1, P1, Psum

  if (.not.rand_radiation_PP_init) then
     if (.not. allocated(cosx)) allocate(cosx(0))
     if (.not. allocated(Px))   allocate(Px(0))
     open(newunit=unit,file=trim(par%radiation_angular_PDF_file),status='old')
     do while(.true.)
        read(unit,*,iostat=info) cost1, P1
        if (info /= 0) exit
        cosx = [cosx, cost1]
        Px   = [Px,   P1]
     enddo
     close(unit)
     !-- calculate the normalization factor.
     nx   = size(cosx)
     !Psum = 0.0_wp
     !do i=2, nx
     !   Psum = Psum + 0.5_wp * (Px(i) + Px(i-1))*abs(cosx(i) - cosx(i-1))
     !enddo
     Psum  = 0.5_wp * sum( (Px(2:nx)+Px(1:nx-1)) * abs(cosx(2:nx)-cosx(1:nx-1)) )
     Px(:) = Px(:) / Psum

     !-- setup alias
     if (.not. allocated(Palias)) allocate(Palias(nx-1))
     if (.not. allocated(alias))  allocate(alias(nx-1))
     do i=1, nx-1
        Palias(i) = (Px(i) + Px(i+1))/2.0_wp
     enddo
     Palias(:) = Palias(:)/sum(Palias)
     call random_alias_setup(Palias, alias)
     rand_radiation_PP_init = .true.

     !+++ Note that peeling_direct_exteranl_sph1 is not appropriate for external_illumination_sph2.
     peeling_direct_external_sph => peeling_direct_external_sph2
  endif

  cost = rand_alias_linear(Palias, alias, cosx, Px)
end function rand_radiation_angle_PP
!--------------------------------------------------
subroutine external_illumination_sph1(photon,grid)
   implicit none
   type(photon_type), intent(inout) :: photon
   type(grid_type),   intent(in)    :: grid
   real(kind=wp) :: cost,sint,phi,cosp,sinp
   real(kind=wp) :: kx,ky,kz,mx,my,mz,nx,ny,nz
   !--- random number for incident location
   cost = 2.0_wp*rand_number() - 1.0_wp
   sint = sqrt(1.0_wp-cost*cost)
   phi  = twopi*rand_number()
   cosp = cos(phi)
   sinp = sin(phi)

   !--- reference vectors for isotropically incident ray
   kx = sint*cosp
   ky = sint*sinp
   kz = cost
   mx = cost*cosp
   my = cost*sinp
   mz = -sint
   nx = -sinp
   ny =  cosp
   nz =  0.0_wp

   !--- position of the incident photon
   !--- please note that par%rmax should be positive
   photon%x = par%rmax*kx
   photon%y = par%rmax*ky
   photon%z = par%rmax*kz

   !--- random direction for propagation direction
   !--- note that pi/2 < theta < pi, i.e., -1 < cost < 0.
   cost = -sqrt(rand_number())
   sint = sqrt(1.0_wp-cost*cost)
   phi  = twopi*rand_number()
   cosp = cos(phi)
   sinp = sin(phi)

   photon%kx = sint * (cosp * mx + sinp * nx) + cost * kx
   photon%ky = sint * (cosp * my + sinp * ny) + cost * ky
   photon%kz = sint * (cosp * mz + sinp * nz) + cost * kz

   photon%wgt = 1.0_wp

   if (len_trim(par%radiation_angular_PDF_file) > 0) then
      if (.not. ISRF%initialized) call setup_radiation_angular_PDF(par%radiation_angular_PDF_file)
      call interp_eq(ISRF%cost,ISRF%PDF,photon%kz,photon%wgt)
   endif

   if (par%use_stokes) then
     !--- Find theta and phi in the fixed lab frame.
     cost = photon%kz
     sint = sqrt(1.0_wp - cost*cost)
     cosp = photon%kx/sint
     sinp = photon%ky/sint

     !--- Set the reference normal perpendicular to the propagation direction.
     photon%mx =  cost * cosp
     photon%my =  cost * sinp
     photon%mz = -sint
     photon%nx = -sinp
     photon%ny =  cosp
     photon%nz =  0.0_wp

     !--- Set the Stokes parameters (assume unpolarized light)
     photon%I = photon%wgt
     photon%Q = 0.0_wp
     photon%U = 0.0_wp
     photon%V = 0.0_wp
   endif
end subroutine external_illumination_sph1
!--------------------------------------------------
subroutine external_illumination_sph2(photon,grid)
   implicit none
   type(photon_type), intent(inout) :: photon
   type(grid_type),   intent(in)    :: grid
   real(wp) :: cost,sint,phi,cosp,sinp
   integer  :: jj, j0, icell, jcell, kcell
   real(wp) :: delt(6),x0,y0,z0,r0,dist,x,y,z
   real(wp) :: kx,ky,kz

   !--- random number for photon direction
   if (len_trim(par%radiation_angular_PDF_file) == 0) then
      cost = 2.0_wp*rand_number() - 1.0_wp
   else
      cost = rand_radiation_angle_PP()
      call setup_radiation_angular_PDF(par%radiation_angular_PDF_file)
   endif
   sint = sqrt(1.0_wp-cost*cost)
   phi  = twopi*rand_number()
   cosp = cos(phi)
   sinp = sin(phi)

   !--- Set the propagation vector and the reference normal vectors perpendicular to the propagation direction.
   photon%kx =  sint * cosp
   photon%ky =  sint * sinp
   photon%kz =  cost
   photon%mx =  cost * cosp
   photon%my =  cost * sinp
   photon%mz = -sint
   photon%nx = -sinp
   photon%ny =  cosp
   photon%nz =  0.0_wp

   !--- photon weight
   photon%wgt = 1.0_wp

   !--- starting position of the incident photon
   !--- first, find a point on the surface, which is perpendicular to the photon direction and passes through the origin..
   phi  = twopi*rand_number()
   cosp = cos(phi)
   sinp = sin(phi)

   r0 = par%rmax * sqrt(rand_number())
   x0 = r0 * (cosp * photon%mx + sinp * photon%nx)
   y0 = r0 * (cosp * photon%my + sinp * photon%ny)
   z0 = r0 * (cosp * photon%mz + sinp * photon%nz)

   !-- Find the starting boundary where the ray touches the grid system.
   kx = -photon%kx
   ky = -photon%ky
   kz = -photon%kz
   if (kx == 0.0_wp) then
      delt(1) = hugest
      delt(2) = hugest
   else
      delt(1) = (grid%xmax-x0)/kx
      delt(2) = (grid%xmin-x0)/kx
   endif
   if (ky == 0.0_wp) then
      delt(3) = hugest
      delt(4) = hugest
   else
      delt(3) = (grid%ymax-y0)/ky
      delt(4) = (grid%ymin-y0)/ky
   endif
   if (kz == 0.0_wp) then
      delt(5) = hugest
      delt(6) = hugest
   else
      delt(5) = (grid%zmax-z0)/kz
      delt(6) = (grid%zmin-z0)/kz
   endif

   j0   = minloc(delt, mask = delt >= 0.0_wp, dim=1)
   dist = delt(j0)

   photon%x     = x0 + kx * dist
   photon%y     = y0 + ky * dist
   photon%z     = z0 + kz * dist
   photon%icell = floor((photon%x-grid%xmin)/grid%dx)+1
   photon%jcell = floor((photon%y-grid%ymin)/grid%dy)+1
   photon%kcell = floor((photon%z-grid%zmin)/grid%dz)+1

   !if (j0 == 1) then
   !   photon%icell = grid%nx+1
   !   photon%x     = grid%xface(grid%nx+1)
   !else if (j0 == 2) then
   !   photon%icell = 1
   !   photon%x     = grid%xface(1)
   !else if (j0 == 3) then
   !   photon%jcell = grid%ny+1
   !   photon%y     = grid%yface(grid%ny+1)
   !else if (j0 == 4) then
   !   photon%jcell = 1
   !   photon%y     = grid%yface(1)
   !else if (j0 == 5) then
   !   photon%kcell = grid%nz+1
   !   photon%z     = grid%zface(grid%nz+1)
   !else if (j0 == 6) then
   !   photon%kcell = 1
   !   photon%z     = grid%zface(1)
   !endif
   if (dist == delt(1)) then
      photon%icell = grid%nx+1
      photon%x     = grid%xface(grid%nx+1)
   endif
   if (dist == delt(2)) then
      photon%icell = 1
      photon%x     = grid%xface(1)
   endif
   if (dist == delt(3)) then
      photon%jcell = grid%ny+1
      photon%y     = grid%yface(grid%ny+1)
   endif
   if (dist == delt(4)) then
      photon%jcell = 1
      photon%y     = grid%yface(1)
   endif
   if (dist == delt(5)) then
      photon%kcell = grid%nz+1
      photon%z     = grid%zface(grid%nz+1)
   endif
   if (dist == delt(6)) then
      photon%kcell = 1
      photon%z     = grid%zface(1)
   endif

   if (photon%icell == grid%nx+1) photon%icell = grid%nx
   if (photon%jcell == grid%ny+1) photon%jcell = grid%ny
   if (photon%kcell == grid%nz+1) photon%kcell = grid%nz

   if (par%use_stokes) then
      !--- Set the Stokes parameters (assume unpolarized light)
      photon%I = photon%wgt
      photon%Q = 0.0_wp
      photon%U = 0.0_wp
      photon%V = 0.0_wp
   endif
end subroutine external_illumination_sph2
!--------------------------------------------------
subroutine external_illumination_cyl(photon,grid)
   implicit none
   type(photon_type), intent(inout) :: photon
   type(grid_type),   intent(in)    :: grid
   real(kind=wp) :: cost,sint,phi,cosp,sinp,rr
   real(kind=wp) :: area1,area2,frac
   real(kind=wp) :: kx,ky,kz,mx,my,mz,nx,ny,nz
   !---
   area1 = 2.0_wp*pi*par%rmax**2
   area2 = 4.0_wp*pi*par%rmax*(2.0_wp*par%zmax)
   frac  = area1/(area1 + area2)

   !--- random number for incident location
   if (rand_number() < frac) then
      rr  = sqrt(rand_number())*par%rmax
      phi = twopi*rand_number()
      photon%x = rr*cos(phi)
      photon%y = rr*sin(phi)

      kx = 1.0
      ky = 0.0
      kz = 0.0
      mx = 0.0
      my = 1.0
      mz = 0.0
      nx = 0.0
      ny = 0.0
      select case(floor(2.0*rand_number()))
      case (0)
         photon%z = grid%xface(grid%nx+1)
         nz       = 1.0
      case (1)
         photon%z = grid%xface(1)
         nz       = -1.0
      end select
   else
      phi      = twopi*rand_number()
      cosp     = cos(phi)
      sinp     = sin(phi)
      photon%x = par%rmax * cosp
      photon%y = par%rmax * sinp
      photon%z = grid%zrange * rand_number() + grid%zmin

      !--- reference vectors for isotropically incident ray
      !--- theta = pi/2
      kx = cosp
      ky = sinp
      kz = 0.0_wp
      mx = 0.0_wp
      my = 0.0_wp
      mz = -1.0_wp
      nx = -sinp
      ny =  cosp
      nz =  0.0_wp
   endif

   !--- random direction for propagation direction
   !--- note that pi/2 < theta < pi, i.e., -1 < cost < 0.
   cost = -sqrt(rand_number())
   sint = sqrt(1.0_wp-cost*cost)
   phi  = twopi*rand_number()
   cosp = cos(phi)
   sinp = sin(phi)

   photon%kx = sint * (cosp * mx + sinp * nx) + cost * kx
   photon%ky = sint * (cosp * my + sinp * ny) + cost * ky
   photon%kz = sint * (cosp * mz + sinp * nz) + cost * kz

   photon%wgt = 1.0_wp

   if (par%use_stokes) then
     !--- Find theta and phi in the fixed lab frame.
     cost = photon%kz
     sint = sqrt(1.0_wp - cost*cost)
     cosp = photon%kx/sint
     sinp = photon%ky/sint

     !--- Set the reference normal perpendicular to the propagation direction.
     photon%mx =  cost * cosp
     photon%my =  cost * sinp
     photon%mz = -sint
     photon%nx = -sinp
     photon%ny =  cosp
     photon%nz =  0.0_wp

     !--- Set the Stokes parameters (assume unpolarized light)
     photon%I = 1.0_wp
     photon%Q = 0.0_wp
     photon%U = 0.0_wp
     photon%V = 0.0_wp
   endif
end subroutine external_illumination_cyl
!--------------------------------------------------
subroutine external_illumination_rec(photon,grid)
   implicit none
   type(photon_type), intent(inout) :: photon
   type(grid_type),   intent(in)    :: grid
   real(kind=wp) :: cost,sint,phi,cosp,sinp
   real(kind=wp) :: kx,ky,kz,mx,my,mz,nx,ny,nz
   !--- random direction for propagation direction
   !--- note that pi/2 < theta < pi, i.e., -1 < cost < 0.
   cost = -sqrt(rand_number())
   sint = sqrt(1.0_wp-cost*cost)
   phi  = twopi*rand_number()
   cosp = cos(phi)
   sinp = sin(phi)

   kx = 0.0
   ky = 0.0
   kz = 0.0
   mx = 0.0
   my = 0.0
   mz = 0.0
   nx = 0.0
   ny = 0.0
   nz = 0.0
   !--- random number for incident location
   select case(floor(6.0*rand_number() + 1.0))
   case (1)
      photon%x = grid%xface(grid%nx+1)
      photon%y = grid%yrange * rand_number() + grid%ymin
      photon%z = grid%zrange * rand_number() + grid%zmin
      kx = 1.0
      my = 1.0
      nz = 1.0
   case (2)
      photon%x = grid%xface(1)
      photon%y = grid%yrange * rand_number() + grid%ymin
      photon%z = grid%zrange * rand_number() + grid%zmin
      kx = -1.0
      my = 1.0
      nz = 1.0
   case (3)
      photon%y = grid%yface(grid%ny+1)
      photon%z = grid%zrange * rand_number() + grid%zmin
      photon%x = grid%xrange * rand_number() + grid%xmin
      ky = 1.0
      mx = 1.0
      nz = 1.0
   case (4)
      photon%y = grid%yface(1)
      photon%z = grid%zrange * rand_number() + grid%zmin
      photon%x = grid%xrange * rand_number() + grid%xmin
      ky = -1.0
      mx = 1.0
      nz = 1.0
   case (5)
      photon%z = grid%zface(grid%nz+1)
      photon%x = grid%xrange * rand_number() + grid%xmin
      photon%y = grid%yrange * rand_number() + grid%ymin
      kz = 1.0
      mx = 1.0
      ny = 1.0
   case (6)
      photon%z = grid%zface(1)
      photon%x = grid%xrange * rand_number() + grid%xmin
      photon%y = grid%yrange * rand_number() + grid%ymin
      kz = -1.0
      mx = 1.0
      ny = 1.0
   end select

   photon%kx = sint * (cosp * mx + sinp * nx) + cost * kx
   photon%ky = sint * (cosp * my + sinp * ny) + cost * ky
   photon%kz = sint * (cosp * mz + sinp * nz) + cost * kz

   photon%wgt = 1.0_wp

   if (par%use_stokes) then
     !--- Find theta and phi in the fixed lab frame.
     cost = photon%kz
     sint = sqrt(1.0_wp - cost*cost)
     cosp = photon%kx/sint
     sinp = photon%ky/sint

     !--- Set the reference normal perpendicular to the propagation direction.
     photon%mx =  cost * cosp
     photon%my =  cost * sinp
     photon%mz = -sint
     photon%nx = -sinp
     photon%ny =  cosp
     photon%nz =  0.0_wp

     !--- Set the Stokes parameters (assume unpolarized light)
     photon%I = 1.0_wp
     photon%Q = 0.0_wp
     photon%U = 0.0_wp
     photon%V = 0.0_wp
   endif
end subroutine external_illumination_rec
!--------------------------------------------------
subroutine peeling_direct_external_sph1(photon,grid)
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  !-- local variables
  type(photon_type) :: pobs
  real(kind=wp) :: r2,r,kx,ky,kz,wgt,wgt0,tau,scal
  real(kind=wp) :: kx0,ky0,kz0,kr0
  real(kind=wp) :: cost
  integer :: ix,iy
  integer :: k

  do k=1, par%nobs
    !-- Define the peeled-off photon
    pobs = photon

    !--- Photon propagation direction.
    pobs%kx = (observer(k)%x-pobs%x)
    pobs%ky = (observer(k)%y-pobs%y)
    pobs%kz = (observer(k)%z-pobs%z)
    r2      = pobs%kx**2 + pobs%ky**2 + pobs%kz**2
    r       = sqrt(r2)
    pobs%kx = pobs%kx/r
    pobs%ky = pobs%ky/r
    pobs%kz = pobs%kz/r

    !-- (kx0, ky0, kz0) = normal vector perpendicular to the surface at the photon position.
    kx0 = pobs%x
    ky0 = pobs%y
    kz0 = pobs%z
    kr0 = sqrt(kx0**2 + ky0**2 + kz0**2)
    kx0 = kx0/kr0
    ky0 = ky0/kr0
    kz0 = kz0/kr0

    !-- angle between the normal vector and the propagation direction.
    cost = pobs%kx * kx0 + pobs%ky * ky0 + pobs%kz * kz0

    !-- if the angle is less than pi/2, then ignore the photon.
    if (cost > 0.0) cycle

    !-- Weighting factor when the isotropic illumination is assumed.
    !-- Normalized distribution function for mu and phi (-1 < mu < 0 and 0 < phi < 2pi).
    pobs%wgt = abs(cost)/pi

    !--- weighting factor to take into account the radiation angular distribution.
    if (len_trim(par%radiation_angular_PDF_file) > 0) then
       call interp_eq(ISRF%cost,ISRF%PDF,pobs%kz,scal)
       pobs%wgt = pobs%wgt * scal
    endif

    !-- Transform the peeling-off vector to the observer's frame.
    kx = observer(k)%rmatrix(1,1)*pobs%kx + observer(k)%rmatrix(1,2)*pobs%ky + observer(k)%rmatrix(1,3)*pobs%kz
    ky = observer(k)%rmatrix(2,1)*pobs%kx + observer(k)%rmatrix(2,2)*pobs%ky + observer(k)%rmatrix(2,3)*pobs%kz
    kz = observer(k)%rmatrix(3,1)*pobs%kx + observer(k)%rmatrix(3,2)*pobs%ky + observer(k)%rmatrix(3,3)*pobs%kz

    !-- Location in the TAN image.
    ix = floor(atan2(-kx,kz)*rad2deg/observer(k)%dxim+observer(k)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(k)%dyim+observer(k)%nyim/2.0_wp) + 1

    if (ix >= 1 .and. ix <= observer(k)%nxim .and. iy >= 1 .and. iy <= observer(k)%nyim) then
       call raytrace_to_edge_car(pobs,grid,tau)
       wgt0 = 1.0_wp/r2 * pobs%wgt
       wgt  = exp(-tau) * wgt0
       observer(k)%direc(ix,iy) = observer(k)%direc(ix,iy) + wgt
       if (par%save_direc0) then
          observer(k)%direc0(ix,iy) = observer(k)%direc0(ix,iy) + wgt0
       endif
       ! The input photon is assumed to be unpolarized.
       if (par%use_stokes) then
          observer(k)%I(ix,iy) = observer(k)%I(ix,iy) + wgt
       endif
    endif
  enddo
end subroutine peeling_direct_external_sph1
!--------------------------------------------------
subroutine peeling_direct_external_sph2(photon,grid)
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  ! local variables
  type(photon_type) :: pobs
  real(kind=wp) :: r2,r,kx,ky,kz,wgt,wgt0,tau,scal
  real(kind=wp) :: kx0,ky0,kz0,kr0,xx,yy,zz
  real(kind=wp) :: cost, cosvt0, cosvt, sinvt, vphi, cosvp, sinvp, lum_fac
  integer :: ix,iy
  integer :: k

  !-- Assume the isotropic intensity.
  cost   = -sqrt(rand_number())

  cosvt0 = par%rmax / par%distance
  cosvt  = cost*sqrt(1.0_wp - cosvt0**2 + (cosvt0*cost)**2) + cosvt0*(1.0_wp - cost**2)
  sinvt  = sqrt(1.0_wp - cosvt**2)
  vphi   = twopi*rand_number()
  cosvp  = cos(vphi)
  sinvp  = sin(vphi)
  lum_fac = (1.0+cosvt0)/(twopi*3.0*cosvt0**2)*( 2.0*cosvt0**2 - (1.0-cosvt0)*(1.0 - sqrt(1.0-cosvt0**2)) )
  !--- Note that lum_fac = 1/(4pi) as distance -> infinity. (comment added on 2024.04.30)
  !print*,lum_fac

  do k=1, par%nobs
    !-- (kx0, ky0, kz0) = unit vector connecting the cloud center and observer.
    !                     the vector defining a lightcone toward the observer.
    kx0 = observer(k)%x
    ky0 = observer(k)%y
    kz0 = observer(k)%z
    kr0 = sqrt(kx0**2 + ky0**2 + kz0**2)
    kx0 = kx0/kr0
    ky0 = ky0/kr0
    kz0 = kz0/kr0

    !-- (xx, yy, zz) = a direction vector of a point on the cloud surface
    !--                defined by (vartheta, varphi) about the cloud-observer sightline.
    if (abs(kz0) >= 0.99999999999_wp) then
       xx  = sinvt*cosvp
       yy  = sinvt*sinvp
       zz  = cosvt
    else
       kr0 = sqrt(kx0**2 + ky0**2)
       xx  = cosvt*kx0 + sinvt*(kz0*kx0*cosvp - ky0*sinvp)/kr0
       yy  = cosvt*ky0 + sinvt*(kz0*ky0*cosvp + kx0*sinvp)/kr0
       zz  = cosvt*kz0 - sinvt*cosvp*kr0
    endif

    !-- Define the peeled-off photon
    pobs     = photon
    pobs%wgt = lum_fac

    !-- pobs%(x, y, z) = photon location on the cloud surface.
    pobs%x = par%rmax * xx
    pobs%y = par%rmax * yy
    pobs%z = par%rmax * zz

    !--- Photon propagation direction.
    pobs%kx = (observer(k)%x-pobs%x)
    pobs%ky = (observer(k)%y-pobs%y)
    pobs%kz = (observer(k)%z-pobs%z)
    r2      = pobs%kx*pobs%kx + pobs%ky*pobs%ky + pobs%kz*pobs%kz
    r       = sqrt(r2)
    pobs%kx = pobs%kx/r
    pobs%ky = pobs%ky/r
    pobs%kz = pobs%kz/r

    !--- weighting factor to take into account the radiation angular distribution.
    if (len_trim(par%radiation_angular_PDF_file) > 0) then
       call interp_eq(ISRF%cost,ISRF%PDF,pobs%kz,scal)
       pobs%wgt = pobs%wgt * scal
    endif

    !--- Cell Index
    pobs%icell = floor((pobs%x-grid%xmin)/grid%dx)+1
    pobs%jcell = floor((pobs%y-grid%ymin)/grid%dy)+1
    pobs%kcell = floor((pobs%z-grid%zmin)/grid%dz)+1
    !write(*,*) pobs%x,pobs%y,pobs%z,pobs%icell,pobs%jcell,pobs%kcell

    !-- Transform the peeling-off vector to the observer's frame.
    kx = observer(k)%rmatrix(1,1)*pobs%kx + observer(k)%rmatrix(1,2)*pobs%ky + observer(k)%rmatrix(1,3)*pobs%kz
    ky = observer(k)%rmatrix(2,1)*pobs%kx + observer(k)%rmatrix(2,2)*pobs%ky + observer(k)%rmatrix(2,3)*pobs%kz
    kz = observer(k)%rmatrix(3,1)*pobs%kx + observer(k)%rmatrix(3,2)*pobs%ky + observer(k)%rmatrix(3,3)*pobs%kz

    !-- Location in the TAN image.
    ix = floor(atan2(-kx,kz)*rad2deg/observer(k)%dxim+observer(k)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(k)%dyim+observer(k)%nyim/2.0_wp) + 1

    if (ix >= 1 .and. ix <= observer(k)%nxim .and. iy >= 1 .and. iy <= observer(k)%nyim) then
       call raytrace_to_edge_car(pobs,grid,tau)
       wgt0 = 1.0_wp/r2 * pobs%wgt
       wgt  = exp(-tau) * wgt0
       observer(k)%direc(ix,iy) = observer(k)%direc(ix,iy) + wgt
       if (par%save_direc0) then
          observer(k)%direc0(ix,iy) = observer(k)%direc0(ix,iy) + wgt0
       endif
       ! The input photon is assumed to be unpolarized.
       if (par%use_stokes) then
          observer(k)%I(ix,iy) = observer(k)%I(ix,iy) + wgt
       endif
    endif
  enddo
end subroutine peeling_direct_external_sph2
!--------------------------------------------------
subroutine peeling_direct_external_rec1(photon,grid)
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  !-- local variables
  type(photon_type) :: pobs
  real(kind=wp) :: r2,r,kx,ky,kz,wgt,wgt0,tau
  real(kind=wp) :: kx0,ky0,kz0,kr0
  real(kind=wp) :: cost
  integer :: ix,iy
  integer :: k

  do k=1, par%nobs
    !-- Define the peeled-off photon
    pobs = photon

    !--- Photon propagation direction.
    pobs%kx = (observer(k)%x-pobs%x)
    pobs%ky = (observer(k)%y-pobs%y)
    pobs%kz = (observer(k)%z-pobs%z)
    r2      = pobs%kx**2 + pobs%ky**2 + pobs%kz**2
    r       = sqrt(r2)
    pobs%kx = pobs%kx/r
    pobs%ky = pobs%ky/r
    pobs%kz = pobs%kz/r

    !-- (kx0, ky0, kz0) = normal vector perpendicular to the surface at the photon position.
    kx0 = 0.0
    ky0 = 0.0
    kz0 = 0.0
    if (pobs%x == grid%xface(1))         kx0 = -1.0
    if (pobs%x == grid%xface(grid%nx+1)) kx0 =  1.0
    if (pobs%y == grid%yface(1))         ky0 = -1.0
    if (pobs%y == grid%yface(grid%ny+1)) ky0 =  1.0
    if (pobs%z == grid%zface(1))         kz0 = -1.0
    if (pobs%z == grid%zface(grid%nz+1)) kz0 =  1.0

    !-- angle between the normal vector and the propagation direction.
    cost = pobs%kx * kx0 + pobs%ky * ky0 + pobs%kz * kz0

    !-- if the angle is less than pi/2, then ignore the photon.
    if (cost > 0.0) cycle

    !-- Weighting factor when the isotropic illumination is assumed.
    !-- Normalized distribution function for mu and phi (-1 < mu < 0 and 0 < phi < 2pi).
    pobs%wgt = abs(cost)/pi

    !-- Transform the peeling-off vector to the observer's frame.
    kx = observer(k)%rmatrix(1,1)*pobs%kx + observer(k)%rmatrix(1,2)*pobs%ky + observer(k)%rmatrix(1,3)*pobs%kz
    ky = observer(k)%rmatrix(2,1)*pobs%kx + observer(k)%rmatrix(2,2)*pobs%ky + observer(k)%rmatrix(2,3)*pobs%kz
    kz = observer(k)%rmatrix(3,1)*pobs%kx + observer(k)%rmatrix(3,2)*pobs%ky + observer(k)%rmatrix(3,3)*pobs%kz

    !-- Location in the TAN image.
    ix = floor(atan2(-kx,kz)*rad2deg/observer(k)%dxim+observer(k)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(k)%dyim+observer(k)%nyim/2.0_wp) + 1

    if (ix >= 1 .and. ix <= observer(k)%nxim .and. iy >= 1 .and. iy <= observer(k)%nyim) then
       call raytrace_to_edge_car(pobs,grid,tau)
       wgt0 = 1.0_wp/r2 * pobs%wgt
       wgt  = exp(-tau) * wgt0
       observer(k)%direc(ix,iy) = observer(k)%direc(ix,iy) + wgt
       if (par%save_direc0) then
          observer(k)%direc0(ix,iy) = observer(k)%direc0(ix,iy) + wgt0
       endif
       ! The input photon is assumed to be unpolarized.
       if (par%use_stokes) then
          observer(k)%I(ix,iy) = observer(k)%I(ix,iy) + wgt
       endif
    endif
  enddo
end subroutine peeling_direct_external_rec1
!--------------------------------------------------
subroutine peeling_direct_external_cyl1(photon,grid)
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  !-- local variables
  type(photon_type) :: pobs
  real(kind=wp) :: r2,r,kx,ky,kz,wgt,wgt0,tau
  real(kind=wp) :: kx0,ky0,kz0,kr0
  real(kind=wp) :: cost
  integer :: ix,iy
  integer :: k

  do k=1, par%nobs
    !-- Define the peeled-off photon
    pobs = photon

    !--- Photon propagation direction.
    pobs%kx = (observer(k)%x-pobs%x)
    pobs%ky = (observer(k)%y-pobs%y)
    pobs%kz = (observer(k)%z-pobs%z)
    r2      = pobs%kx**2 + pobs%ky**2 + pobs%kz**2
    r       = sqrt(r2)
    pobs%kx = pobs%kx/r
    pobs%ky = pobs%ky/r
    pobs%kz = pobs%kz/r

    !-- (kx0, ky0, kz0) = normal vector perpendicular to the surface at the photon position.
    if (pobs%z /= grid%zface(1) .and. pobs%z /= grid%zface(grid%nz+1)) then
       kx0 = pobs%x
       ky0 = pobs%y
       kr0 = sqrt(kx0**2 + ky0**2)
       kx0 = kx0/kr0
       ky0 = ky0/kr0
       kz0 = 0.0_wp
    else
       kx0 = 0.0
       ky0 = 0.0
       if (pobs%z == grid%zface(1))         kz0 = -1.0
       if (pobs%z == grid%zface(grid%nz+1)) kz0 =  1.0
    endif

    !-- angle between the normal vector and the propagation direction.
    cost = pobs%kx * kx0 + pobs%ky * ky0 + pobs%kz * kz0

    !-- if the angle is less than pi/2, then ignore the photon.
    if (cost > 0.0) cycle

    !-- Weighting factor when the isotropic illumination is assumed.
    !-- Normalized distribution function for mu and phi (-1 < mu < 0 and 0 < phi < 2pi).
    pobs%wgt = abs(cost)/pi

    !-- Transform the peeling-off vector to the observer's frame.
    kx = observer(k)%rmatrix(1,1)*pobs%kx + observer(k)%rmatrix(1,2)*pobs%ky + observer(k)%rmatrix(1,3)*pobs%kz
    ky = observer(k)%rmatrix(2,1)*pobs%kx + observer(k)%rmatrix(2,2)*pobs%ky + observer(k)%rmatrix(2,3)*pobs%kz
    kz = observer(k)%rmatrix(3,1)*pobs%kx + observer(k)%rmatrix(3,2)*pobs%ky + observer(k)%rmatrix(3,3)*pobs%kz

    !-- Location in the TAN image.
    ix = floor(atan2(-kx,kz)*rad2deg/observer(k)%dxim+observer(k)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(k)%dyim+observer(k)%nyim/2.0_wp) + 1

    if (ix >= 1 .and. ix <= observer(k)%nxim .and. iy >= 1 .and. iy <= observer(k)%nyim) then
       call raytrace_to_edge_car(pobs,grid,tau)
       wgt0 = 1.0_wp/r2 * pobs%wgt
       wgt  = exp(-tau) * wgt0
       observer(k)%direc(ix,iy) = observer(k)%direc(ix,iy) + wgt
       if (par%save_direc0) then
          observer(k)%direc0(ix,iy) = observer(k)%direc0(ix,iy) + wgt0
       endif
       ! The input photon is assumed to be unpolarized.
       if (par%use_stokes) then
          observer(k)%I(ix,iy) = observer(k)%I(ix,iy) + wgt
       endif
    endif
  enddo
end subroutine peeling_direct_external_cyl1
!--------------------------------------------------
end module external_radiation
