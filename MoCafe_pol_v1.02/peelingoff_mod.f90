module peelingoff_mod
  use define
  use utility
  use mathlib
  use raytrace, only : raytrace_to_edge_car
  use memory_mod
contains
!--------------------------------------------------
subroutine sightline_tau(grid)
  !--- calculate the optical depth that are projected to a detector-plane pixel (2020/08/09).
  use mpi
  implicit none
  type(grid_type),     intent(in)    :: grid
  !--- local variables
  type(photon_type) :: pobs
  real(kind=wp) :: delt(6),dist,tau
  real(kind=wp) :: kx,ky,kz,kr
  integer       :: ix,iy
  integer       :: loop,loop1,loop2,nsize
  integer       :: ierr
  integer       :: jj

  nsize = observer%nxim * observer%nyim
  call loop_divide(nsize,mpar%nproc,mpar%p_rank,loop1,loop2)
  do loop=loop1,loop2
     call array_2D_indices(observer%nxim,observer%nyim,loop,ix,iy)
     kx = tan((ix - (observer%nxim+1.0_wp)/2.0_wp) * observer%dxim/rad2deg)
     ky = tan((iy - (observer%nyim+1.0_wp)/2.0_wp) * observer%dyim/rad2deg)
     kz = -1.0_wp
     kr = sqrt(kx*kx + ky*ky + kz*kz)
     kx = kx/kr
     ky = ky/kr
     kz = kz/kr
     pobs%kx = observer%rmatrix(1,1)*kx + observer%rmatrix(2,1)*ky + observer%rmatrix(3,1)*kz
     pobs%ky = observer%rmatrix(1,2)*kx + observer%rmatrix(2,2)*ky + observer%rmatrix(3,2)*kz
     pobs%kz = observer%rmatrix(1,3)*kx + observer%rmatrix(2,3)*ky + observer%rmatrix(3,3)*kz
     if (pobs%kx == 0.0_wp) then
        delt(1) = hugest
        delt(2) = hugest
     else
        delt(1) = (grid%xmax-observer%x)/pobs%kx
        delt(2) = (grid%xmin-observer%x)/pobs%kx
     endif
     if (pobs%ky == 0.0_wp) then
        delt(3) = hugest
        delt(4) = hugest
     else
        delt(3) = (grid%ymax-observer%y)/pobs%ky
        delt(4) = (grid%ymin-observer%y)/pobs%ky
     endif
     if (pobs%kz == 0.0_wp) then
        delt(5) = hugest
        delt(6) = hugest
     else
        delt(5) = (grid%zmax-observer%z)/pobs%kz
        delt(6) = (grid%zmin-observer%z)/pobs%kz
     endif
     !-- Find the boundary where the ray touches the grid system.
     dist = hugest
     do jj=1,6
        if (delt(jj) > 0.0_wp .and. delt(jj) < hugest) then
           pobs%x     = observer%x + pobs%kx * delt(jj)
           pobs%y     = observer%y + pobs%ky * delt(jj)
           pobs%z     = observer%z + pobs%kz * delt(jj)
           pobs%icell = floor((pobs%x-grid%xmin)/grid%dx)+1
           pobs%jcell = floor((pobs%y-grid%ymin)/grid%dy)+1
           pobs%kcell = floor((pobs%z-grid%zmin)/grid%dz)+1
           if (pobs%icell >=0 .and. pobs%icell <= grid%nx+1 .and. &
               pobs%jcell >=0 .and. pobs%jcell <= grid%ny+1 .and. &
               pobs%kcell >=0 .and. pobs%kcell <= grid%nz+1) then
              if (delt(jj) < dist) dist = delt(jj)
           endif
        endif
     enddo

     if (dist < hugest) then
        pobs%x     = observer%x + pobs%kx * dist
        pobs%y     = observer%y + pobs%ky * dist
        pobs%z     = observer%z + pobs%kz * dist
        pobs%icell = floor((pobs%x-grid%xmin)/grid%dx)+1
        pobs%jcell = floor((pobs%y-grid%ymin)/grid%dy)+1
        pobs%kcell = floor((pobs%z-grid%zmin)/grid%dz)+1

        if (pobs%icell == 0) pobs%icell = 1
        if (pobs%jcell == 0) pobs%jcell = 1
        if (pobs%kcell == 0) pobs%kcell = 1
        if (pobs%icell == grid%nx+1) pobs%icell = grid%nx
        if (pobs%jcell == grid%ny+1) pobs%jcell = grid%ny
        if (pobs%kcell == grid%nz+1) pobs%kcell = grid%nz

        if (pobs%icell >=1 .and. pobs%icell <= grid%nx .and. &
            pobs%jcell >=1 .and. pobs%jcell <= grid%ny .and. &
            pobs%kcell >=1 .and. pobs%kcell <= grid%nz) then
           !--- note that there should be no limitation on tau (tau_huge) in this routine (2020.11.08).
           call raytrace_to_edge_car(pobs,grid,tau)
           observer%tau(ix,iy) = tau
        endif
     endif
  enddo

  call reduce_mem(observer%tau)
end subroutine sightline_tau
!--------------------------------------------------
!-- if mode /= 0, capture only the input photons that are injected toward the observer.
!-- update on 2022.10.07
subroutine peeling_direct_photon_nostokes(photon,grid,mode)
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  integer, intent(in), optional :: mode
  ! local variables
  type(photon_type) :: pobs
  real(kind=wp) :: r2,r,vdet(3),wgt,tau
  integer :: peeloff_mode
  integer :: ix,iy
  integer :: i

  if (present(mode)) then
     peeloff_mode = mode
  else
     peeloff_mode = 0
  endif

  pobs = photon
  if (peeloff_mode == 0) then
     pobs%kx = (observer%x-photon%x)
     pobs%ky = (observer%y-photon%y)
     pobs%kz = (observer%z-photon%z)
     r2      = pobs%kx*pobs%kx + pobs%ky*pobs%ky + pobs%kz*pobs%kz
     r       = sqrt(r2)
     pobs%kx = pobs%kx/r
     pobs%ky = pobs%ky/r
     pobs%kz = pobs%kz/r
  endif

  do i=1,3
     vdet(i) = observer%rmatrix(i,1) * pobs%kx + &
               observer%rmatrix(i,2) * pobs%ky + &
               observer%rmatrix(i,3) * pobs%kz
  enddo

  ix = floor(atan2(-vdet(1),vdet(3))*rad2deg/observer%dxim+observer%nxim/2.0_wp) + 1
  iy = floor(atan2(-vdet(2),vdet(3))*rad2deg/observer%dyim+observer%nyim/2.0_wp) + 1

  if (ix >= 1 .and. ix <= observer%nxim .and. iy >= 1 .and. iy <= observer%nyim) then
     call raytrace_to_edge_car(pobs,grid,tau)
     if (peeloff_mode == 0) then
        wgt = exp(-tau)/(fourpi*r2) * photon%wgt
     else
        wgt = exp(-tau) * photon%wgt
     endif
     observer%direc(ix,iy) = observer%direc(ix,iy) + wgt
     ! The input photon is assumed to be unpolarized.
     if (par%use_stokes) then
        observer%I(ix,iy) = observer%I(ix,iy) + wgt
     endif
  endif
end subroutine peeling_direct_photon_nostokes
!--------------------------------------------------
subroutine peeling_scattered_photon_nostokes(photon,grid)
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  ! local variables
  real(kind=wp), parameter :: three_over_16pi = 3.0_wp/(16.0_wp*pi)
  type (photon_type) :: pobs
  real(kind=wp) :: r2,r,vdet(3)
  real(kind=wp) :: cosa,wgt,peel,tau
  integer       :: ix,iy,i

  pobs    = photon
  pobs%kx = (observer%x-photon%x)
  pobs%ky = (observer%y-photon%y)
  pobs%kz = (observer%z-photon%z)
  r2      = pobs%kx*pobs%kx + pobs%ky*pobs%ky + pobs%kz*pobs%kz
  r       = sqrt(r2)
  pobs%kx = pobs%kx/r
  pobs%ky = pobs%ky/r
  pobs%kz = pobs%kz/r

  do i=1,3
     vdet(i) = observer%rmatrix(i,1) * pobs%kx + &
               observer%rmatrix(i,2) * pobs%ky + &
               observer%rmatrix(i,3) * pobs%kz
  enddo

  !--- Bin the photon into TAN image
  ix = floor(atan2(-vdet(1),vdet(3))*rad2deg/observer%dxim+observer%nxim/2.0_wp) + 1
  iy = floor(atan2(-vdet(2),vdet(3))*rad2deg/observer%dyim+observer%nyim/2.0_wp) + 1

  if (ix >= 1 .and. ix <= observer%nxim .and. iy >= 1 .and. iy <= observer%nyim) then
     call raytrace_to_edge_car(pobs,grid,tau)
     cosa = photon%kx*pobs%kx+photon%ky*pobs%ky+photon%kz*pobs%kz
     peel = (1.0_wp - par%hgg**2)/((1.0_wp + par%hgg**2)-2.0_wp*par%hgg*cosa)**1.5_wp/fourpi
     wgt  = peel/r2 * exp(-tau) * photon%wgt
     observer%scatt(ix,iy) = observer%scatt(ix,iy) + wgt
  endif
end subroutine peeling_scattered_photon_nostokes
!--------------------------------------------------
subroutine peeling_scattered_photon_stokes(photon,grid)
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  ! local variables
  type(photon_type) :: pobs
  real(kind=wp) :: r2,r
  real(kind=wp) :: cost,sint,cosp,sinp,cos2p,sin2p,cosg,sing,cos2g,sin2g
  real(kind=wp) :: S11,S12,S33,S34
  real(kind=wp) :: Q0,U0
  real(kind=wp) :: Iobs,Qobs,Uobs,Vobs,Idet,Qdet,Udet,Vdet
  real(kind=wp) :: kx,ky,kz,nx,ny
  real(kind=wp) :: tau,wgt
  integer :: ix,iy,i

  !--- Calculate the propagation vector for the peeling-off.
  pobs    = photon
  pobs%kx = (observer%x-photon%x)
  pobs%ky = (observer%y-photon%y)
  pobs%kz = (observer%z-photon%z)
  r2      = pobs%kx**2 + pobs%ky**2 + pobs%kz**2
  r       = sqrt(r2)
  pobs%kx = pobs%kx/r
  pobs%ky = pobs%ky/r
  pobs%kz = pobs%kz/r

  !--- Transform the peeling-off vector to the observer's frame.
  kx = observer%rmatrix(1,1)*pobs%kx + observer%rmatrix(1,2)*pobs%ky + observer%rmatrix(1,3)*pobs%kz
  ky = observer%rmatrix(2,1)*pobs%kx + observer%rmatrix(2,2)*pobs%ky + observer%rmatrix(2,3)*pobs%kz
  kz = observer%rmatrix(3,1)*pobs%kx + observer%rmatrix(3,2)*pobs%ky + observer%rmatrix(3,3)*pobs%kz

  !--- Location in TAN image.
  ix = floor(atan2(-kx,kz)*rad2deg/observer%dxim+observer%nxim/2.0_wp) + 1
  iy = floor(atan2(-ky,kz)*rad2deg/observer%dyim+observer%nyim/2.0_wp) + 1

  if (ix >= 1 .and. ix <= observer%nxim .and. iy >= 1 .and. iy <= observer%nyim) then
     !--- Calculate scattering angle toward the observer.
     cost = photon%kx * pobs%kx + photon%ky * pobs%ky + photon%kz * pobs%kz
     sint = sqrt(1.0_wp - cost*cost)

     !--- Calculate the reference normal vector for peeling-off.
     if (sint /= 0.0_wp) then
        pobs%nx = (photon%ky * pobs%kz - photon%kz * pobs%ky)/sint
        pobs%ny = (photon%kz * pobs%kx - photon%kx * pobs%kz)/sint
        pobs%nz = (photon%kx * pobs%ky - photon%ky * pobs%kx)/sint
     else
        pobs%nx = photon%nx
        pobs%ny = photon%ny
        pobs%nz = photon%nz
     endif

     !--- Calculate azimuthal scattering angle toward the observer.
     cosp  =   pobs%nx * photon%nx + pobs%ny * photon%ny + pobs%nz * photon%nz
     sinp  = -(pobs%nx * photon%mx + pobs%ny * photon%my + pobs%nz * photon%mz)
     cos2p = 2.0_wp*cosp*cosp - 1.0_wp
     sin2p = 2.0_wp*cosp*sinp

     !--- Calculate scattering matrix for the peeling-off vector for the scattering angles.
     call interp_eq(scatt_mat%coss,scatt_mat%S11,cost,S11)
     call interp_eq(scatt_mat%coss,scatt_mat%S12,cost,S12)
     call interp_eq(scatt_mat%coss,scatt_mat%S33,cost,S33)
     call interp_eq(scatt_mat%coss,scatt_mat%S34,cost,S34)

     !--- Calculate Stokes parameters for the peeling-off vector.
     !--- Stokes S11 is normalized by Integral(S11(cos\theta)) d(cos\theta).
     !--- Therefore, 1/2pi should be multiplied for 4pi factor.
     Q0 =  cos2p*photon%Q + sin2p*photon%U
     U0 = -sin2p*photon%Q + cos2p*photon%U

     !Iobs = ( S11*photon%I + S12*Q0  )/twopi
     !Qobs = ( S12*photon%I + S11*Q0  )/twopi
     Idet = ( S11    + S12*Q0      )/twopi
     Qobs = ( S12    + S11*Q0      )/twopi
     Uobs = ( S33*U0 + S34*photon%V)/twopi
     Vdet = (-S34*U0 + S33*photon%V)/twopi

     !--- Rotate the Stokes parameters to the detector plane.
     !--- The Stokes parameters are defined by counter-clock rotation from the North to the South.
     !--- IAU 1974 recommendation.
     !nx = observer%rmatrix(1,1)*pobs%nx + observer%rmatrix(1,2)*pobs%ny + observer%rmatrix(1,3)*pobs%nz
     !ny = observer%rmatrix(2,1)*pobs%nx + observer%rmatrix(2,2)*pobs%ny + observer%rmatrix(2,3)*pobs%nz
     !nz = observer%rmatrix(3,1)*pobs%nx + observer%rmatrix(3,2)*pobs%ny + observer%rmatrix(3,3)*pobs%nz
     !cosg  = -nx
     !sing  =  ny
     cosg  = -(observer%rmatrix(1,1)*pobs%nx + observer%rmatrix(1,2)*pobs%ny + observer%rmatrix(1,3)*pobs%nz)
     sing  =   observer%rmatrix(2,1)*pobs%nx + observer%rmatrix(2,2)*pobs%ny + observer%rmatrix(2,3)*pobs%nz
     cos2g = 2.0_wp*cosg*cosg - 1.0_wp
     sin2g = 2.0_wp*cosg*sing

     !Idet = Iobs
     Qdet =  cos2g*Qobs + sin2g*Uobs
     Udet = -sin2g*Qobs + cos2g*Uobs
     !Vdet = Vobs

     !--- Place a fraction to be peeled off to the output Stokes images.
     call raytrace_to_edge(pobs,grid,tau)
     wgt  = 1.0_wp/r2 * exp(-tau) * photon%wgt

     observer%scatt(ix,iy) = observer%scatt(ix,iy) + wgt * Idet
     observer%I(ix,iy)     = observer%I(ix,iy)     + wgt * Idet
     observer%Q(ix,iy)     = observer%Q(ix,iy)     + wgt * Qdet
     observer%U(ix,iy)     = observer%U(ix,iy)     + wgt * Udet
     observer%V(ix,iy)     = observer%V(ix,iy)     + wgt * Vdet
  endif
end subroutine peeling_scattered_photon_stokes
end module peelingoff_mod
