module peelingoff_mod
  use define
  use utility
  use mathlib
  use memory_mod
  use scan_mod
  !--- Peel-off integrates optical depth to the observer via the
  !--- raytrace_to_edge procedure pointer (bound in setup_procedure):
  !--- raytrace_to_edge_car for the Cartesian grid, raytrace_to_edge_clump for
  !--- the clump medium.  Using the pointer (rather than calling
  !--- raytrace_to_edge_car directly) lets the clump grid reuse these routines
  !--- unchanged; the default binding keeps the Cartesian path bit-identical.
contains
!--------------------------------------------------
!-- update on 2022.10.07
subroutine peeling_direct_photon_nostokes(photon,grid)
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  ! local variables
  type(photon_type) :: pobs
  real(kind=wp) :: r2,r,vdet(3),wgt,wgt0,tau
  integer :: ix,iy
  integer :: i, k

  do k=1, par%nobs
     pobs = photon
     pobs%kx = (observer(k)%x-photon%x)
     pobs%ky = (observer(k)%y-photon%y)
     pobs%kz = (observer(k)%z-photon%z)
     r2      = pobs%kx*pobs%kx + pobs%ky*pobs%ky + pobs%kz*pobs%kz
     r       = sqrt(r2)
     pobs%kx = pobs%kx/r
     pobs%ky = pobs%ky/r
     pobs%kz = pobs%kz/r

     do i=1,3
        vdet(i) = observer(k)%rmatrix(i,1) * pobs%kx + &
                  observer(k)%rmatrix(i,2) * pobs%ky + &
                  observer(k)%rmatrix(i,3) * pobs%kz
     enddo

     ix = floor(atan2(-vdet(1),vdet(3))*rad2deg/observer(k)%dxim+observer(k)%nxim/2.0_wp) + 1
     iy = floor(atan2(-vdet(2),vdet(3))*rad2deg/observer(k)%dyim+observer(k)%nyim/2.0_wp) + 1

     if (ix >= 1 .and. ix <= observer(k)%nxim .and. iy >= 1 .and. iy <= observer(k)%nyim) then
        call raytrace_to_edge(pobs,grid,tau)
        wgt0 = 1.0_wp/(fourpi*r2) * photon%wgt
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
  integer       :: ix,iy,i,k

  do k=1, par%nobs
     pobs    = photon
     pobs%kx = (observer(k)%x-photon%x)
     pobs%ky = (observer(k)%y-photon%y)
     pobs%kz = (observer(k)%z-photon%z)
     r2      = pobs%kx*pobs%kx + pobs%ky*pobs%ky + pobs%kz*pobs%kz
     r       = sqrt(r2)
     pobs%kx = pobs%kx/r
     pobs%ky = pobs%ky/r
     pobs%kz = pobs%kz/r

     do i=1,3
        vdet(i) = observer(k)%rmatrix(i,1) * pobs%kx + &
                  observer(k)%rmatrix(i,2) * pobs%ky + &
                  observer(k)%rmatrix(i,3) * pobs%kz
     enddo

     !--- Bin the photon into TAN image
     ix = floor(atan2(-vdet(1),vdet(3))*rad2deg/observer(k)%dxim+observer(k)%nxim/2.0_wp) + 1
     iy = floor(atan2(-vdet(2),vdet(3))*rad2deg/observer(k)%dyim+observer(k)%nyim/2.0_wp) + 1

     if (ix >= 1 .and. ix <= observer(k)%nxim .and. iy >= 1 .and. iy <= observer(k)%nyim) then
        call raytrace_to_edge(pobs,grid,tau)
        cosa = photon%kx*pobs%kx+photon%ky*pobs%ky+photon%kz*pobs%kz
        peel = (1.0_wp - par%hgg**2)/((1.0_wp + par%hgg**2)-2.0_wp*par%hgg*cosa)**1.5_wp/fourpi
        wgt  = peel/r2 * exp(-tau) * photon%wgt
        observer(k)%scatt(ix,iy) = observer(k)%scatt(ix,iy) + wgt
     endif
  enddo
end subroutine peeling_scattered_photon_nostokes
!--------------------------------------------------
subroutine peeling_scattered_photon_nostokes_scan(photon,grid)
  implicit none
  !--- (albedo, asymmetry-factor) scan peel-off.  Mirrors
  !--- peeling_scattered_photon_nostokes, but the expensive parts (tau to the
  !--- observer, 1/r^2, pixel binning) are computed ONCE and reused across the
  !--- whole n_a x n_g grid; only the closed-form factors vary inside the loop.
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  ! local variables
  type (photon_type) :: pobs
  real(kind=wp) :: r2,r,vdet(3)
  real(kind=wp) :: cosa,base,peelg,tau
  integer       :: ix,iy,i,k,ia,ig

  do k=1, par%nobs
     pobs    = photon
     pobs%kx = (observer(k)%x-photon%x)
     pobs%ky = (observer(k)%y-photon%y)
     pobs%kz = (observer(k)%z-photon%z)
     r2      = pobs%kx*pobs%kx + pobs%ky*pobs%ky + pobs%kz*pobs%kz
     r       = sqrt(r2)
     pobs%kx = pobs%kx/r
     pobs%ky = pobs%ky/r
     pobs%kz = pobs%kz/r

     do i=1,3
        vdet(i) = observer(k)%rmatrix(i,1) * pobs%kx + &
                  observer(k)%rmatrix(i,2) * pobs%ky + &
                  observer(k)%rmatrix(i,3) * pobs%kz
     enddo

     ix = floor(atan2(-vdet(1),vdet(3))*rad2deg/observer(k)%dxim+observer(k)%nxim/2.0_wp) + 1
     iy = floor(atan2(-vdet(2),vdet(3))*rad2deg/observer(k)%dyim+observer(k)%nyim/2.0_wp) + 1

     if (ix >= 1 .and. ix <= observer(k)%nxim .and. iy >= 1 .and. iy <= observer(k)%nyim) then
        call raytrace_to_edge(pobs,grid,tau)
        cosa = photon%kx*pobs%kx + photon%ky*pobs%ky + photon%kz*pobs%kz
        !--- (a,g)-independent part, computed once.  photon%wgt holds the base
        !--- (forced-first-scatter) weight only — no albedo factor in scan mode.
        base = exp(-tau)/r2 * photon%wgt
        do ig = 1, scan_ng
           peelg = hg_kernel(cosa, scan_glist(ig)) / fourpi   ! Phi(theta_obs | g_ig)
           do ia = 1, scan_na
              observer(k)%scatt_ag(ix,iy,ia,ig) = observer(k)%scatt_ag(ix,iy,ia,ig) &
                   + base * peelg * scan_Aalb(ia) * scan_Wg(ig)
           enddo
        enddo
     endif
  enddo
end subroutine peeling_scattered_photon_nostokes_scan
!--------------------------------------------------
subroutine peeling_direct_photon_nostokes_tau(photon,grid)
  implicit none
  !--- Polychromatic (tau) scan variant of peeling_direct_photon_nostokes.
  !--- The direct beam is attenuated per optical-depth scaling s_t, so the
  !--- direct image gains a tau axis: direc_t(ix,iy,it) += exp(-s_t*tau)*wgt0.
  !--- The unattenuated direc0 (if saved) stays 2-D (tau-independent).
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  ! local variables
  type(photon_type) :: pobs
  real(kind=wp) :: r2,r,vdet(3),wgt0,tau
  integer :: ix,iy,i,k,it

  do k=1, par%nobs
     pobs    = photon
     pobs%kx = (observer(k)%x-photon%x)
     pobs%ky = (observer(k)%y-photon%y)
     pobs%kz = (observer(k)%z-photon%z)
     r2      = pobs%kx*pobs%kx + pobs%ky*pobs%ky + pobs%kz*pobs%kz
     r       = sqrt(r2)
     pobs%kx = pobs%kx/r
     pobs%ky = pobs%ky/r
     pobs%kz = pobs%kz/r

     do i=1,3
        vdet(i) = observer(k)%rmatrix(i,1) * pobs%kx + &
                  observer(k)%rmatrix(i,2) * pobs%ky + &
                  observer(k)%rmatrix(i,3) * pobs%kz
     enddo

     ix = floor(atan2(-vdet(1),vdet(3))*rad2deg/observer(k)%dxim+observer(k)%nxim/2.0_wp) + 1
     iy = floor(atan2(-vdet(2),vdet(3))*rad2deg/observer(k)%dyim+observer(k)%nyim/2.0_wp) + 1

     if (ix >= 1 .and. ix <= observer(k)%nxim .and. iy >= 1 .and. iy <= observer(k)%nyim) then
        call raytrace_to_edge(pobs,grid,tau)
        wgt0 = 1.0_wp/(fourpi*r2) * photon%wgt
        do it = 1, scan_nt
           observer(k)%direc_t(ix,iy,it) = observer(k)%direc_t(ix,iy,it) + exp(-scan_s(it)*tau) * wgt0
        enddo
        if (par%save_direc0) then
           observer(k)%direc0(ix,iy) = observer(k)%direc0(ix,iy) + wgt0
        endif
     endif
  enddo
end subroutine peeling_direct_photon_nostokes_tau
!--------------------------------------------------
subroutine peeling_scattered_photon_nostokes_scan_tau(photon,grid)
  implicit none
  !--- (albedo, asymmetry-factor, optical-depth) scan peel-off.  Mirrors
  !--- peeling_scattered_photon_nostokes_scan, but the optical depth to the
  !--- observer scales with s_t, so the attenuation is moved out of "base" into
  !--- the tau loop.  The contribution factorizes into a (a,g,tau) outer product;
  !--- the expensive parts (raytrace, 1/r^2, pixel binning) are computed ONCE.
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  ! local variables
  type (photon_type) :: pobs
  real(kind=wp) :: r2,r,vdet(3)
  real(kind=wp) :: cosa,base0,vg,atten,tau
  integer       :: ix,iy,i,k,ia,ig,it

  do k=1, par%nobs
     pobs    = photon
     pobs%kx = (observer(k)%x-photon%x)
     pobs%ky = (observer(k)%y-photon%y)
     pobs%kz = (observer(k)%z-photon%z)
     r2      = pobs%kx*pobs%kx + pobs%ky*pobs%ky + pobs%kz*pobs%kz
     r       = sqrt(r2)
     pobs%kx = pobs%kx/r
     pobs%ky = pobs%ky/r
     pobs%kz = pobs%kz/r

     do i=1,3
        vdet(i) = observer(k)%rmatrix(i,1) * pobs%kx + &
                  observer(k)%rmatrix(i,2) * pobs%ky + &
                  observer(k)%rmatrix(i,3) * pobs%kz
     enddo

     ix = floor(atan2(-vdet(1),vdet(3))*rad2deg/observer(k)%dxim+observer(k)%nxim/2.0_wp) + 1
     iy = floor(atan2(-vdet(2),vdet(3))*rad2deg/observer(k)%dyim+observer(k)%nyim/2.0_wp) + 1

     if (ix >= 1 .and. ix <= observer(k)%nxim .and. iy >= 1 .and. iy <= observer(k)%nyim) then
        call raytrace_to_edge(pobs,grid,tau)
        cosa  = photon%kx*pobs%kx + photon%ky*pobs%ky + photon%kz*pobs%kz
        !--- (a,g,tau)-independent part, computed once.  No attenuation here:
        !--- the optical depth to the observer is s_t-scaled inside the tau loop.
        base0 = photon%wgt / r2
        do it = 1, scan_nt
           atten = exp(-scan_s(it)*tau) * scan_Wtau(it)
           do ig = 1, scan_ng
              vg = hg_kernel(cosa, scan_glist(ig)) / fourpi * scan_Wg(ig)   ! Phi(theta_obs|g)*W_g
              do ia = 1, scan_na
                 observer(k)%scatt_agt(ix,iy,ia,ig,it) = observer(k)%scatt_agt(ix,iy,ia,ig,it) &
                      + base0 * scan_Aalb(ia) * vg * atten
              enddo
           enddo
        enddo
     endif
  enddo
end subroutine peeling_scattered_photon_nostokes_scan_tau
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
  integer :: ix,iy,i,k

  do k=1,par%nobs
     !--- Calculate the propagation vector for the peeling-off.
     pobs    = photon
     pobs%kx = (observer(k)%x-photon%x)
     pobs%ky = (observer(k)%y-photon%y)
     pobs%kz = (observer(k)%z-photon%z)
     r2      = pobs%kx**2 + pobs%ky**2 + pobs%kz**2
     r       = sqrt(r2)
     pobs%kx = pobs%kx/r
     pobs%ky = pobs%ky/r
     pobs%kz = pobs%kz/r

     !--- Transform the peeling-off vector to the observer's frame.
     kx = observer(k)%rmatrix(1,1)*pobs%kx + observer(k)%rmatrix(1,2)*pobs%ky + observer(k)%rmatrix(1,3)*pobs%kz
     ky = observer(k)%rmatrix(2,1)*pobs%kx + observer(k)%rmatrix(2,2)*pobs%ky + observer(k)%rmatrix(2,3)*pobs%kz
     kz = observer(k)%rmatrix(3,1)*pobs%kx + observer(k)%rmatrix(3,2)*pobs%ky + observer(k)%rmatrix(3,3)*pobs%kz

     !--- Location in TAN image.
     ix = floor(atan2(-kx,kz)*rad2deg/observer(k)%dxim+observer(k)%nxim/2.0_wp) + 1
     iy = floor(atan2(-ky,kz)*rad2deg/observer(k)%dyim+observer(k)%nyim/2.0_wp) + 1

     if (ix >= 1 .and. ix <= observer(k)%nxim .and. iy >= 1 .and. iy <= observer(k)%nyim) then
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
        cosg  = -(observer(k)%rmatrix(1,1)*pobs%nx + observer(k)%rmatrix(1,2)*pobs%ny + observer(k)%rmatrix(1,3)*pobs%nz)
        sing  =   observer(k)%rmatrix(2,1)*pobs%nx + observer(k)%rmatrix(2,2)*pobs%ny + observer(k)%rmatrix(2,3)*pobs%nz
        cos2g = 2.0_wp*cosg*cosg - 1.0_wp
        sin2g = 2.0_wp*cosg*sing

        !Idet = Iobs
        Qdet =  cos2g*Qobs + sin2g*Uobs
        Udet = -sin2g*Qobs + cos2g*Uobs
        !Vdet = Vobs

        !--- Place a fraction to be peeled off to the output Stokes images.
        call raytrace_to_edge(pobs,grid,tau)
        wgt  = 1.0_wp/r2 * exp(-tau) * photon%wgt

        observer(k)%scatt(ix,iy) = observer(k)%scatt(ix,iy) + wgt * Idet
        observer(k)%I(ix,iy)     = observer(k)%I(ix,iy)     + wgt * Idet
        observer(k)%Q(ix,iy)     = observer(k)%Q(ix,iy)     + wgt * Qdet
        observer(k)%U(ix,iy)     = observer(k)%U(ix,iy)     + wgt * Udet
        observer(k)%V(ix,iy)     = observer(k)%V(ix,iy)     + wgt * Vdet
     endif
  enddo
end subroutine peeling_scattered_photon_stokes
end module peelingoff_mod
