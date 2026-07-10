module mrw_mod
!--- Modified Random Walk (MoCafe v2.00, Stage 8) for very optically thick
!--- cells (Min et al. 2009, A&A 497, 155; Robitaille 2010, A&A 520, A70).
!--- When a photon sits deep inside an optically thick cell, a normal walk
!--- takes O(tau^2) tiny scatter/absorb steps.  MRW replaces them with a single
!--- diffusion step to the surface of the largest sphere inscribed at the
!--- photon position: it samples the total path length travelled inside the
!--- sphere from the diffusion first-passage distribution, deposits that path
!--- into the J tally (Lucy energy estimator), and moves the photon to a random
!--- point on the sphere with a fresh isotropic direction.
!---
!--- Cartesian grid only; used in the Lucy energy passes (par%use_mrw).  The
!--- opacity is taken at the photon wavelength (rho*kappa*s_ext).
  use define
  use random,     only : rand_number
  use jtally_mod, only : jt_on, jt_sum
  implicit none
  private
  public :: setup_mrw, mrw_try_step, mrw_on

  logical :: mrw_on = .false.
  real(kind=wp) :: mrw_gamma = 2.0_wp        ! trigger: R0*rho*kappa > gamma
  integer, parameter :: NCDF = 200
  real(kind=wp) :: cdf_xi(NCDF), cdf_y(NCDF) ! inverse first-passage CDF: xi -> y

contains
  !---------------------------------------------------------------
  subroutine setup_mrw()
  implicit none
  integer :: i, n
  real(kind=wp) :: y, p, term
  mrw_gamma = par%mrw_gamma
  !--- tabulate the diffusion first-passage CDF P(y) = 2 Sum (-1)^(n+1) y^(n^2)
  !--- on a y grid, then store it as (xi=P) -> y for inversion by search.
  do i = 1, NCDF
     y = dble(i)/dble(NCDF+1)          ! y in (0,1)
     p = 0.0_wp
     do n = 1, 100
        term = ((-1.0_wp)**(n+1)) * y**(n*n)
        p = p + term
        if (abs(term) < 1.0e-12_wp) exit
     enddo
     p = 2.0_wp*p
     cdf_y(i)  = y
     cdf_xi(i) = max(0.0_wp, min(1.0_wp, p))
  enddo
  mrw_on = .true.
  if (mpar%p_rank == 0) write(*,'(a,f5.2,a)') &
     '--- Modified Random Walk enabled (gamma = ', mrw_gamma, ') ---'
  end subroutine setup_mrw

  !---------------------------------------------------------------
  !--- attempt one MRW step at the photon's current position.  Returns .true.
  !--- if a step was taken (photon moved, energy deposited); .false. if the
  !--- cell is not optically thick enough (caller does the normal walk).
  logical function mrw_try_step(photon, grid) result(stepped)
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(in)    :: grid
  real(kind=wp) :: asca, aabs, R0, dx0, dx1, dy0, dy1, dz0, dz1
  real(kind=wp) :: xi, y, ct, deposit, ux, uy, uz, cost, sint, phi, alb, rhk

  stepped = .false.
  if (.not. mrw_on) return
  rhk = grid%rhokap(photon%icell,photon%jcell,photon%kcell)*photon%s_ext
  if (rhk <= 0.0_wp) return
  !--- split into scattering and absorption opacity.  MRW models the diffusion
  !--- of the SCATTERING random walk; absorption decays the weight along the
  !--- path (so it works correctly in absorbing dust, where it simply rarely
  !--- triggers because the scattering optical depth stays small).
  alb  = photon%albedo
  asca = rhk*alb
  aabs = rhk*(1.0_wp - alb)

  !--- R0 = distance to the nearest cell wall.
  dx0 = photon%x - grid%xface(photon%icell);  dx1 = grid%xface(photon%icell+1) - photon%x
  dy0 = photon%y - grid%yface(photon%jcell);  dy1 = grid%yface(photon%jcell+1) - photon%y
  dz0 = photon%z - grid%zface(photon%kcell);  dz1 = grid%zface(photon%kcell+1) - photon%z
  R0  = min(dx0, dx1, dy0, dy1, dz0, dz1)
  if (R0 <= 0.0_wp) return
  if (asca*R0 <= mrw_gamma) return           ! not scattering-thick: normal walk

  !--- sample the first-passage variable y from the tabulated CDF.
  xi = rand_number()
  if (xi <= cdf_xi(1)) then
     y = cdf_y(1)
  else if (xi >= cdf_xi(NCDF)) then
     y = cdf_y(NCDF)
  else
     block
       integer :: lo2, hi2, mid2
       lo2 = 1;  hi2 = NCDF
       do while (hi2 - lo2 > 1)
          mid2 = (lo2+hi2)/2
          if (cdf_xi(mid2) < xi) then;  lo2 = mid2;  else;  hi2 = mid2;  endif
       enddo
       y = cdf_y(lo2) + (cdf_y(hi2)-cdf_y(lo2))*(xi-cdf_xi(lo2))/(cdf_xi(hi2)-cdf_xi(lo2))
     end block
  endif

  !--- total scattering path length inside the inscribed sphere (Robitaille
  !--- 2010): ct = -ln(y) * 3 * asca * (R0/pi)^2.  The J-tally deposit is the
  !--- absorption-weighted path integral int_0^ct wgt*exp(-aabs*l) dl, and the
  !--- weight decays by exp(-aabs*ct) over the diffusion path.
  ct = -log(max(y, tinest)) * 3.0_wp*asca*(R0/pi)**2
  if (aabs*ct > 1.0e-12_wp) then
     deposit = photon%wgt*photon%Lpacket*(1.0_wp - exp(-aabs*ct))/aabs
  else
     deposit = photon%wgt*photon%Lpacket*ct
  endif
  if (jt_on) jt_sum(photon%il,photon%icell,photon%jcell,photon%kcell) = &
       jt_sum(photon%il,photon%icell,photon%jcell,photon%kcell) + deposit
  photon%wgt = photon%wgt * exp(-aabs*ct)

  !--- move the photon to a random point on the sphere of radius R0 and give it
  !--- a fresh isotropic propagation direction.
  cost = 2.0_wp*rand_number()-1.0_wp;  sint = sqrt(1.0_wp-cost*cost);  phi = twopi*rand_number()
  ux = sint*cos(phi);  uy = sint*sin(phi);  uz = cost
  photon%x = photon%x + R0*ux
  photon%y = photon%y + R0*uy
  photon%z = photon%z + R0*uz
  photon%icell = min(max(floor((photon%x-grid%xmin)/grid%dx)+1, 1), grid%nx)
  photon%jcell = min(max(floor((photon%y-grid%ymin)/grid%dy)+1, 1), grid%ny)
  photon%kcell = min(max(floor((photon%z-grid%zmin)/grid%dz)+1, 1), grid%nz)
  cost = 2.0_wp*rand_number()-1.0_wp;  sint = sqrt(1.0_wp-cost*cost);  phi = twopi*rand_number()
  photon%kx = sint*cos(phi);  photon%ky = sint*sin(phi);  photon%kz = cost
  photon%nscatt = photon%nscatt + 1
  stepped = .true.
  end function mrw_try_step

end module mrw_mod
