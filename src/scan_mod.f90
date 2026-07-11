module scan_mod
!------------------------------------------------------------------------------
! Single-run (albedo, asymmetry-factor) scan.
!
! Implements the reweighting algorithm of Seon 2010, PKAS, 25, 177:
! photons are tracked with one simulated asymmetry factor g0 (= par%hgg) and an
! immortal/continuous-absorption forward pass, while per-(a,g) weights are
! accumulated so that one simulation yields the scattered image for a whole grid
! of albedos and asymmetry factors at once:
!
!     scatt(x, y, a, g)     a = scan_alist(1:scan_na)   (default 0.1..1.0)
!                           g = scan_glist(1:scan_ng)   (default 0.0..0.9)
!
! Per scattering event k:
!   albedo:  Aalb(ia) <- a_ia    * Aalb(ia)        => Aalb(ia) = a_ia^k
!   g     :  Wg(ig)   <- [Phi(mu|g_ig)/Phi(mu|g0)] * Wg(ig)   (paper eq. 2)
! The Aalb update happens before the peel-off (so the peel at scattering k sees
! a^k); the Wg update happens after the actual scattering cosine mu is sampled
! (so the peel at scattering k sees the previous k-1 actual scatters).
!
! MoCafe is pure MPI with a single photon in flight per rank, so the
! photon accumulators below are safe as module-level state (no OpenMP-within-photon).
!------------------------------------------------------------------------------
  use define
  use utility, only : is_finite
  implicit none
  public

  !--- scan grids (filled by scan_setup)
  integer       :: scan_na = 0
  integer       :: scan_ng = 0
  real(kind=wp) :: scan_g0 = 0.0_wp
  real(kind=wp) :: scan_alist(MAX_SCAN) = 0.0_wp
  real(kind=wp) :: scan_glist(MAX_SCAN) = 0.0_wp

  !--- polychromatic optical-depth (tau) scan grids (Jonsson 2006, eqs. 29-32)
  integer       :: scan_nt = 0
  real(kind=wp) :: scan_taumax_ref = 0.0_wp        ! simulated reference taumax (s = 1)
  real(kind=wp) :: scan_tlist(MAX_SCAN) = 0.0_wp   ! target taumax values
  real(kind=wp) :: scan_s(MAX_SCAN)     = 1.0_wp   ! s_t = tlist / taumax_ref

  !--- accumulators for each photon (reset at each photon birth)
  real(kind=wp) :: scan_Aalb(MAX_SCAN) = 1.0_wp
  real(kind=wp) :: scan_Wg(MAX_SCAN)   = 1.0_wp
  real(kind=wp) :: scan_Wtau(MAX_SCAN) = 1.0_wp

  !--- carrier for each segment: the free-path optical depth tau_k just drawn for the
  !--- segment leading to the current scattering (set by run_simulation, read by
  !--- the scan scatter routine).  Only meaningful when par%use_tau_list = .true.
  real(kind=wp) :: scan_tau_step = 0.0_wp

contains
  !---------------------------------------------------------------------------
  !--- 4*pi * Phi(mu|g) for the Henyey-Greenstein phase function.
  !--- The 4*pi factor cancels in the Wg ratio; for peel-off divide by fourpi.
  pure function hg_kernel(mu, g) result(p)
    implicit none
    real(kind=wp), intent(in) :: mu, g
    real(kind=wp) :: p
    p = (1.0_wp - g*g) / (1.0_wp + g*g - 2.0_wp*g*mu)**1.5_wp
  end function hg_kernel

  !---------------------------------------------------------------------------
  subroutine scan_setup()
    implicit none
    integer :: i
    logical :: have_ref
    real(kind=wp), parameter :: stol = 1.0e-6_wp
    real(kind=wp) :: tau_ref

    scan_g0 = par%hgg

    !--- scan reference optical depth: follow the grid-normalization choice
    !--- (taumax takes precedence, else tauhomo), so tau_list is expressed in
    !--- the same measure that normalizes the medium (cf. grid_create).
    if (par%taumax > 0.0_wp) then
       tau_ref = par%taumax
    else if (par%tauhomo > 0.0_wp) then
       tau_ref = par%tauhomo
    else
       tau_ref = 0.0_wp
    endif

    !--- albedo / asymmetry-factor lists
    if (par%use_ag_list) then
       !--- albedo list: use finite entries verbatim, else default 0.1..1.0
       scan_na = count(is_finite(par%albedo_list))
       if (scan_na > 0) then
          scan_alist(1:scan_na) = pack(par%albedo_list, is_finite(par%albedo_list))
       else
          scan_na = 10
          do i = 1, scan_na
             scan_alist(i) = 0.1_wp * real(i, wp)        ! 0.1, 0.2, ..., 1.0
          enddo
       endif

       !--- asymmetry list: use finite entries verbatim, else default 0.0..0.9
       scan_ng = count(is_finite(par%hgg_list))
       if (scan_ng > 0) then
          scan_glist(1:scan_ng) = pack(par%hgg_list, is_finite(par%hgg_list))
       else
          scan_ng = 10
          do i = 1, scan_ng
             scan_glist(i) = 0.1_wp * real(i-1, wp)      ! 0.0, 0.1, ..., 0.9
          enddo
       endif
    else
       !--- (a,g) axes collapse to the single simulated (albedo, g0).
       !--- This lets the tau scan run on its own as scatt(x,y,1,1,tau).
       scan_na = 1; scan_alist(1) = par%albedo
       scan_ng = 1; scan_glist(1) = par%hgg
    endif

    !--- polychromatic optical-depth (tau) list (target taumax values)
    if (par%use_tau_list) then
       scan_taumax_ref = tau_ref
       scan_nt = count(is_finite(par%tau_list))
       if (scan_nt > 0) then
          scan_tlist(1:scan_nt) = pack(par%tau_list, is_finite(par%tau_list))
       else
          !--- default: sqrt(2)-spaced tau/taumax in [0.5, 2] about the reference
          scan_nt = 5
          scan_tlist(1) = 0.5_wp       * scan_taumax_ref
          scan_tlist(2) = sqrt(0.5_wp) * scan_taumax_ref
          scan_tlist(3) = 1.0_wp       * scan_taumax_ref
          scan_tlist(4) = sqrt(2.0_wp) * scan_taumax_ref
          scan_tlist(5) = 2.0_wp       * scan_taumax_ref
       endif
       !--- auto-insert the reference (s = 1) slice if absent (bit-identical to a
       !--- single run; physically the simulated medium).
       have_ref = .false.
       do i = 1, scan_nt
          if (abs(scan_tlist(i) - scan_taumax_ref) <= stol*max(scan_taumax_ref, 1.0_wp)) have_ref = .true.
       enddo
       if (.not. have_ref .and. scan_nt < MAX_SCAN) then
          scan_nt = scan_nt + 1
          scan_tlist(scan_nt) = scan_taumax_ref
       endif
       if (scan_taumax_ref > 0.0_wp) then
          scan_s(1:scan_nt) = scan_tlist(1:scan_nt) / scan_taumax_ref
       else
          scan_s(1:scan_nt) = 1.0_wp
       endif
    else
       scan_nt         = 1
       scan_taumax_ref = tau_ref
       scan_tlist(1)   = tau_ref
       scan_s(1)       = 1.0_wp
    endif
  end subroutine scan_setup

  !---------------------------------------------------------------------------
  subroutine scan_reset_photon()
    implicit none
    scan_Aalb(1:scan_na) = 1.0_wp
    scan_Wg(1:scan_ng)   = 1.0_wp
    scan_Wtau(1:scan_nt) = 1.0_wp
  end subroutine scan_reset_photon

  !---------------------------------------------------------------------------
  !--- Advance the weight for each tau by one free-path segment of optical depth tau
  !--- (drawn at the reference): Wtau(it) <- s_it * exp((1-s_it)*tau) * Wtau(it).
  !--- Jonsson 2006, eq. 31; the forced first segment uses the same factor because
  !--- the (1-e^{-tau^e}) normalization cancels (see POLYCHROMATIC_TAU_SCAN_PLAN.md).
  subroutine scan_step_tau(tau)
    implicit none
    real(kind=wp), intent(in) :: tau
    integer :: it
    do it = 1, scan_nt
       scan_Wtau(it) = scan_Wtau(it) * scan_s(it) * exp((1.0_wp - scan_s(it)) * tau)
    enddo
  end subroutine scan_step_tau
  !---------------------------------------------------------------------------
end module scan_mod
