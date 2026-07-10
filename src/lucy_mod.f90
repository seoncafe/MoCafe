module lucy_mod
!--- Lucy (1999) iteration driver for dust self-absorption (MoCafe v2.00,
!--- Stage 3 follow-up #1).  Energy-only passes (no observer peel) transport
!--- stellar photons once and dust-emission photons each iteration; the
!--- per-cell J_lambda tally feeds SEDust (dustemis_mod) to update the
!--- emission, iterating until the total emitted (=absorbed) luminosity
!--- converges.  The stellar tally is computed once and reused; each iteration
!--- adds a fresh dust-photon tally on top.
!---
!--- After convergence jt_sum holds the total (stellar + dust) J_lambda; the
!--- caller then runs the normal imaging pass (peel on) and writes the dust SED.
  use define
  use photon_mod,   only : gen_photon
  use dustemis_mod, only : gen_dust_photon, compute_dustemis, Labs_total, dustemis_ready
  use jtally_mod,   only : jt_on, jt_first, jt_sum, jt_eabs, jtally_reduce
  use peelingoff_mod, only : peel_enabled
  use mrw_mod,      only : mrw_on, mrw_try_step
  use random,       only : rand_number
  use mpi
  implicit none
  private
  public :: run_lucy_iteration

contains
  !---------------------------------------------------------------
  subroutine run_lucy_iteration(grid)
  implicit none
  type(grid_type), intent(inout) :: grid
  type(photon_type) :: photon
  real(kind=wp), allocatable :: jt_star(:,:,:,:)
  real(kind=wp) :: eabs_star, Lstar_packet, Ldust_packet, Lprev, drel
  integer(kind=int64) :: n_star, n_dust, ip
  integer :: iter, ierr

  peel_enabled = .false.          ! energy-only passes
  n_star = par%nphotons
  n_dust = int(par%dust_no_photons, int64)

  !--- stellar energy pass (once): tally J_star.
  jt_sum(:,:,:,:) = 0.0_wp;  jt_eabs = 0.0_wp
  do ip = mpar%p_rank+1, n_star, mpar%nproc
     call gen_photon(grid, photon)        ! Lpacket = luminosity/nphotons
     call transport(photon, grid)
  enddo
  call jtally_reduce()                     ! jt_sum, jt_eabs now full (all ranks)
  allocate(jt_star(size(jt_sum,1), size(jt_sum,2), size(jt_sum,3), size(jt_sum,4)))
  jt_star   = jt_sum
  eabs_star = jt_eabs

  call compute_dustemis(grid)              ! emission from stellar J
  Lprev = Labs_total
  if (mpar%p_rank == 0) write(*,'(a,i3,a,es13.6)') 'Lucy iter ', 1, ': L_emit = ', Labs_total

  !--- dust self-absorption iterations.
  do iter = 2, par%dust_niter
     if (.not. dustemis_ready) exit
     Ldust_packet = Labs_total/dble(n_dust)

     !--- dust-photon energy pass: tally the dust contribution alone.
     jt_sum(:,:,:,:) = 0.0_wp;  jt_eabs = 0.0_wp
     do ip = mpar%p_rank+1, n_dust, mpar%nproc
        call gen_dust_photon(grid, photon, Ldust_packet)
        call transport(photon, grid)
     enddo
     call jtally_reduce()                  ! jt_sum = full dust J, jt_eabs = full dust absorbed

     !--- total field = stellar + dust; recompute emission.
     jt_sum  = jt_sum + jt_star
     jt_eabs = jt_eabs + eabs_star
     call compute_dustemis(grid)

     drel = abs(Labs_total - Lprev)/max(Labs_total, tinest)
     if (mpar%p_rank == 0) write(*,'(a,i3,a,es13.6,a,es10.3)') &
        'Lucy iter ', iter, ': L_emit = ', Labs_total, '  rel.change = ', drel
     Lprev = Labs_total
     if (drel < par%dust_tol) then
        if (mpar%p_rank == 0) write(*,'(a,i0,a)') 'Lucy iteration converged after ', iter, ' iterations.'
        exit
     endif
  enddo

  deallocate(jt_star)
  peel_enabled = .true.           ! restore for the imaging pass
  !--- disable further tallying: jt_sum already holds the converged total J
  !--- (the imaging pass must not overwrite it).
  jt_on = .false.
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end subroutine run_lucy_iteration

  !---------------------------------------------------------------
  !--- energy-only photon transport: forced first scattering (with the
  !--- analytic first-flight J tally), then the standard scatter/free-path
  !--- loop.  Identical to the run_simulation body but peel routines are
  !--- no-ops (peel_enabled = .false.).  s_ext converts reference-wavelength
  !--- optical depths to the photon wavelength (see run_simulation_mod).
  subroutine transport(photon, grid)
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  real(kind=wp) :: tau, tau0, wgt1
  real(kind=wp), parameter :: RR_CUT = 1.0e-4_wp   ! Russian-roulette weight cutoff
  real(kind=wp), parameter :: RR_SURV = 0.1_wp     ! survival probability
  integer(kind=int64) :: nstep
  integer(kind=int64), parameter :: NSTEP_MAX = 100000000_int64

  nstep = 0
  do while (photon%inside)
     nstep = nstep + 1
     if (nstep > NSTEP_MAX) exit          ! safety cap against pathological walks
     !--- Modified Random Walk: in a very optically thick cell take one big
     !--- diffusion step instead of the normal free path (deposits into the J
     !--- tally and moves the photon to the inscribed-sphere surface).
     if (mrw_on .and. photon%nscatt > 0) then
        if (mrw_try_step(photon, grid)) then
           !--- Russian roulette after the (absorption-decayed) MRW weight.
           if (photon%wgt < RR_CUT) then
              if (rand_number() > RR_SURV) exit
              photon%wgt = photon%wgt/RR_SURV
           endif
           cycle
        endif
     endif
     if (photon%nscatt == 0) then
        jt_first = jt_on
        call raytrace_to_edge(photon, grid, tau0)
        jt_first = .false.
        tau0       = tau0 * photon%s_ext
        wgt1       = 1.0_wp - exp(-tau0)
        photon%wgt = photon%wgt * wgt1
        if (tau0 > 0.0_wp) then
           tau = -log(1.0_wp - rand_number()*wgt1)
        else
           tau = hugest
        endif
     else
        tau = -log(rand_number())
     endif
     call raytrace_to_tau(photon, grid, tau/photon%s_ext)
     if (photon%inside) then
        call scattering(photon, grid)
        !--- Russian roulette: unbiasedly terminate low-weight photons so
        !--- high-albedo, high-tau clouds do not diffuse forever.
        if (photon%wgt < RR_CUT) then
           if (rand_number() > RR_SURV) exit
           photon%wgt = photon%wgt/RR_SURV
        endif
     endif
  enddo
  end subroutine transport

end module lucy_mod
