module dustemis_mod
!--- Dust thermal emission via the SEDust library (MoCafe v2.00, Stage 3).
!--- Mode 1 (Lucy 1999): from the per-cell mean intensity J_lambda tallied by
!--- jtally_mod, compute each cell's dust emission spectrum with SEDust
!--- (equilibrium + stochastically heated grains + PAHs), then produce the
!--- emergent dust thermal SED.
!---
!--- This first version is NON-ITERATIVE: a single stellar pass sets J_lambda,
!--- and the dust emission is computed once.  It is exact in the limit that the
!--- dust is optically thin to its own (re-emitted) radiation.  Dust
!--- self-absorption (the Lucy iteration proper) is a later refinement.
!---
!--- Units.  MoCafe J_lambda is assembled in CGS [erg s^-1 cm^-2 sr^-1 um^-1]
!--- from the raw pathlength tally jt_sum = Sum(wgt*dl):
!---   J_cgs(il) = E_p * jt_sum(il,cell) / (4 pi V dlam) / dist_cm^2 ,
!--- with E_p = par%luminosity/nphotons [erg/s].  SEDust's Planck (bbody) is SI
!--- [W m^-2 sr^-1 m^-1], so J_SI = 1.0e3 * J_cgs.  dust_emission returns
!--- lamI_total(lambda) = lambda * (emission per H) whose SHAPE gives the
!--- spectral PDF; the absolute per-cell luminosity is set to the locally
!--- absorbed power (radiative equilibrium), which is exact by construction and
!--- sidesteps the SEDust absolute-normalization convention.
  use define
  use cellinfo_mod, only : ncell_total, cell_rhokap, cell_volume, cell_center, &
                           cell_random_position, cell_center_photon, car_ijk
  use dust_lib, only : dust_model_t, dust_emis_table_t, &
                       build_astrodust, build_dl07, build_zubko, &
                       dust_emission, dust_emission_single_teq, &
                       dust_build_table, dust_emission_interp, &
                       dust_nlam, dust_lambda
  implicit none
  private
  public :: setup_dustemis, compute_dustemis, write_dustemis, gen_dust_photon
  public :: Labs_total, dustemis_ready

  type(dust_model_t)         :: dmodel
  integer                    :: nl_sed = 0            ! SEDust wavelength grid length
  real(kind=wp), allocatable :: lam_sed(:)           ! SEDust lambda grid [um]
  !--- fast-path emission table (par%dust_fast_table): emission spectrum vs the
  !--- intensity scaling U for a fixed reference field shape J_ref.  Built once
  !--- from the stellar radiation field; per cell an O(us) log-log interpolation
  !--- replaces the O(0.1 s) exact SEDust solve.  Approximate when the local
  !--- spectral SHAPE departs from J_ref (e.g. hardened fields near a star).
  type(dust_emis_table_t)    :: emtab
  logical                    :: table_built = .false.
  real(kind=wp)              :: Jref_bol = 0.0_wp     ! bolometric integral of J_ref (SI) for U scaling
  real(kind=wp), allocatable :: cell_Lemit(:)        ! per-cell emitted (=absorbed) L [erg/s]
  real(kind=wp), allocatable :: cell_pdf(:,:)        ! (nl_mc, ncell) emission energy fraction per bin
  real(kind=wp), allocatable :: cell_cdfw(:,:)       ! (nl_mc, ncell) cumulative of cell_pdf (wavelength sampling)
  real(kind=wp), allocatable :: cell_Teq(:)          ! per-cell luminosity-weighted equilibrium-ish T [K] (diagnostic)
  real(kind=wp), allocatable :: cell_Lcdf(:)         ! (ncell) cumulative of cell_Lemit (cell sampling), normalized
  real(kind=wp) :: Labs_total = 0.0_wp               ! total absorbed L over the grid [erg/s]
  logical       :: dustemis_ready = .false.          ! .true. once a nonzero emission field exists
  integer :: ncell_tot = 0

contains
  !---------------------------------------------------------------
  subroutine setup_dustemis(grid)
  use mpi
  use ifport, only : chdir, getcwd
  implicit none
  type(grid_type), intent(in) :: grid
  integer :: ierr, cstat
  character(len=512) :: cwd_save

  !--- SEDust reads its dielectric tables via paths hard-coded relative to its
  !--- sed/ directory ('../data/dielectric/...'), so build the model from
  !--- par%sed_workdir (default SEDust/sed) and restore the working directory.
  !--- par%sed_qtable / par%sed_sizedist are given as absolute paths.
  !--- chdir/getcwd are the Intel IFPORT integer functions (return 0 on success).
  cstat = getcwd(cwd_save)
  if (len_trim(par%sed_workdir) > 0) then
     cstat = chdir(trim(par%sed_workdir))
     if (cstat /= 0 .and. mpar%p_rank == 0) write(*,'(3a)') &
        'WARNING: could not chdir to par%sed_workdir = ''', trim(par%sed_workdir), &
        ''' (SEDust dielectric files may not be found).'
  endif

  select case (trim(par%dust_model_sed))
  case ('astrodust')
     call build_astrodust(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
                          par%sed_NT, par%sed_Tlo, par%sed_Thi)
  case default
     cstat = chdir(trim(cwd_save))
     if (mpar%p_rank == 0) write(*,'(3a)') &
        'ERROR: par%dust_model_sed = ''', trim(par%dust_model_sed), &
        ''' not wired yet (use ''astrodust'').'
     call MPI_FINALIZE(ierr);  stop
  end select

  cstat = chdir(trim(cwd_save))

  nl_sed = dust_nlam(dmodel)
  lam_sed = dust_lambda(dmodel)          ! allocatable assignment
  ncell_tot = ncell_total(grid)

  if (mpar%p_rank == 0) then
     write(*,'(a)')      '--- Dust thermal emission (SEDust, Mode 1 Lucy, non-iterative) ---'
     write(*,'(2a)')     'dust model                : ', trim(dmodel%name)
     write(*,'(a,i6)')   'SEDust wavelength grid    : ', nl_sed
     write(*,'(a,2es11.3)') 'SEDust lambda range (um)  : ', lam_sed(1), lam_sed(nl_sed)
  endif
  end subroutine setup_dustemis

  !---------------------------------------------------------------
  !--- compute per-cell emission from the reduced J tally.
  subroutine compute_dustemis(grid)
  use mpi
  use sed_mod,    only : sed_nlam, sed_wave, sed_dwave, sed_sext, sed_albedo
  use jtally_mod, only : jt_sum
  implicit none
  type(grid_type), intent(in) :: grid
  real(kind=wp), allocatable :: Jcgs(:), Jsi(:), Jsed(:), lamI_sed(:), lamI_mc(:), emis(:)
  real(kind=wp) :: vol, dist2, jnorm, Labs, esum, lam_mean, csum, cell_Teq_true, rhk
  integer :: nl_mc, ic, il, ierr, ndone, nmine

  nl_mc = sed_nlam
  ndone = 0
  nmine = (ncell_tot - mpar%p_rank + mpar%nproc - 1)/mpar%nproc
  dist2 = par%distance2cm**2

  if (.not. allocated(cell_Lemit)) allocate(cell_Lemit(ncell_tot))
  if (.not. allocated(cell_pdf))   allocate(cell_pdf(nl_mc, ncell_tot))
  if (.not. allocated(cell_cdfw))  allocate(cell_cdfw(nl_mc, ncell_tot))
  if (.not. allocated(cell_Teq))   allocate(cell_Teq(ncell_tot))
  if (.not. allocated(cell_Lcdf))  allocate(cell_Lcdf(ncell_tot))
  cell_Lemit(:)  = 0.0_wp
  cell_pdf(:,:)  = 0.0_wp
  cell_Teq(:)    = 0.0_wp

  allocate(Jcgs(nl_mc), Jsi(nl_mc), Jsed(nl_sed), lamI_sed(nl_sed), lamI_mc(nl_mc), emis(nl_mc))

  !--- fast path: build the emission table once from the current (stellar)
  !--- field, then interpolate per cell instead of solving exactly.
  if (par%dust_fast_table .and. .not. table_built) call build_dust_table(grid)

  if (mpar%p_rank == 0) write(*,'(a,i0,a,a)') 'computing dust emission for ~', nmine, ' cells/rank ', &
     merge('(table interp)', '(exact solve) ', par%dust_fast_table)

  !--- distribute cells across MPI ranks (round-robin on the linear index).
  do ic = mpar%p_rank+1, ncell_tot, mpar%nproc
     ndone = ndone + 1
     if (mpar%p_rank == 0 .and. nmine >= 500 .and. mod(ndone, nmine/5) == 0) &
        write(*,'(a,i0,a,i0,a)') '   dust emission: rank0 ', ndone, '/', nmine, ' cells'
     rhk = cell_rhokap(grid, ic)
     if (rhk <= 0.0_wp) cycle
     if (all(jt_sum(:,ic) <= 0.0_wp)) cycle
     vol = cell_volume(grid, ic)

     !--- physical mean intensity in this cell (SI units for SEDust).
     !--- jt_sum already carries Lpacket [erg/s], so no E_p factor.
     jnorm = 1.0_wp/(fourpi*vol*dist2)
     do il = 1, nl_mc
        Jcgs(il) = jt_sum(il,ic)*jnorm/sed_dwave(il)        ! erg/s/cm^2/sr/um
     enddo
     Jsi(:) = Jcgs(:)*1.0e3_wp                              ! -> W/m^2/sr/m

     !--- absorbed power in this cell [erg/s] (radiative equilibrium => emitted).
     Labs = rhk * sum(jt_sum(:,ic)*sed_sext(:)*(1.0_wp - sed_albedo(:)))
     if (Labs <= 0.0_wp) cycle

     !--- SEDust emission spectrum (shape): resample J onto the SEDust grid,
     !--- then either interpolate the table at this cell's intensity scaling U
     !--- (fast) or solve exactly.  Resample lamI back onto the MoCafe grid.
     call resample_to_sed(sed_wave, Jsi, lam_sed, Jsed)
     if (par%dust_fast_table .and. table_built) then
        block
          real(kind=wp) :: Ucell
          Ucell = sum(Jsed(:))/max(Jref_bol, tinest)
          call dust_emission_interp(emtab, max(Ucell, tinest), lamI_sed)
        end block
     else if (par%dust_single_teq) then
        !--- fast equilibrium-only path (SEDust single mixture Teq): more
        !--- accurate than the B&W mixture-mean and gives the true Teq, but no
        !--- stochastic-heating / PAH features.
        block
          real(kind=wp) :: Teq
          call dust_emission_single_teq(dmodel, Jsed, lamI_sed, Teq)
          cell_Teq_true = Teq
        end block
     else
        call dust_emission(dmodel, Jsed, lamI_sed)
        cell_Teq_true = -1.0_wp
     endif
     call resample_from_sed(lam_sed, lamI_sed, sed_wave, lamI_mc)

     !--- emission energy fraction per MoCafe bin: j_lam ~ lamI/lam, energy ~ j_lam*dlam.
     do il = 1, nl_mc
        emis(il) = max(lamI_mc(il), 0.0_wp)/sed_wave(il) * sed_dwave(il)
     enddo
     esum = sum(emis)
     if (esum <= 0.0_wp) cycle

     cell_Lemit(ic)  = Labs
     cell_pdf(:,ic)  = emis(:)/esum
     if (cell_Teq_true > 0.0_wp) then
        !--- true equilibrium temperature from the single-Teq solver.
        cell_Teq(ic) = cell_Teq_true
     else
        !--- diagnostic "colour temperature": energy-weighted mean wavelength -> Wien.
        lam_mean     = sum(sed_wave(:)*cell_pdf(:,ic))
        cell_Teq(ic) = 2897.77_wp/max(lam_mean, 1.0e-30_wp)   ! Wien [K], lam in um
     endif
  enddo

  !--- combine partial per-cell results across ranks.
  call MPI_ALLREDUCE(MPI_IN_PLACE, cell_Lemit, ncell_tot,       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, cell_pdf,   nl_mc*ncell_tot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, cell_Teq,   ncell_tot,       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  Labs_total = sum(cell_Lemit)

  !--- build sampling tables for dust-emission photons (Lucy iteration):
  !--- cell selection CDF (prob. proportional to cell_Lemit) and per-cell
  !--- wavelength CDF (from cell_pdf).  Identical on every rank (inputs are
  !--- ALLREDUCEd), so all ranks sample from the same distributions.
  csum = 0.0_wp
  do ic = 1, ncell_tot
     csum         = csum + cell_Lemit(ic)
     cell_Lcdf(ic) = csum
  enddo
  if (csum > 0.0_wp) then
     cell_Lcdf(:) = cell_Lcdf(:)/csum
     dustemis_ready = .true.
  else
     dustemis_ready = .false.
  endif
  do ic = 1, ncell_tot
     if (cell_Lemit(ic) > 0.0_wp) then
        esum = 0.0_wp
        do il = 1, nl_mc
           esum = esum + cell_pdf(il,ic)
           cell_cdfw(il,ic) = esum
        enddo
     endif
  enddo

  if (mpar%p_rank == 0) then
     !--- Labs_total = SUM over cells of absorbed=emitted power.  With dust
     !--- self-absorption it counts each re-absorption, so it can exceed the
     !--- stellar luminosity; the EMERGENT dust luminosity (write_dustemis)
     !--- stays <= stellar (energy conservation).
     write(*,'(a,es14.6)') 'total cell emission (absorbed) L : ', Labs_total
     write(*,'(a,es14.6)') 'input stellar luminosity         : ', par%luminosity
     write(*,'(a,f10.5)')  'L_emit/L_star (>1 w/ self-abs)   : ', Labs_total/par%luminosity
  endif

  deallocate(Jcgs, Jsi, Jsed, lamI_sed, lamI_mc, emis)
  end subroutine compute_dustemis

  !---------------------------------------------------------------
  !--- build the emission table (fast path): reference field shape J_ref =
  !--- the mean physical J over illuminated cells, and a log-spaced intensity
  !--- grid U bracketing the per-cell heating range.  All ranks build the same
  !--- table (identical inputs) so per-cell interpolation is purely local.
  subroutine build_dust_table(grid)
  use mpi
  use sed_mod,    only : sed_nlam, sed_wave, sed_dwave
  use jtally_mod, only : jt_sum
  implicit none
  type(grid_type), intent(in) :: grid
  real(kind=wp), allocatable :: Jsi(:), Jsed(:), Jref_sed(:), Jsum(:), Ugrid(:)
  real(kind=wp) :: jnorm, dist2, Ucell, Umin, Umax, dlnU
  integer :: nl_mc, ic, il, ierr, ncnt, iu

  nl_mc = sed_nlam
  dist2 = par%distance2cm**2
  allocate(Jsi(nl_mc), Jsed(nl_sed), Jref_sed(nl_sed), Jsum(nl_sed))

  !--- pass 1: accumulate the mean physical field shape on the SEDust grid.
  Jsum(:) = 0.0_wp;  ncnt = 0
  do ic = mpar%p_rank+1, ncell_tot, mpar%nproc
     if (cell_rhokap(grid,ic) <= 0.0_wp) cycle
     if (all(jt_sum(:,ic) <= 0.0_wp)) cycle
     jnorm = 1.0e3_wp/(fourpi*cell_volume(grid,ic)*dist2)   ! -> SI (incl. cgs->SI 1e3)
     do il = 1, nl_mc
        Jsi(il) = jt_sum(il,ic)*jnorm/sed_dwave(il)
     enddo
     call resample_to_sed(sed_wave, Jsi, lam_sed, Jsed)
     Jsum(:) = Jsum(:) + Jsed(:)
     ncnt = ncnt + 1
  enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE, Jsum, nl_sed, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, ncnt, 1,      MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
  if (ncnt == 0) then
     deallocate(Jsi, Jsed, Jref_sed, Jsum);  return
  endif
  Jref_sed(:) = Jsum(:)/dble(ncnt)
  Jref_bol    = sum(Jref_sed)

  !--- pass 2: per-cell intensity scaling U = bolometric(J_cell)/bolometric(J_ref).
  Umin = hugest;  Umax = tinest
  do ic = mpar%p_rank+1, ncell_tot, mpar%nproc
     if (cell_rhokap(grid,ic) <= 0.0_wp) cycle
     if (all(jt_sum(:,ic) <= 0.0_wp)) cycle
     jnorm = 1.0e3_wp/(fourpi*cell_volume(grid,ic)*dist2)
     do il = 1, nl_mc
        Jsi(il) = jt_sum(il,ic)*jnorm/sed_dwave(il)
     enddo
     call resample_to_sed(sed_wave, Jsi, lam_sed, Jsed)
     Ucell = sum(Jsed)/max(Jref_bol, tinest)
     if (Ucell > 0.0_wp) then
        Umin = min(Umin, Ucell);  Umax = max(Umax, Ucell)
     endif
  enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE, Umin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, Umax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  if (.not. (Umax > Umin)) then
     Umin = 0.5_wp;  Umax = 2.0_wp
  endif

  !--- log-spaced U grid, padded 1.5x each side so all cells interpolate.
  allocate(Ugrid(par%dust_nU))
  Umin = Umin/1.5_wp;  Umax = Umax*1.5_wp
  dlnU = log(Umax/Umin)/dble(par%dust_nU-1)
  do iu = 1, par%dust_nU
     Ugrid(iu) = Umin*exp(dble(iu-1)*dlnU)
  enddo

  call dust_build_table(dmodel, Jref_sed, Ugrid, emtab)
  table_built = .true.
  if (mpar%p_rank == 0) then
     write(*,'(a,i0,a,2es11.3)') 'dust emission table built: NU=', par%dust_nU, '  U range ', Umin, Umax
  endif
  deallocate(Jsi, Jsed, Jref_sed, Jsum, Ugrid)
  end subroutine build_dust_table

  !---------------------------------------------------------------
  !--- generate a dust-emission photon (Lucy iteration): pick an emitting cell
  !--- with probability proportional to its luminosity, a uniform position in
  !--- the cell, an isotropic direction, and a wavelength from the cell's
  !--- emission spectrum.  Sets photon%Lpacket so the pass injects L_dust_total.
  subroutine gen_dust_photon(grid, photon, Lpacket)
  use random,  only : rand_number
  use sed_mod, only : sed_wave, sed_sext, sed_albedo, sed_hgg
  implicit none
  type(grid_type),   intent(in)  :: grid
  type(photon_type), intent(out) :: photon
  real(kind=wp),     intent(in)  :: Lpacket
  real(kind=wp) :: u, cost, sint, phi
  integer :: ic, lo, hi, mid, il

  !--- select cell by binary search on the (normalized) cumulative luminosity.
  u = rand_number()
  lo = 1;  hi = ncell_tot
  do while (lo < hi)
     mid = (lo+hi)/2
     if (u <= cell_Lcdf(mid)) then
        hi = mid
     else
        lo = mid+1
     endif
  enddo
  ic = lo

  !--- uniform position within the cell (sets photon cell/leaf index).
  call cell_random_position(grid, ic, photon)

  !--- wavelength from the cell emission CDF (binary search).
  u = rand_number()
  lo = 1;  hi = size(cell_cdfw,1)
  do while (lo < hi)
     mid = (lo+hi)/2
     if (u <= cell_cdfw(mid,ic)) then
        hi = mid
     else
        lo = mid+1
     endif
  enddo
  il = lo
  photon%il     = il
  photon%lambda = sed_wave(il)
  photon%s_ext  = sed_sext(il)
  photon%albedo = sed_albedo(il)
  photon%hgg    = sed_hgg(il)

  !--- isotropic emission direction.
  cost = 2.0_wp*rand_number() - 1.0_wp
  sint = sqrt(1.0_wp - cost*cost)
  phi  = twopi*rand_number()
  photon%kx = sint*cos(phi)
  photon%ky = sint*sin(phi)
  photon%kz = cost

  photon%nscatt  = 0
  photon%inside  = .true.
  photon%wgt     = 1.0_wp
  photon%Lpacket = Lpacket
  end subroutine gen_dust_photon

  !---------------------------------------------------------------
  !--- emergent dust thermal SED: raytrace each emitting cell to each
  !--- observer, attenuate at each wavelength, accumulate a per-observer
  !--- 1-D emergent SED and a 3-D (x,y,lambda) image.  Also writes the
  !--- intrinsic (unattenuated) SED and the per-cell T/L maps.
  subroutine write_dustemis(grid)
  use mpi
  use sed_mod,    only : sed_nlam, sed_wave, sed_dwave, sed_sext
  use iofile_mod
  use utility,    only : get_base_name
  implicit none
  type(grid_type), intent(in) :: grid
  type(io_file_type) :: file
  character(len=192) :: filename
  type(photon_type)  :: p
  real(kind=wp), allocatable :: sed_emergent(:,:), sed_intrinsic(:)
  real(kind=wp), allocatable :: dust_image(:,:,:)     ! (nxim,nyim,nl) surface brightness, observer 1
  real(kind=wp), allocatable :: Tcube(:,:,:), Lcube(:,:,:), Tleaf(:), Lleaf(:), leafxyz(:,:)
  real(kind=wp) :: tau, r2, atten, wcell, vdet(3), r, cx, cy, cz
  integer :: nl_mc, ic, i, j, k, il, kobs, ierr, status, ix, iy
  logical :: is_amr

  nl_mc = sed_nlam
  allocate(sed_emergent(nl_mc, par%nobs), sed_intrinsic(nl_mc))
  allocate(dust_image(par%nxim, par%nyim, nl_mc))
  sed_emergent(:,:)   = 0.0_wp
  sed_intrinsic(:)    = 0.0_wp
  dust_image(:,:,:)   = 0.0_wp

  !--- intrinsic (whole-grid) dust SED: sum of cell emission, no attenuation.
  do ic = 1, ncell_tot
     if (cell_Lemit(ic) > 0.0_wp) sed_intrinsic(:) = sed_intrinsic(:) + cell_Lemit(ic)*cell_pdf(:,ic)
  enddo

  !--- emergent SED and a pixel-resolved dust-emission image: attenuate each
  !--- cell's isotropic emission to each observer (deterministic, no MC noise;
  !--- IR scattering is negligible).  Distribute cells across ranks; raytrace
  !--- at the reference wavelength and scale tau per bin by s_ext.  The image
  !--- (observer 1) is a physical surface brightness [erg/s/cm^2/sr per bin].
  do ic = mpar%p_rank+1, ncell_tot, mpar%nproc
     if (cell_Lemit(ic) <= 0.0_wp) cycle
     call cell_center_photon(grid, ic, p)
     do kobs = 1, par%nobs
        p%kx = observer(kobs)%x - p%x
        p%ky = observer(kobs)%y - p%y
        p%kz = observer(kobs)%z - p%z
        r2   = p%kx**2 + p%ky**2 + p%kz**2
        r    = sqrt(r2)
        p%kx = p%kx/r;  p%ky = p%ky/r;  p%kz = p%kz/r
        call raytrace_to_edge(p, grid, tau)
        wcell = cell_Lemit(ic)/(fourpi*r2*par%distance2cm**2)
        do il = 1, nl_mc
           atten = exp(-sed_sext(il)*tau)
           sed_emergent(il,kobs) = sed_emergent(il,kobs) + wcell*cell_pdf(il,ic)*atten
        enddo
        !--- pixel-resolved image for observer 1 only.
        if (kobs == 1) then
           vdet(1) = observer(1)%rmatrix(1,1)*p%kx + observer(1)%rmatrix(1,2)*p%ky + observer(1)%rmatrix(1,3)*p%kz
           vdet(2) = observer(1)%rmatrix(2,1)*p%kx + observer(1)%rmatrix(2,2)*p%ky + observer(1)%rmatrix(2,3)*p%kz
           vdet(3) = observer(1)%rmatrix(3,1)*p%kx + observer(1)%rmatrix(3,2)*p%ky + observer(1)%rmatrix(3,3)*p%kz
           ix = floor(atan2(-vdet(1),vdet(3))*rad2deg/observer(1)%dxim+observer(1)%nxim/2.0_wp) + 1
           iy = floor(atan2(-vdet(2),vdet(3))*rad2deg/observer(1)%dyim+observer(1)%nyim/2.0_wp) + 1
           if (ix >= 1 .and. ix <= observer(1)%nxim .and. iy >= 1 .and. iy <= observer(1)%nyim) then
              do il = 1, nl_mc
                 dust_image(ix,iy,il) = dust_image(ix,iy,il) + &
                    wcell*cell_pdf(il,ic)*exp(-sed_sext(il)*tau)/observer(1)%steradian_pix
              enddo
           endif
        endif
     enddo
  enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE, sed_emergent, nl_mc*par%nobs, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  do il = 1, nl_mc
     call MPI_ALLREDUCE(MPI_IN_PLACE, dust_image(:,:,il), par%nxim*par%nyim, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  enddo

  if (mpar%p_rank /= 0) then
     deallocate(sed_emergent, sed_intrinsic, dust_image)
     return
  endif

  status = 0
  filename = trim(get_base_name(par%out_file))//'_dustsed'//trim(io_file_extension(par%file_format))
  call io_open_new(file, trim(filename), status)
  call io_append_image(file, sed_emergent, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','SED_emergent','emergent dust SED L_lam-like (nlam,nobs)',status)
  call io_put_keyword(file,'SED_UNIT','luminosity/(4pi dist_cm^2)/bin','emergent SED unit',status)
  call io_put_keyword(file,'TOT_LUM', par%luminosity, 'input stellar luminosity',status)
  call io_put_keyword(file,'L_ABS',   Labs_total,     'total absorbed = emitted L',status)
  call io_put_keyword(file,'DIST_CM', par%distance2cm,'distance unit (cm)',status)
  call io_append_image(file, sed_intrinsic, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','SED_intrinsic','intrinsic (unattenuated) dust SED [erg/s per bin]',status)
  call io_append_image(file, sed_wave, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','Wavelength','bin centers [um]',status)
  call io_append_image(file, sed_dwave, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','Dwavelength','bin widths [um]',status)
  !--- per-cell T and L maps: Cartesian cube (nx,ny,nz) or per-leaf arrays (AMR).
  is_amr = trim(par%grid_type) == 'amr'
  if (is_amr) then
     allocate(Tleaf(ncell_tot), Lleaf(ncell_tot), leafxyz(ncell_tot,3))
     do ic = 1, ncell_tot
        Tleaf(ic) = cell_Teq(ic);  Lleaf(ic) = cell_Lemit(ic)
        call cell_center(grid, ic, cx, cy, cz)
        leafxyz(ic,1) = cx;  leafxyz(ic,2) = cy;  leafxyz(ic,3) = cz
     enddo
     call io_append_image(file, Tleaf, status, bitpix=-64)
     call io_put_keyword(file,'EXTNAME','Tdust','per-leaf dust temperature [K]',status)
     call io_append_image(file, Lleaf, status, bitpix=-64)
     call io_put_keyword(file,'EXTNAME','Ldust','per-leaf emitted luminosity [erg/s]',status)
     call io_append_image(file, leafxyz, status, bitpix=-64)
     call io_put_keyword(file,'EXTNAME','LeafXYZ','leaf center x,y,z (code units)',status)
  else
     allocate(Tcube(grid%nx,grid%ny,grid%nz), Lcube(grid%nx,grid%ny,grid%nz))
     do ic = 1, ncell_tot
        call car_ijk(grid, ic, i, j, k)
        Tcube(i,j,k) = cell_Teq(ic);  Lcube(i,j,k) = cell_Lemit(ic)
     enddo
     call io_append_image(file, Tcube, status, bitpix=-64)
     call io_put_keyword(file,'EXTNAME','Tdust','per-cell dust temperature [K]',status)
     call io_append_image(file, Lcube, status, bitpix=-64)
     call io_put_keyword(file,'EXTNAME','Ldust','per-cell emitted luminosity [erg/s]',status)
  endif
  !--- pixel-resolved emergent dust-emission image (observer 1), surface
  !--- brightness [erg/s/cm^2/sr per bin]: dust emission in the observer image.
  call io_append_image(file, dust_image, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','DustEmis_image','emergent dust emission SB (x,y,lambda), observer 1',status)
  call io_put_keyword(file,'SB_UNIT','erg/s/cm^2/sr/bin','surface-brightness unit',status)
  call io_close(file, status)
  write(*,'(2a)') 'dust SED written to: ', trim(filename)

  !--- energy-conservation report: the emergent bolometric dust luminosity
  !--- (4 pi dist_cm^2 * integral of the emergent SED over the sky, estimated
  !--- from observer 1 as isotropic) should approach the stellar energy
  !--- absorbed in the FIRST generation.  Print the intrinsic and emergent
  !--- bolometric sums for a quick check.
  write(*,'(a,es14.6)') 'intrinsic dust L (sum cells)   : ', sum(sed_intrinsic)
  write(*,'(a,es14.6)') 'emergent dust L (4pi d^2 * SED): ', &
     fourpi*(par%distance*par%distance2cm)**2*sum(sed_emergent(:,1))

  deallocate(sed_emergent, sed_intrinsic, dust_image)
  if (allocated(Tcube)) deallocate(Tcube, Lcube)
  if (allocated(Tleaf)) deallocate(Tleaf, Lleaf, leafxyz)
  end subroutine write_dustemis

  !---------------------------------------------------------------
  !--- log-log resample of y(xin) onto xout (both ascending), clamped.
  subroutine resample_to_sed(xin, yin, xout, yout)
  implicit none
  real(kind=wp), intent(in)  :: xin(:), yin(:), xout(:)
  real(kind=wp), intent(out) :: yout(:)
  call loglog_interp(xin, yin, xout, yout)
  end subroutine resample_to_sed

  subroutine resample_from_sed(xin, yin, xout, yout)
  implicit none
  real(kind=wp), intent(in)  :: xin(:), yin(:), xout(:)
  real(kind=wp), intent(out) :: yout(:)
  call loglog_interp(xin, yin, xout, yout)
  end subroutine resample_from_sed

  subroutine loglog_interp(x, y, xnew, ynew)
  implicit none
  real(kind=wp), intent(in)  :: x(:), y(:), xnew(:)
  real(kind=wp), intent(out) :: ynew(:)
  integer :: n, m, i, iq
  real(kind=wp) :: lx, t, ly
  n = size(x);  m = size(xnew)
  do iq = 1, m
     lx = log(xnew(iq))
     if (xnew(iq) <= x(1)) then
        ynew(iq) = max(y(1), 0.0_wp)
     else if (xnew(iq) >= x(n)) then
        ynew(iq) = max(y(n), 0.0_wp)
     else
        do i = 2, n
           if (xnew(iq) <= x(i)) exit
        enddo
        t = (lx - log(x(i-1)))/(log(x(i)) - log(x(i-1)))
        if (y(i-1) > 0.0_wp .and. y(i) > 0.0_wp) then
           ly = log(y(i-1)) + t*(log(y(i)) - log(y(i-1)))
           ynew(iq) = exp(ly)
        else
           ynew(iq) = max((1.0_wp-t)*y(i-1) + t*y(i), 0.0_wp)
        endif
     endif
  enddo
  end subroutine loglog_interp

end module dustemis_mod
