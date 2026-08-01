module grain_model_mod
!---------------------------------------------------------------
!--- The grain population model of the run (SEDust dust_model_t): shared
!--- wavelength/size/temperature grids plus one population for each grain
!--- species and charge state, each carrying its size distribution dn(a),
!--- its optics C_abs/C_sca/<cos>(lambda,a), and its calorimetry.
!---
!--- The same object supplies both halves of the dust physics:
!---   emission   dust_emission()   -- dustemis_mod (Lucy 1999 iteration)
!---   extinction dust_extinction() -- sed_mod (the transport cross sections)
!--- Taking C_ext, albedo and <cos> from here rather than from a separate
!--- opacity file is what makes the absorbed power and the reemitted spectrum
!--- refer to the same grains: the energy a cell absorbs is set by C_abs of
!--- the same size distribution that then radiates it away.  Naming the model
!--- (par%dust_model) is therefore enough to fix the dust physics of a run.
!---
!--- Built once for each run; both callers go through build_grain_model, which
!--- returns immediately after the first call.
!---------------------------------------------------------------
  use define
  use dust_lib, only : dust_model_t, build_astrodust, build_dl07, build_zubko, &
                       dust_extinction, dust_nlam, dust_lambda
  implicit none
  private

  type(dust_model_t), save, public :: dmodel
  logical,            save, public :: grain_model_ready = .false.

  public :: build_grain_model, grain_extinction_table

contains
  !---------------------------------------------------------------
  !--- Build the model named by par%dust_model.  Every rank builds its own
  !--- copy from the same files, so the result is identical without a broadcast.
  subroutine build_grain_model()
  use mpi
  use ifport, only : chdir, getcwd
  implicit none
  integer :: ierr, cstat, st
  character(len=512) :: cwd_save

  if (grain_model_ready) return
  st = 0

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

  select case (trim(par%dust_model))
  case ('astrodust')
     call build_astrodust(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
                          par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st)
  case ('dl07')
     call build_dl07(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
                     par%sed_dl07_sdindex, par%sed_dl07_uisrf, &
                     par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st)
  case ('zubko')
     call build_zubko(dmodel, trim(par%sed_zubko_config), trim(par%sed_zubko_dir), &
                      par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st)
  case default
     cstat = chdir(trim(cwd_save))
     if (mpar%p_rank == 0) write(*,'(3a)') &
        'ERROR: par%dust_model = ''', trim(par%dust_model), &
        ''' unknown (use ''astrodust'', ''dl07'', or ''zubko'').'
     call MPI_FINALIZE(ierr);  stop
  end select

  cstat = chdir(trim(cwd_save))

  !--- SEDust reports a missing or malformed input through status instead of
  !--- stopping itself, so a bad par%sed_* path fails cleanly on every rank
  !--- here rather than aborting mid-build.
  if (st /= 0) then
     if (mpar%p_rank == 0) write(*,'(3a,i0,a)') &
        'ERROR: SEDust failed to build the ''', trim(par%dust_model), &
        ''' dust model (status=', st, '); check the par%sed_* input paths.'
     call MPI_FINALIZE(ierr);  stop
  endif

  grain_model_ready = .true.
  end subroutine build_grain_model

  !---------------------------------------------------------------
  !--- Size-integrated extinction of the model, on the model's own wavelength
  !--- grid, in the (lambda, albedo, <cos>, C_ext) form the transport reads:
  !---
  !---   C_ext(lambda) = sum_pop sum_a dn(a) [C_abs(lambda,a) + C_sca(lambda,a)]
  !---   albedo        = C_sca / C_ext
  !---   <cos>         = sum dn C_sca g / sum dn C_sca
  !---
  !--- Cross sections are per H atom [cm^2/H], the same normalization as the
  !--- kext table read from a file.  Populations without scattering optics (the
  !--- PAHs) enter through absorption only, so the albedo falls where they
  !--- dominate.  albedo and <cos> are zero wherever nothing scatters.
  subroutine grain_extinction_table(lam, alb, cosg, cext, n)
  use mpi
  implicit none
  real(kind=wp), allocatable, intent(out) :: lam(:), alb(:), cosg(:), cext(:)
  integer,                    intent(out) :: n
  real(kind=wp), allocatable :: cabs(:), csca(:)
  integer :: ierr, st, il

  call build_grain_model()

  n = dust_nlam(dmodel)
  if (n < 2) then
     if (mpar%p_rank == 0) write(*,'(a,i0,a)') &
        'ERROR: the dust model has ', n, ' wavelength points; at least 2 are needed.'
     call MPI_FINALIZE(ierr);  stop
  endif

  allocate(lam(n), alb(n), cosg(n), cext(n), cabs(n), csca(n))
  lam = dust_lambda(dmodel)
  call dust_extinction(dmodel, cext, cabs, csca, gbar=cosg, status=st)
  if (st /= 0) then
     if (mpar%p_rank == 0) write(*,'(a,i0,a)') &
        'ERROR: dust_extinction failed (status=', st, ').'
     call MPI_FINALIZE(ierr);  stop
  endif

  !--- Every model SEDust builds carries scattering optics, so a C_sca that
  !--- vanishes at every wavelength means the optics failed to load rather than
  !--- that the dust does not scatter.  An albedo of zero would turn grains
  !--- whose true albedo is ~0.5-0.7 in the optical into a purely absorbing
  !--- medium, so the run stops rather than transport that silently.
  if (maxval(csca) <= 0.0_wp) then
     if (mpar%p_rank == 0) write(*,'(3a)') &
        'ERROR: the ''', trim(par%dust_model), &
        ''' model returned no scattering at any wavelength;'
     if (mpar%p_rank == 0) write(*,'(a)') &
        '       check the par%sed_* optics paths.'
     call MPI_FINALIZE(ierr);  stop
  endif

  do il = 1, n
     if (cext(il) > 0.0_wp) then
        alb(il) = csca(il)/cext(il)
     else
        alb(il) = 0.0_wp
     endif
  enddo
  deallocate(cabs, csca)
  end subroutine grain_extinction_table

end module grain_model_mod
