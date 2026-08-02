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
!--- The emission is computed from the populations above.  The extinction is
!--- the size-integrated curve precomputed for the same model, which the
!--- builder loads from SEDust/data/kext_*.dat and dust_extinction interpolates
!--- onto the model's own wavelength grid; the size integral behind that file
!--- is SEDust's size_integrated_extinction, run by its calc_kext.x.  Both
!--- halves therefore refer to the same grains -- the energy a cell absorbs is
!--- set by C_abs of the size distribution that then radiates it away -- for as
!--- long as the table matches the model actually built.  It does for the
!--- shipped tables and the default par%sed_* inputs; change the size
!--- distribution (par%sed_sizedist, par%sed_dl07_sdindex) and the table has to
!--- be regenerated with calc_kext.x, or the two halves part company.  Naming
!--- the model (par%dust_model) is otherwise enough to fix the dust physics.
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
  logical :: kext_named

  if (grain_model_ready) return
  st = 0

  !--- A blank par%sed_kext leaves the extinction table to SEDust's own default
  !--- for the model, so the standard table names stay in one place (SEDust).
  !--- A named file is instead mandatory: the build fails if it cannot be read.
  kext_named = len_trim(par%sed_kext) > 0

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
     if (kext_named) then
        call build_astrodust(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
                             par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st, &
                             kext_path=trim(par%sed_kext))
     else
        call build_astrodust(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
                             par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st)
     endif
  case ('dl07')
     if (kext_named) then
        call build_dl07(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
                        par%sed_dl07_sdindex, par%sed_dl07_uisrf, &
                        par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st, &
                        kext_path=trim(par%sed_kext))
     else
        call build_dl07(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
                        par%sed_dl07_sdindex, par%sed_dl07_uisrf, &
                        par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st)
     endif
  case ('zubko')
     if (kext_named) then
        call build_zubko(dmodel, trim(par%sed_zubko_config), trim(par%sed_zubko_dir), &
                         par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st, &
                         kext_path=trim(par%sed_kext))
     else
        call build_zubko(dmodel, trim(par%sed_zubko_config), trim(par%sed_zubko_dir), &
                         par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st)
     endif
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
     !--- A named extinction table is required, and is the last thing a builder
     !--- reads, so it is the likeliest source of a status from a build whose
     !--- optics paths are otherwise the defaults.
     if (mpar%p_rank == 0 .and. kext_named) write(*,'(3a)') &
        '       par%sed_kext = ''', trim(par%sed_kext), &
        ''' must be readable from par%sed_workdir.'
     call MPI_FINALIZE(ierr);  stop
  endif

  grain_model_ready = .true.
  end subroutine build_grain_model

  !---------------------------------------------------------------
  !--- Size-integrated extinction of the model, on the model's own wavelength
  !--- grid, in the (lambda, albedo, <cos>, C_ext) form the transport reads.
  !--- dust_extinction serves the precomputed curve of the model's kext table,
  !--- interpolated onto that grid -- a power law in lambda for the cross
  !--- sections, linear in ln(lambda) for <cos>.  The size integral behind the
  !--- table is
  !---
  !---   C_ext(lambda) = sum_pop sum_a dn(a) [C_abs(lambda,a) + C_sca(lambda,a)]
  !---   albedo        = C_sca / C_ext
  !---   <cos>         = sum dn C_sca g / sum dn C_sca
  !---
  !--- Cross sections are per H atom [cm^2/H], the same normalization as a
  !--- table read through par%kext_file.  Populations without scattering optics
  !--- (the PAHs) enter through absorption only, so the albedo falls where they
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
     if (mpar%p_rank == 0) then
        write(*,'(a,i0,a)') 'ERROR: dust_extinction failed (status=', st, ').'
        !--- The two statuses a run can actually provoke: no table was loaded
        !--- (the model's default file is missing under SEDust/data), or the
        !--- model grid reaches past the table, which is never extrapolated.
        if (st == 2) write(*,'(3a)') &
           '       no extinction table was loaded for the ''', trim(par%dust_model), &
           ''' model; set par%sed_kext, or restore SEDust/data/kext_*.dat.'
        if (st == 3) write(*,'(a)') &
           '       the model wavelength grid reaches outside that table.'
     endif
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
