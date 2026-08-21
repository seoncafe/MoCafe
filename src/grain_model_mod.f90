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
!--- builder loads from /kext of SEDust/data/<model>/sedust_<model>.h5 and
!--- dust_extinction interpolates onto the model's own wavelength grid; the
!--- size integral behind that curve is SEDust's size_integrated_extinction,
!--- run by its calc_kext.x.  Both halves therefore refer to the same grains --
!--- the energy a cell absorbs is set by C_abs of the size distribution that
!--- then radiates it away -- because the optics and the curve come out of the
!--- one file.  Change the size distribution (par%sed_dl07_sdindex, or a
!--- par%sed_kext naming a curve computed for another one) and the two halves
!--- part company.  Naming the model (par%dust_model) is otherwise enough to
!--- fix the dust physics.
!---
!--- Built once for each run; both callers go through build_grain_model, which
!--- returns immediately after the first call.
!---------------------------------------------------------------
  use define
  use dust_lib, only : dust_model_t, build_dust, &
                       dust_extinction, dust_nlam, dust_lambda
  implicit none
  private

  type(dust_model_t), save, public :: dmodel
  logical,            save, public :: grain_model_ready = .false.

  !--- Lyman limit [um].  Every SEDust product stores one wavelength axis with
  !--- this point marked in it, and the grain model is built on the axis above
  !--- it or on the whole of it.
  real(kind=wp), parameter :: LAM_LYMAN_UM = 0.0912_wp

  public :: build_grain_model, grain_extinction_table

contains
  !---------------------------------------------------------------
  !--- Build the model named by par%dust_model.  Every rank builds its own
  !--- copy from the same files, so the result is identical without a broadcast.
  subroutine build_grain_model()
  use mpi
  implicit none
  integer :: ierr, st
  character(len=256) :: why
  logical :: kext_named, euv

  if (grain_model_ready) return
  st = 0

  !--- A blank par%sed_kext leaves the extinction curve to /kext of the model's
  !--- own product, so the standard names stay in one place (SEDust).  A named
  !--- file is instead mandatory: the build fails if it cannot be read.
  kext_named = len_trim(par%sed_kext) > 0

  !--- Which view of the model's stored wavelength axis to build on.  The axis
  !--- is one array with the Lyman limit marked in it, and dust_extinction is
  !--- refused rather than extrapolated outside the grid the model was built
  !--- on, so the grid has to reach at least as far as the transport does --
  !--- and no further, the ionizing part being what makes every emission solve
  !--- roughly 2.4x more expensive.  par%lambda_min is that boundary.  The cut
  !--- is an index cut at the last node at or below the Lyman limit, the same
  !--- for every model, so the grid covers that limit whichever model is named.
  euv = par%lambda_min < LAM_LYMAN_UM

  !--- One entry point for every model SEDust codes: the model name selects the
  !--- builder, and par%sed_datadir is the data root for the whole build -- the
  !--- optics product, the size distribution, the dielectric functions and the
  !--- default extinction curve all resolve inside it -- so the run names one
  !--- directory and stays where it is.  sd_index/u_isrf are DL07's alone and
  !--- are ignored by the other models.
  if (kext_named) then
     call build_dust(dmodel, trim(par%dust_model), trim(par%sed_datadir), &
                     par%sed_NT, par%sed_Tlo, par%sed_Thi, include_euv=euv, &
                     status=st, message=why, sd_index=par%sed_dl07_sdindex, &
                     u_isrf=par%sed_dl07_uisrf, kext_path=trim(par%sed_kext))
  else
     call build_dust(dmodel, trim(par%dust_model), trim(par%sed_datadir), &
                     par%sed_NT, par%sed_Tlo, par%sed_Thi, include_euv=euv, &
                     status=st, message=why, sd_index=par%sed_dl07_sdindex, &
                     u_isrf=par%sed_dl07_uisrf)
  endif

  !--- SEDust reports a missing or malformed input through status instead of
  !--- stopping itself, so a bad par%sed_* path fails cleanly on every rank
  !--- here rather than aborting mid-build.  message names the stage in words,
  !--- in one vocabulary shared by every model.
  if (st /= 0) then
     if (mpar%p_rank == 0) then
        write(*,'(3a,i0,a)') 'ERROR: SEDust failed to build the ''', &
           trim(par%dust_model), ''' dust model (status=', st, ').'
        write(*,'(2a)') '       ', trim(why)
        if (st == 90) write(*,'(3a)') &
           '       par%dust_model = ''', trim(par%dust_model), &
           ''' is not one of ''astrodust'', ''dl07'', ''mrn'', ''zubko''.'
        if (st /= 90) write(*,'(3a)') &
           '       par%sed_datadir = ''', trim(par%sed_datadir), &
           ''' must hold <dust_model>/sedust_<dust_model>.h5.'
        !--- A named extinction curve is required, unlike the model's own, so
        !--- it is the likeliest source of a status from a build whose other
        !--- inputs are the defaults.
        if (kext_named) write(*,'(3a)') &
           '       par%sed_kext = ''', trim(par%sed_kext), ''' must be readable.'
     endif
     call MPI_FINALIZE(ierr);  stop
  endif

  grain_model_ready = .true.
  end subroutine build_grain_model

  !---------------------------------------------------------------
  !--- Size-integrated extinction of the model, on the model's own wavelength
  !--- grid, in the (lambda, albedo, <cos>, C_ext) form the transport reads.
  !--- dust_extinction serves the precomputed curve of the model's own product,
  !--- interpolated onto that grid -- a power law in lambda for the cross
  !--- sections, linear in ln(lambda) for <cos>.  The size integral behind the
  !--- curve is
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
        !--- The two statuses a run can actually provoke: no curve was loaded
        !--- (the model's product is missing under par%sed_datadir), or the
        !--- model grid reaches past the curve, which is never extrapolated.
        if (st == 2) write(*,'(3a)') &
           '       no extinction curve was loaded for the ''', trim(par%dust_model), &
           ''' model; set par%sed_kext, or restore its sedust_*.h5.'
        if (st == 3) write(*,'(a)') &
           '       the model wavelength grid reaches outside that curve.'
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
