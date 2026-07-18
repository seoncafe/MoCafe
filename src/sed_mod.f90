module sed_mod
!--- Multi-wavelength (SED) infrastructure for MoCafe v2.00 (Stage 1).
!--- Provides: a log-spaced wavelength grid; wavelength-dependent dust
!--- properties C_ext(lambda), albedo(lambda), g(lambda) read from an
!--- extinction table (e.g. SEDust calc_kext_astrodust.x output); and a
!--- stellar source spectrum (Planck or 2-column file) sampled per photon
!--- with an alias table.
!---
!--- Transport of a photon at wavelength bin il uses the grey-rescaling
!--- factor sed_sext(il) = C_ext(il)/C_ext(lambda_ref) (photon%s_ext):
!--- the grid rhokap is the opacity at the reference wavelength, and every
!--- optical depth is scaled by s_ext (same mathematics as the Jonsson 2006
!--- tau scan), so the raytrace routines (car/clump/amr) need no changes.
  use define
  use random,  only : random_alias_setup, rand_alias_choise
  implicit none
  public

  integer :: sed_nlam = 0
  real(kind=wp), allocatable :: sed_wave(:)     ! bin centers [um]
  real(kind=wp), allocatable :: sed_dwave(:)    ! bin widths  [um]
  real(kind=wp), allocatable :: sed_cext(:)     ! C_ext/H [cm^2/H] at bin centers
  real(kind=wp), allocatable :: sed_albedo(:)   ! dust albedo(lambda)
  real(kind=wp), allocatable :: sed_hgg(:)      ! asymmetry g(lambda)
  real(kind=wp), allocatable :: sed_sext(:)     ! C_ext(lambda)/C_ext(lambda_ref)
  real(kind=wp), allocatable :: sed_lum(:)      ! source luminosity fraction per bin (sum = 1)
  real(kind=wp), allocatable :: sed_src_pdf(:)  ! alias probability table
  integer,       allocatable :: sed_src_alias(:)
  real(kind=wp), allocatable :: sed_src_cdf(:)  ! cumulative luminosity (quasi-random inverse-CDF path)
  real(kind=wp) :: sed_cext_ref = 0.0_wp        ! C_ext/H at par%lambda_ref

  !--- external-field spectrum (SED mode with external illumination).
  logical :: sed_ext_on = .false.
  real(kind=wp), allocatable :: sed_ext_lum(:)   ! external spectrum fraction per bin (sum = 1)
  real(kind=wp), allocatable :: sed_ext_pdf(:)   ! alias probability table
  integer,       allocatable :: sed_ext_alias(:)
  real(kind=wp), allocatable :: sed_ext_cdf(:)   ! cumulative external intensity (inverse-CDF path)

  !--- unit-conversion constants for physical par%spectrum_type files.
  real(kind=wp), parameter :: c_um     = 2.99792458e14_wp  ! speed of light [um/s]
  real(kind=wp), parameter :: hc_evum  = 1.23984193_wp     ! h*c [eV um]

contains
  !---------------------------------------------------------------
  subroutine setup_sed()
  use mpi
  implicit none
  ! local variables
  real(kind=wp), allocatable :: tb_lam(:), tb_alb(:), tb_cos(:), tb_cext(:)
  real(kind=wp), allocatable :: sp_lam(:), sp_lum(:)
  real(kind=wp), allocatable :: edge(:)
  real(kind=wp) :: dlnlam, lum_sum
  integer       :: ntab, nsp, il, ierr
  logical       :: spec_is_absolute, spec_derived

  if (par%nlambda < 2) then
     if (mpar%p_rank == 0) write(*,'(a)') 'ERROR: par%nlambda must be >= 2 in SED mode.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (.not. (par%lambda_max > par%lambda_min .and. par%lambda_min > 0.0_wp)) then
     if (mpar%p_rank == 0) write(*,'(a)') 'ERROR: require 0 < lambda_min < lambda_max in SED mode.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (len_trim(par%kext_file) == 0) then
     if (mpar%p_rank == 0) write(*,'(a)') 'ERROR: SED mode requires par%kext_file (lambda, albedo, <cos>, C_ext/H table).'
     call MPI_FINALIZE(ierr);  stop
  endif

  !--- log-spaced wavelength grid (bin edges and geometric bin centers).
  sed_nlam = par%nlambda
  allocate(edge(sed_nlam+1))
  allocate(sed_wave(sed_nlam), sed_dwave(sed_nlam))
  allocate(sed_cext(sed_nlam), sed_albedo(sed_nlam), sed_hgg(sed_nlam), sed_sext(sed_nlam))
  allocate(sed_lum(sed_nlam), sed_src_pdf(sed_nlam), sed_src_alias(sed_nlam), sed_src_cdf(sed_nlam))
  dlnlam = log(par%lambda_max/par%lambda_min)/sed_nlam
  do il = 1, sed_nlam+1
     edge(il) = par%lambda_min * exp((il-1)*dlnlam)
  enddo
  do il = 1, sed_nlam
     sed_wave(il)  = sqrt(edge(il)*edge(il+1))
     sed_dwave(il) = edge(il+1) - edge(il)
  enddo

  !--- read the dust extinction table (rank 0) and broadcast.
  if (mpar%p_rank == 0) call read_kext_table(trim(par%kext_file), tb_lam, tb_alb, tb_cos, tb_cext, ntab)
  call MPI_BCAST(ntab, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if (mpar%p_rank /= 0) allocate(tb_lam(ntab), tb_alb(ntab), tb_cos(ntab), tb_cext(ntab))
  call MPI_BCAST(tb_lam,  ntab, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tb_alb,  ntab, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tb_cos,  ntab, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tb_cext, ntab, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  !--- interpolate onto the RT wavelength grid: C_ext log-log, albedo/g
  !--- linear in ln(lambda); clamp outside the table range.
  do il = 1, sed_nlam
     sed_cext(il)   = exp(interp_clamped(log(tb_lam), log(tb_cext), log(sed_wave(il))))
     sed_albedo(il) = interp_clamped(log(tb_lam), tb_alb, log(sed_wave(il)))
     sed_hgg(il)    = interp_clamped(log(tb_lam), tb_cos, log(sed_wave(il)))
  enddo
  sed_cext_ref  = exp(interp_clamped(log(tb_lam), log(tb_cext), log(par%lambda_ref)))
  sed_sext(:)   = sed_cext(:)/sed_cext_ref
  !--- the grid opacity (rhokap) is defined at the reference wavelength.
  par%cext_dust = sed_cext_ref
  par%albedo    = interp_clamped(log(tb_lam), tb_alb, log(par%lambda_ref))
  par%hgg       = interp_clamped(log(tb_lam), tb_cos, log(par%lambda_ref))
  par%lambda0   = par%lambda_ref

  !--- source spectrum.  With multiple source components (par%nsource > 1) the
  !--- spectra for each source are set up in sources_mod, so the global single-source
  !--- spectrum is optional here: build a flat placeholder and skip the checks.
  if (par%nsource > 1 .and. len_trim(par%source_spectrum) == 0 .and. par%tstar <= 0.0_wp) then
     sed_lum(:) = sed_dwave(:)
  else if (len_trim(par%source_spectrum) > 0) then
     if (mpar%p_rank == 0) then
        call read_spectrum_file(trim(par%source_spectrum), sp_lam, sp_lum, nsp)
        call convert_spectrum_units(sp_lam, sp_lum)
     endif
     call MPI_BCAST(nsp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     if (mpar%p_rank /= 0) allocate(sp_lam(nsp), sp_lum(nsp))
     call MPI_BCAST(sp_lam, nsp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(sp_lum, nsp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     do il = 1, sed_nlam
        !--- outside the tabulated spectrum: zero luminosity (no clamping).
        if (sed_wave(il) < sp_lam(1) .or. sed_wave(il) > sp_lam(nsp)) then
           sed_lum(il) = 0.0_wp
        else
           sed_lum(il) = interp_clamped(log(sp_lam), sp_lum, log(sed_wave(il))) * sed_dwave(il)
        endif
     enddo
     deallocate(sp_lam, sp_lum)
  else if (par%tstar > 0.0_wp) then
     do il = 1, sed_nlam
        sed_lum(il) = planck_shape(sed_wave(il), par%tstar) * sed_dwave(il)
     enddo
  else if (trim(par%source_geometry(1:8)) == 'external') then
     !--- external illumination: the internal (point/extended) source spectrum
     !--- is unused; the external-field spectrum is set in setup_sed_external.
     !--- Build a flat placeholder so the checks below stay valid.
     sed_lum(:) = sed_dwave(:)
  else
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: SED mode requires a source spectrum: par%source_spectrum (file) or par%tstar (Planck, K).'
     call MPI_FINALIZE(ierr);  stop
  endif

  lum_sum = sum(sed_lum)
  if (.not. (lum_sum > 0.0_wp)) then
     if (mpar%p_rank == 0) write(*,'(a)') 'ERROR: source spectrum has zero luminosity on the wavelength grid.'
     call MPI_FINALIZE(ierr);  stop
  endif

  !--- absolute (physical-type) file spectrum: lum_sum is the [erg/s] luminosity
  !--- integrated over the wavelength grid.  Derive par%luminosity from it when
  !--- the scale is unset (sentinel); a set scale rescales the file to it (the
  !--- normalization below is unchanged, so only par%luminosity carries the scale).
  spec_is_absolute = (trim(par%spectrum_type) /= 'shape') .and. len_trim(par%source_spectrum) > 0
  spec_derived     = .false.
  if (spec_is_absolute .and. par%luminosity < -900.0_wp) then
     par%luminosity = lum_sum
     spec_derived   = .true.
  endif

  sed_lum(:) = sed_lum(:)/lum_sum

  !--- alias table for luminosity-weighted wavelength-bin sampling.
  sed_src_pdf(:) = sed_lum(:)
  call random_alias_setup(sed_src_pdf, sed_src_alias)

  !--- monotone cumulative distribution for the quasi-random inverse-CDF
  !--- wavelength sampler (sample_sed_lambda_u); the alias table above stays the
  !--- default sampling path so its stream is unchanged.
  sed_src_cdf(1) = sed_lum(1)
  do il = 2, sed_nlam
     sed_src_cdf(il) = sed_src_cdf(il-1) + sed_lum(il)
  enddo
  sed_src_cdf(sed_nlam) = 1.0_wp

  if (mpar%p_rank == 0) then
     write(*,'(a)')          '--- SED (multi-wavelength) mode ---'
     write(*,'(a,i6)')       'N wavelength bins         : ', sed_nlam
     write(*,'(a,2es12.4)')  'lambda_min, lambda_max(um): ', par%lambda_min, par%lambda_max
     write(*,'(a,es12.4)')   'reference lambda (um)     : ', par%lambda_ref
     write(*,'(a,es12.4)')   'C_ext/H at lambda_ref     : ', sed_cext_ref
     write(*,'(a,2f8.4)')    'albedo, g at lambda_ref   : ', par%albedo, par%hgg
     if (par%tstar > 0.0_wp .and. len_trim(par%source_spectrum) == 0) then
        write(*,'(a,es12.4)') 'source: Planck T_star (K) : ', par%tstar
     else
        write(*,'(2a)')       'source spectrum file      : ', trim(par%source_spectrum)
     endif
     write(*,'(2a)')          'spectrum_type             : ', trim(par%spectrum_type)
     if (spec_derived) &
        write(*,'(a,es12.4)') 'derived luminosity (erg/s): ', par%luminosity
  endif

  deallocate(tb_lam, tb_alb, tb_cos, tb_cext, edge)

  !--- external illumination: build the external-field spectrum and, when the
  !--- mean intensity J is known, the absolute luminosity pi*J*A_surface.
  if (trim(par%source_geometry(1:8)) == 'external') call setup_sed_external()
  end subroutine setup_sed

  !---------------------------------------------------------------
  !--- External-field spectrum for SED mode with external illumination.
  !--- Spectrum-source priority: ext_spectrum file > ext_tstar Planck > the
  !--- global source spectrum (par%source_spectrum / par%tstar).  A physical
  !--- par%spectrum_type + ext_spectrum file carries the mean-intensity density
  !--- J_lambda; its bin integral is the band mean intensity J [erg/s/cm^2/sr].
  !--- When J is known (derived from the file, or set via par%ext_intensity) the
  !--- luminosity entering the grid is par%luminosity = pi*J*A_surface.
  subroutine setup_sed_external(compose_mode, lum_out)
  use mpi
  implicit none
  !--- compose_mode = .true.: the external field composes with an internal
  !--- source (par%source_geometry is internal, not 'external_*').  The
  !--- boundary surface area comes from par%ext_geometry and the derived
  !--- luminosity pi*J*A is returned in lum_out WITHOUT overwriting
  !--- par%luminosity (which the caller keeps as L_tot = L_int + L_ext).
  logical,       intent(in),  optional :: compose_mode
  real(kind=wp), intent(out), optional :: lum_out
  real(kind=wp), allocatable :: sp_lam(:), sp_lum(:)
  real(kind=wp) :: ext_sum, J_band, J_use, A_surf, lum_J, Lx, Ly, Lz
  integer       :: nsp, il, ierr
  logical       :: is_phys, J_defined, planck_ext, is_compose
  character(len=32) :: src_label

  is_compose = .false.
  if (present(compose_mode)) is_compose = compose_mode
  if (present(lum_out)) lum_out = 0.0_wp

  sed_ext_on = .true.
  allocate(sed_ext_lum(sed_nlam), sed_ext_pdf(sed_nlam), sed_ext_alias(sed_nlam), sed_ext_cdf(sed_nlam))
  sed_ext_lum(:) = 0.0_wp
  is_phys        = trim(par%spectrum_type) /= 'shape'
  J_defined      = .false.
  J_band         = 0.0_wp
  J_use          = 0.0_wp
  planck_ext     = .false.

  if (len_trim(par%ext_spectrum) > 0) then
     !--- explicit external spectrum file (columns in par%spectrum_type units).
     if (mpar%p_rank == 0) then
        call read_spectrum_file(trim(par%ext_spectrum), sp_lam, sp_lum, nsp)
        call convert_spectrum_units(sp_lam, sp_lum)
     endif
     call MPI_BCAST(nsp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     if (mpar%p_rank /= 0) allocate(sp_lam(nsp), sp_lum(nsp))
     call MPI_BCAST(sp_lam, nsp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(sp_lum, nsp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     do il = 1, sed_nlam
        if (sed_wave(il) < sp_lam(1) .or. sed_wave(il) > sp_lam(nsp)) then
           sed_ext_lum(il) = 0.0_wp
        else
           sed_ext_lum(il) = interp_clamped(log(sp_lam), sp_lum, log(sed_wave(il))) * sed_dwave(il)
        endif
     enddo
     deallocate(sp_lam, sp_lum)
     src_label = 'ext_spectrum file'
     if (is_phys) then
        !--- file columns are J_lambda densities: the bin integral is the band J.
        J_band = sum(sed_ext_lum)
        if (par%ext_intensity < -900.0_wp) then
           J_use = J_band          ! derive J from the file integral
        else
           J_use = par%ext_intensity  ! set scale: file used as shape only
        endif
        J_defined = .true.
     else if (par%ext_intensity > -900.0_wp) then
        J_use = par%ext_intensity;  J_defined = .true.
     endif
  else if (par%ext_tstar > 0.0_wp) then
     !--- Planck external field (a shape); J from par%ext_intensity if set.
     do il = 1, sed_nlam
        sed_ext_lum(il) = planck_shape(sed_wave(il), par%ext_tstar) * sed_dwave(il)
     enddo
     planck_ext = .true.
     src_label  = 'ext_tstar Planck'
     if (par%ext_intensity > -900.0_wp) then
        J_use = par%ext_intensity;  J_defined = .true.
     endif
  else
     !--- fall back to the global source spectrum (already built as sed_lum, a
     !--- normalized shape); J only from par%ext_intensity.
     if (len_trim(par%source_spectrum) > 0 .or. par%tstar > 0.0_wp) then
        sed_ext_lum(:) = sed_lum(:)
        src_label      = 'global source spectrum'
        if (par%ext_intensity > -900.0_wp) then
           J_use = par%ext_intensity;  J_defined = .true.
        endif
     else
        if (mpar%p_rank == 0) write(*,'(a)') &
           'ERROR: external SED illumination requires a spectrum: par%ext_spectrum (file), '// &
           'par%ext_tstar (Planck, K), or the global par%source_spectrum / par%tstar.'
        call MPI_FINALIZE(ierr);  stop
     endif
  endif

  ext_sum = sum(sed_ext_lum)
  if (.not. (ext_sum > 0.0_wp)) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: external SED spectrum has zero intensity on the wavelength grid.'
     call MPI_FINALIZE(ierr);  stop
  endif

  !--- surface area (cm^2) of the illuminated boundary; the luminosity entering
  !--- the grid for an isotropic external field of mean intensity J is pi*J*A.
  if (J_defined) then
     if (is_compose) then
        A_surf = ext_surface_area(par%ext_geometry)
     else
        select case (trim(par%source_geometry))
        case ('external_sph', 'external_sph_radial')
           A_surf = fourpi*(par%rmax*par%distance2cm)**2
        case ('external_rec')
           Lx = 2.0_wp*par%xmax*par%distance2cm
           Ly = 2.0_wp*par%ymax*par%distance2cm
           Lz = 2.0_wp*par%zmax*par%distance2cm
           A_surf = 2.0_wp*(Lx*Ly + Ly*Lz + Lx*Lz)
        case ('external_cyl')
           Lz = 2.0_wp*par%zmax*par%distance2cm
           A_surf = twopi*par%rmax*par%distance2cm*Lz + twopi*(par%rmax*par%distance2cm)**2
        case default
           A_surf = fourpi*(par%rmax*par%distance2cm)**2
        end select
     endif
     lum_J = pi * J_use * A_surf
     if (is_compose) then
        !--- compose mode: return L_ext; the caller owns par%luminosity (= L_tot).
        if (present(lum_out)) lum_out = lum_J
     else
        if (par%luminosity > -900.0_wp .and. &
            abs(par%luminosity - lum_J) > 1.0e-10_wp*abs(lum_J)) then
           if (mpar%p_rank == 0) write(*,'(a,es12.4,a,es12.4,a)') &
              'WARNING: par%luminosity (', par%luminosity, &
              ') conflicts with the J-based value (', lum_J, '); using the J-based value.'
        endif
        par%luminosity = lum_J
     endif
  else if (is_compose) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: composing an internal source with an external SED field requires a known '// &
        'mean intensity J (set par%ext_intensity, or a physical par%spectrum_type + par%ext_spectrum file).'
     call MPI_FINALIZE(ierr);  stop
  endif

  !--- normalize the external spectrum to a PDF and build the alias table.
  sed_ext_lum(:) = sed_ext_lum(:)/ext_sum
  sed_ext_pdf(:) = sed_ext_lum(:)
  call random_alias_setup(sed_ext_pdf, sed_ext_alias)

  !--- cumulative distribution for the quasi-random inverse-CDF sampler
  !--- (sample_ext_lambda_u); the alias table stays the default path.
  sed_ext_cdf(1) = sed_ext_lum(1)
  do il = 2, sed_nlam
     sed_ext_cdf(il) = sed_ext_cdf(il-1) + sed_ext_lum(il)
  enddo
  sed_ext_cdf(sed_nlam) = 1.0_wp

  if (mpar%p_rank == 0) then
     write(*,'(a)')          '--- external-field SED spectrum ---'
     write(*,'(2a)')         'external spectrum source  : ', trim(src_label)
     if (planck_ext) write(*,'(a,es12.4)') 'ext_tstar (K)             : ', par%ext_tstar
     if (J_defined) then
        if (is_phys .and. len_trim(par%ext_spectrum) > 0 .and. par%ext_intensity < -900.0_wp) then
           write(*,'(a,es12.4)') 'J_band derived (erg/s/cm2/sr): ', J_band
        endif
        write(*,'(a,es12.4)') 'J used   (erg/s/cm^2/sr)  : ', J_use
        write(*,'(a,es12.4)') 'A_surface (cm^2)          : ', A_surf
        write(*,'(a,es12.4)') 'luminosity pi*J*A (erg/s) : ', lum_J
     else
        write(*,'(a)')        'J unset (relative normalization; par%luminosity unchanged).'
     endif
  endif
  end subroutine setup_sed_external

  !---------------------------------------------------------------
  !--- illuminated boundary surface area [cm^2] for an external field placed on
  !--- a 'sph'|'cyl'|'rec' boundary (compose mode).  The luminosity entering the
  !--- grid for an isotropic field of mean intensity J is pi*J*A.
  function ext_surface_area(geom) result(A_surf)
  implicit none
  character(len=*), intent(in) :: geom
  real(kind=wp) :: A_surf, Lx, Ly, Lz
  select case (trim(geom))
  case ('sph')
     A_surf = fourpi*(par%rmax*par%distance2cm)**2
  case ('rec')
     Lx = 2.0_wp*par%xmax*par%distance2cm
     Ly = 2.0_wp*par%ymax*par%distance2cm
     Lz = 2.0_wp*par%zmax*par%distance2cm
     A_surf = 2.0_wp*(Lx*Ly + Ly*Lz + Lx*Lz)
  case ('cyl')
     Lz = 2.0_wp*par%zmax*par%distance2cm
     A_surf = twopi*par%rmax*par%distance2cm*Lz + twopi*(par%rmax*par%distance2cm)**2
  case default
     A_surf = fourpi*(par%rmax*par%distance2cm)**2
  end select
  end function ext_surface_area

  !---------------------------------------------------------------
  !--- sample a wavelength-bin index from the source spectrum.
  function sample_sed_lambda() result(il)
  implicit none
  integer :: il
  il = rand_alias_choise(sed_src_pdf, sed_src_alias)
  end function sample_sed_lambda

  !---------------------------------------------------------------
  !--- sample a wavelength-bin index from the external-field spectrum.
  function sample_ext_lambda() result(il)
  implicit none
  integer :: il
  il = rand_alias_choise(sed_ext_pdf, sed_ext_alias)
  end function sample_ext_lambda

  !---------------------------------------------------------------
  !--- inverse-CDF wavelength-bin sampler for the quasi-random launch: return
  !--- the smallest bin il with uf <= cdf(il).  uf is a uniform in (0,1).  The
  !--- monotone inverse CDF keeps a stratified Sobol coordinate stratified in
  !--- the wavelength bins (the alias table would scramble it).
  function sample_sed_lambda_u(uf) result(il)
  implicit none
  real(kind=wp), intent(in) :: uf
  integer :: il
  il = cdf_search(sed_src_cdf, uf)
  end function sample_sed_lambda_u

  function sample_ext_lambda_u(uf) result(il)
  implicit none
  real(kind=wp), intent(in) :: uf
  integer :: il
  il = cdf_search(sed_ext_cdf, uf)
  end function sample_ext_lambda_u

  !---------------------------------------------------------------
  !--- binary search: smallest index il with uf <= cdf(il) (cdf ascending,
  !--- cdf(n) = 1).
  pure function cdf_search(cdf, uf) result(il)
  implicit none
  real(kind=wp), intent(in) :: cdf(:), uf
  integer :: il, lo, hi, mid
  lo = 1;  hi = size(cdf)
  do while (lo < hi)
     mid = (lo + hi)/2
     if (uf <= cdf(mid)) then;  hi = mid;  else;  lo = mid + 1;  endif
  enddo
  il = lo
  end function cdf_search

  !---------------------------------------------------------------
  !--- Planck function B_lambda (arbitrary normalization), lambda in um.
  pure function planck_shape(lam_um, T) result(b)
  implicit none
  real(kind=wp), intent(in) :: lam_um, T
  real(kind=wp) :: b, x
  real(kind=wp), parameter :: hc_over_k = 1.43877687750393e4_wp  ! [um K]
  x = hc_over_k/(lam_um*T)
  if (x > 700.0_wp) then
     b = 0.0_wp
  else
     b = 1.0_wp/(lam_um**5 * (exp(x) - 1.0_wp))
  endif
  end function planck_shape

  !---------------------------------------------------------------
  !--- linear interpolation with clamping at the table ends.
  pure function interp_clamped(x, y, xnew) result(ynew)
  implicit none
  real(kind=wp), intent(in) :: x(:), y(:), xnew
  real(kind=wp) :: ynew, t
  integer :: n, i
  n = size(x)
  if (xnew <= x(1)) then
     ynew = y(1)
  else if (xnew >= x(n)) then
     ynew = y(n)
  else
     do i = 2, n
        if (xnew <= x(i)) exit
     enddo
     t    = (xnew - x(i-1))/(x(i) - x(i-1))
     ynew = y(i-1) + t*(y(i) - y(i-1))
  endif
  end function interp_clamped

  !---------------------------------------------------------------
  !--- read the extinction table: comment lines start with '#';
  !--- columns: lambda[um]  albedo  <cos>  C_ext/H [cm^2/H]  (extra columns ignored).
  subroutine read_kext_table(fname, lam, alb, cosg, cext, n)
  implicit none
  character(len=*), intent(in) :: fname
  real(kind=wp), allocatable, intent(out) :: lam(:), alb(:), cosg(:), cext(:)
  integer, intent(out) :: n
  character(len=512) :: line
  integer :: unit, ios, i
  open(newunit=unit, file=fname, status='old', action='read')
  n = 0
  do
     read(unit,'(a)',iostat=ios) line
     if (ios /= 0) exit
     line = adjustl(line)
     if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
     n = n + 1
  enddo
  if (n < 2) then
     write(*,'(3a)') 'ERROR: extinction table ', trim(fname), ' has fewer than 2 data rows.'
     stop
  endif
  allocate(lam(n), alb(n), cosg(n), cext(n))
  rewind(unit)
  i = 0
  do
     read(unit,'(a)',iostat=ios) line
     if (ios /= 0) exit
     line = adjustl(line)
     if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
     i = i + 1
     read(line,*) lam(i), alb(i), cosg(i), cext(i)
  enddo
  close(unit)
  end subroutine read_kext_table

  !---------------------------------------------------------------
  !--- read a 2-column spectrum file: lambda[um], L_lambda[arbitrary].
  subroutine read_spectrum_file(fname, lam, lum, n)
  implicit none
  character(len=*), intent(in) :: fname
  real(kind=wp), allocatable, intent(out) :: lam(:), lum(:)
  integer, intent(out) :: n
  character(len=512) :: line
  integer :: unit, ios, i
  open(newunit=unit, file=fname, status='old', action='read')
  n = 0
  do
     read(unit,'(a)',iostat=ios) line
     if (ios /= 0) exit
     line = adjustl(line)
     if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
     n = n + 1
  enddo
  if (n < 2) then
     write(*,'(3a)') 'ERROR: spectrum file ', trim(fname), ' has fewer than 2 data rows.'
     stop
  endif
  allocate(lam(n), lum(n))
  rewind(unit)
  i = 0
  do
     read(unit,'(a)',iostat=ios) line
     if (ios /= 0) exit
     line = adjustl(line)
     if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
     i = i + 1
     read(line,*) lam(i), lum(i)
  enddo
  close(unit)
  end subroutine read_spectrum_file

  !---------------------------------------------------------------
  !--- convert a source spectrum file from par%spectrum_type column units to
  !--- the internal (lambda[um], L_lambda per um) representation.  'shape' and
  !--- 'per_um' are already in the internal units; the frequency/energy types
  !--- invert the abscissa order, so the arrays are reversed to ascending lambda.
  !--- Called on rank 0 immediately after read_spectrum_file (before broadcast).
  subroutine convert_spectrum_units(lam, lum)
  implicit none
  real(kind=wp), allocatable, intent(inout) :: lam(:), lum(:)
  real(kind=wp), allocatable :: tmp(:)
  integer :: n, i
  n = size(lam)
  select case (trim(par%spectrum_type))
  case ('shape', 'per_um')
     return
  case ('per_ang')
     !--- lambda [A] -> [um], L_lambda per A -> per um.
     lam(:) = lam(:)*1.0e-4_wp
     lum(:) = lum(:)*1.0e4_wp
  case ('per_hz')
     !--- nu [Hz] -> lambda [um]; L_nu -> L_lambda = L_nu * c/lambda^2.
     lam(:) = c_um/lam(:)
     lum(:) = lum(:)*c_um/lam(:)**2
  case ('per_ev')
     !--- E [eV] -> lambda [um]; L_E -> L_lambda = L_E * hc/lambda^2.
     lam(:) = hc_evum/lam(:)
     lum(:) = lum(:)*hc_evum/lam(:)**2
  end select
  !--- per_hz/per_ev flip the abscissa order: reverse to ascending lambda.
  if (trim(par%spectrum_type) == 'per_hz' .or. trim(par%spectrum_type) == 'per_ev') then
     allocate(tmp(n))
     tmp(:) = lam(n:1:-1);  lam(:) = tmp(:)
     tmp(:) = lum(n:1:-1);  lum(:) = tmp(:)
     deallocate(tmp)
  endif
  !--- verify the converted grid is strictly ascending in lambda.
  do i = 2, n
     if (lam(i) <= lam(i-1)) then
        write(*,'(2a)') 'ERROR: source spectrum not strictly ascending in lambda ', &
           'after unit conversion (spectrum_type='//trim(par%spectrum_type)//').'
        stop
     endif
  enddo
  end subroutine convert_spectrum_units

end module sed_mod
