module sed_mod
!--- Multi-wavelength (SED) infrastructure for MoCafe v2.00 (Stage 1).
!--- Provides: a log-spaced wavelength grid; wavelength-dependent dust
!--- properties C_ext(lambda), albedo(lambda), g(lambda); and a stellar source
!--- spectrum (Planck or 2-column file) from which every photon draws a
!--- continuous wavelength.
!---
!--- The cross sections come from the grain model named by par%dust_model,
!--- integrated over its size distribution by grain_model_mod.  They then refer
!--- to the same grains that radiate the absorbed energy back out in
!--- dustemis_mod, which is what makes the energy balance of a cell physical.
!--- An explicit table file (par%kext_file) overrides the model when set; the
!--- two agree only if the file was produced from the same model.
!---
!--- The wavelength is sampled by inverting the analytic cumulative distribution
!--- of the spectrum (spectrum_sampler_mod), not by drawing a bin index.  Pinning
!--- a photon to a bin center would erase the spectral structure inside the bin
!--- and freeze the dust cross sections at their bin-center values.  The bin
!--- index sed_bin_of(lambda) is still carried by the photon, because the output
!--- images and the radiation-field tally are binned in wavelength, but every
!--- optical quantity is evaluated at the sampled wavelength itself.
!---
!--- Transport of a photon at wavelength lambda uses the grey-rescaling factor
!--- sed_sext_at(lambda) = C_ext(lambda)/C_ext(lambda_ref) (photon%s_ext):
!--- the grid rhokap is the opacity at the reference wavelength, and every
!--- optical depth is scaled by s_ext (same mathematics as the Jonsson 2006
!--- tau scan), so the raytrace routines (car/clump/amr) need no changes.
  use define
  use random, only : rand_number
  use spectrum_sampler_mod
  use grain_model_mod, only : grain_extinction_table
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
  real(kind=wp), allocatable :: sed_edge(:)     ! bin edges [um], size sed_nlam+1
  real(kind=wp) :: sed_cext_ref = 0.0_wp        ! C_ext/H at par%lambda_ref
  real(kind=wp) :: sed_dlnlam   = 0.0_wp        ! ln(lambda) spacing of the grid

  !--- dust extinction table, kept so the cross sections can be evaluated at an
  !--- arbitrary wavelength and not only at the bin centers.
  integer :: sed_ntab = 0
  real(kind=wp), allocatable :: sed_tab_lam(:), sed_tab_alb(:)
  real(kind=wp), allocatable :: sed_tab_cos(:), sed_tab_cext(:)

  !--- continuous wavelength samplers for the source and external spectra.
  type(spectrum_sampler_type) :: sed_src_spectrum
  type(spectrum_sampler_type) :: sed_ext_spectrum

  !--- external-field spectrum (SED mode with external illumination).
  logical :: sed_ext_on = .false.
  real(kind=wp), allocatable :: sed_ext_lum(:)   ! external spectrum fraction per bin (sum = 1)

  !--- unit-conversion constants for physical par%spectrum_type files.
  real(kind=wp), parameter :: c_um     = 2.99792458e14_wp  ! speed of light [um/s]
  real(kind=wp), parameter :: hc_evum  = 1.23984193_wp     ! h*c [eV um]

contains
  !---------------------------------------------------------------
  subroutine setup_sed()
  use mpi
  implicit none
  ! local variables
  real(kind=wp), allocatable :: sp_lam(:), sp_lum(:)
  real(kind=wp) :: lum_sum
  integer       :: nsp, il, ierr
  logical       :: spec_is_absolute, spec_derived, ok, from_grain_model

  if (par%nlambda < 2) then
     if (mpar%p_rank == 0) write(*,'(a)') 'ERROR: par%nlambda must be >= 2 in SED mode.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (.not. (par%lambda_max > par%lambda_min .and. par%lambda_min > 0.0_wp)) then
     if (mpar%p_rank == 0) write(*,'(a)') 'ERROR: require 0 < lambda_min < lambda_max in SED mode.'
     call MPI_FINALIZE(ierr);  stop
  endif
  !--- The transport cross sections come either from the named grain model
  !--- (par%dust_model, the default) or from an explicit table file.
  !--- Deriving them from the model keeps the absorbed power and the reemitted
  !--- spectrum on the same grains and on the same wavelength grid; a file is
  !--- read only when par%kext_file is set, which then overrides the model.
  from_grain_model = len_trim(par%kext_file) == 0

  !--- log-spaced wavelength grid (bin edges and geometric bin centers).
  sed_nlam = par%nlambda
  allocate(sed_edge(sed_nlam+1))
  allocate(sed_wave(sed_nlam), sed_dwave(sed_nlam))
  allocate(sed_cext(sed_nlam), sed_albedo(sed_nlam), sed_hgg(sed_nlam), sed_sext(sed_nlam))
  allocate(sed_lum(sed_nlam))
  sed_dlnlam = log(par%lambda_max/par%lambda_min)/sed_nlam
  do il = 1, sed_nlam+1
     sed_edge(il) = par%lambda_min * exp((il-1)*sed_dlnlam)
  enddo
  sed_edge(sed_nlam+1) = par%lambda_max
  do il = 1, sed_nlam
     sed_wave(il)  = sqrt(sed_edge(il)*sed_edge(il+1))
     sed_dwave(il) = sed_edge(il+1) - sed_edge(il)
  enddo

  !--- the dust extinction table, kept for the whole run: the cross sections are
  !--- evaluated at each photon's sampled wavelength, not only at the bin centers.
  if (from_grain_model) then
     !--- size-integrated C_ext/albedo/<cos> of the named model.  Every rank
     !--- builds the model from the same files, so no broadcast is needed.
     call grain_extinction_table(sed_tab_lam, sed_tab_alb, sed_tab_cos, sed_tab_cext, sed_ntab)
  else
     !--- explicit table (rank 0 reads and broadcasts).
     if (mpar%p_rank == 0) call read_kext_table(trim(par%kext_file), sed_tab_lam, sed_tab_alb, &
                                                sed_tab_cos, sed_tab_cext, sed_ntab)
     call MPI_BCAST(sed_ntab, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     if (mpar%p_rank /= 0) allocate(sed_tab_lam(sed_ntab), sed_tab_alb(sed_ntab), &
                                    sed_tab_cos(sed_ntab), sed_tab_cext(sed_ntab))
     call MPI_BCAST(sed_tab_lam,  sed_ntab, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(sed_tab_alb,  sed_ntab, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(sed_tab_cos,  sed_ntab, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(sed_tab_cext, sed_ntab, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  endif

  if (mpar%p_rank == 0) then
     if (from_grain_model) then
        write(*,'(3a,i0,a)') 'Dust optics               : ', trim(par%dust_model), &
                             ' (SEDust, ', sed_ntab, ' wavelengths)'
     else
        write(*,'(2a)')      'Dust optics               : ', trim(par%kext_file)
     endif
  endif

  !--- bin-center values of the dust properties, kept for the dust-emission
  !--- arrays and the output writers.  Photons use sed_*_at(lambda) instead.
  do il = 1, sed_nlam
     sed_cext(il)   = sed_cext_at(sed_wave(il))
     sed_albedo(il) = sed_albedo_at(sed_wave(il))
     sed_hgg(il)    = sed_hgg_at(sed_wave(il))
  enddo
  sed_cext_ref  = sed_cext_at(par%lambda_ref)
  sed_sext(:)   = sed_cext(:)/sed_cext_ref
  !--- the grid opacity (rhokap) is defined at the reference wavelength.
  par%cext_dust = sed_cext_ref
  par%albedo    = sed_albedo_at(par%lambda_ref)
  par%hgg       = sed_hgg_at(par%lambda_ref)
  par%lambda0   = par%lambda_ref

  !--- source spectrum.  With multiple source components (par%nsource > 1) the
  !--- spectra for each source are set up in sources_mod, so the global single-source
  !--- spectrum is optional here: build a flat placeholder and skip the checks.
  ok = .false.
  if (par%nsource > 1 .and. len_trim(par%source_spectrum) == 0 .and. par%tstar <= 0.0_wp) then
     call flat_spectrum_over_grid(sed_src_spectrum, ok)
  else if (len_trim(par%source_spectrum) > 0) then
     if (mpar%p_rank == 0) then
        call read_spectrum_file(trim(par%source_spectrum), sp_lam, sp_lum, nsp)
        call convert_spectrum_units(sp_lam, sp_lum)
     endif
     call MPI_BCAST(nsp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     if (mpar%p_rank /= 0) allocate(sp_lam(nsp), sp_lum(nsp))
     call MPI_BCAST(sp_lam, nsp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(sp_lum, nsp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     !--- the support is the overlap of the grid with the tabulated range: the
     !--- spectrum carries no luminosity where it is not tabulated.
     call tabulated_spectrum_sampler(sed_src_spectrum, sp_lam, sp_lum, &
                                     par%lambda_min, par%lambda_max, sed_edge, ok)
     deallocate(sp_lam, sp_lum)
  else if (par%tstar > 0.0_wp) then
     call planck_spectrum_sampler(sed_src_spectrum, par%tstar, &
                                  par%lambda_min, par%lambda_max, sed_edge, ok)
  else if (trim(par%source_geometry(1:8)) == 'external') then
     !--- external illumination: the internal (point/extended) source spectrum
     !--- is unused; the external-field spectrum is set in setup_sed_external.
     !--- Build a flat placeholder so the checks below stay valid.
     call flat_spectrum_over_grid(sed_src_spectrum, ok)
  else
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: SED mode requires a source spectrum: par%source_spectrum (file) or par%tstar (Planck, K).'
     call MPI_FINALIZE(ierr);  stop
  endif

  if (.not. ok) then
     if (mpar%p_rank == 0) write(*,'(a)') 'ERROR: source spectrum has zero luminosity on the wavelength grid.'
     call MPI_FINALIZE(ierr);  stop
  endif

  !--- luminosity fraction of each bin, an exact integral of the spectrum over
  !--- the bin: the bin edges are nodes of the sampler.
  do il = 1, sed_nlam
     sed_lum(il) = band_luminosity_fraction(sed_src_spectrum, sed_edge(il), sed_edge(il+1))
  enddo
  lum_sum = sed_src_spectrum%total

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
     write(*,'(a,i8)')        'spectrum sampling nodes   : ', size(sed_src_spectrum%lam)
  endif

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
  logical       :: is_phys, J_defined, planck_ext, is_compose, ok
  character(len=32) :: src_label

  is_compose = .false.
  if (present(compose_mode)) is_compose = compose_mode
  if (present(lum_out)) lum_out = 0.0_wp

  sed_ext_on = .true.
  allocate(sed_ext_lum(sed_nlam))
  sed_ext_lum(:) = 0.0_wp
  is_phys        = trim(par%spectrum_type) /= 'shape'
  J_defined      = .false.
  J_band         = 0.0_wp
  J_use          = 0.0_wp
  planck_ext     = .false.
  ok             = .false.

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
     call tabulated_spectrum_sampler(sed_ext_spectrum, sp_lam, sp_lum, &
                                     par%lambda_min, par%lambda_max, sed_edge, ok)
     deallocate(sp_lam, sp_lum)
     src_label = 'ext_spectrum file'
     if (is_phys .and. ok) then
        !--- file columns are J_lambda densities: their integral over the grid is
        !--- the band mean intensity J.
        J_band = sed_ext_spectrum%total
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
     call planck_spectrum_sampler(sed_ext_spectrum, par%ext_tstar, &
                                  par%lambda_min, par%lambda_max, sed_edge, ok)
     planck_ext = .true.
     src_label  = 'ext_tstar Planck'
     if (par%ext_intensity > -900.0_wp) then
        J_use = par%ext_intensity;  J_defined = .true.
     endif
  else
     !--- reuse the global source spectrum as the external shape; J only from
     !--- par%ext_intensity.
     if (len_trim(par%source_spectrum) > 0 .or. par%tstar > 0.0_wp) then
        sed_ext_spectrum = sed_src_spectrum
        ok             = sampler_is_ready(sed_ext_spectrum)
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

  if (.not. ok) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: external SED spectrum has zero intensity on the wavelength grid.'
     call MPI_FINALIZE(ierr);  stop
  endif
  ext_sum = sed_ext_spectrum%total
  !--- intensity fraction of each bin, an exact integral over the bin.
  do il = 1, sed_nlam
     sed_ext_lum(il) = band_luminosity_fraction(sed_ext_spectrum, sed_edge(il), sed_edge(il+1))
  enddo

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
  !--- draw a wavelength from the source spectrum [um].
  function sample_sed_wavelength() result(lam)
  implicit none
  real(kind=wp) :: lam
  lam = sample_wavelength(sed_src_spectrum, rand_number())
  end function sample_sed_wavelength

  !---------------------------------------------------------------
  !--- draw a wavelength from the external-field spectrum [um].
  function sample_ext_wavelength() result(lam)
  implicit none
  real(kind=wp) :: lam
  lam = sample_wavelength(sed_ext_spectrum, rand_number())
  end function sample_ext_wavelength

  !---------------------------------------------------------------
  !--- the same draws from a caller-supplied uniform, for the quasi-random
  !--- launch.  Inverting a monotone cumulative distribution maps a stratified
  !--- Sobol coordinate onto a stratified wavelength, which is why the launch
  !--- takes this route rather than any order-scrambling sampler.
  function sample_sed_wavelength_u(uf) result(lam)
  implicit none
  real(kind=wp), intent(in) :: uf
  real(kind=wp) :: lam
  lam = sample_wavelength(sed_src_spectrum, uf)
  end function sample_sed_wavelength_u

  function sample_ext_wavelength_u(uf) result(lam)
  implicit none
  real(kind=wp), intent(in) :: uf
  real(kind=wp) :: lam
  lam = sample_wavelength(sed_ext_spectrum, uf)
  end function sample_ext_wavelength_u

  !---------------------------------------------------------------
  !--- index of the wavelength bin holding lambda.  The grid is log-spaced with
  !--- constant sed_dlnlam, so the bin follows in closed form; lambda_max lands
  !--- in the last bin rather than one past it.
  pure function sed_bin_of(lam) result(il)
  implicit none
  real(kind=wp), intent(in) :: lam
  integer :: il
  il = floor(log(lam/par%lambda_min)/sed_dlnlam) + 1
  if (il < 1)        il = 1
  if (il > sed_nlam) il = sed_nlam
  end function sed_bin_of

  !---------------------------------------------------------------
  !--- dust cross sections at an arbitrary wavelength, read from the extinction
  !--- table with the model the grid values use: C_ext log-log, albedo and g
  !--- linear in ln(lambda), clamped outside the tabulated range.
  pure function sed_cext_at(lam) result(cext)
  implicit none
  real(kind=wp), intent(in) :: lam
  real(kind=wp) :: cext
  cext = exp(interp_clamped(log(sed_tab_lam), log(sed_tab_cext), log(lam)))
  end function sed_cext_at

  pure function sed_sext_at(lam) result(sext)
  implicit none
  real(kind=wp), intent(in) :: lam
  real(kind=wp) :: sext
  sext = sed_cext_at(lam)/sed_cext_ref
  end function sed_sext_at

  pure function sed_albedo_at(lam) result(alb)
  implicit none
  real(kind=wp), intent(in) :: lam
  real(kind=wp) :: alb
  alb = interp_clamped(log(sed_tab_lam), sed_tab_alb, log(lam))
  end function sed_albedo_at

  pure function sed_hgg_at(lam) result(g)
  implicit none
  real(kind=wp), intent(in) :: lam
  real(kind=wp) :: g
  g = interp_clamped(log(sed_tab_lam), sed_tab_cos, log(lam))
  end function sed_hgg_at

  !---------------------------------------------------------------
  !--- a spectrum flat in wavelength across the whole grid, used where the
  !--- luminosity of this component carries no spectral information.
  subroutine flat_spectrum_over_grid(sampler, ok)
  implicit none
  type(spectrum_sampler_type), intent(out) :: sampler
  logical, intent(out) :: ok
  real(kind=wp) :: lam2(2), dens2(2)
  lam2  = [par%lambda_min, par%lambda_max]
  dens2 = [1.0_wp, 1.0_wp]
  call tabulated_spectrum_sampler(sampler, lam2, dens2, &
                                  par%lambda_min, par%lambda_max, sed_edge, ok)
  end subroutine flat_spectrum_over_grid

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
