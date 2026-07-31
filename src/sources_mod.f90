module sources_mod
!--- Multiple stellar source components (MoCafe v2.00, Stage 6).  Each
!--- component has its own spectrum (Planck or 2-column file), luminosity,
!--- geometry and geometry parameters.  Sources are sampled in proportion to
!--- their luminosity; each carries its own continuous wavelength sampler so a
!--- run mixes, e.g., a hot young population and a cool old population with
!--- different spatial distributions.  When par%nsource == 1 this reduces to
!--- the single-source SED path (sed_mod) and is not activated.
  use define
  use random,  only : rand_number, rand_gauss, &
                      rand_zexp, rand_sech2, rand_r1exp
  use sed_mod, only : sed_nlam, sed_edge, sed_bin_of, sed_sext_at, sed_albedo_at, sed_hgg_at, &
                      read_spectrum_file, convert_spectrum_units, cdf_search
  use spectrum_sampler_mod, only : spectrum_sampler_type, tabulated_spectrum_sampler, &
                      planck_spectrum_sampler, sample_wavelength, band_luminosity_fraction
  use random_bulge,  only : rand_sersic, rand_boxy, rand_bar, rand_xbar
  implicit none
  private
  public :: setup_sources, gen_source_photon, gen_source_photon_qmc, use_sources

  logical :: use_sources = .false.
  integer :: nsrc = 0
  real(kind=wp), allocatable :: src_lum_cdf(:)          ! (nsrc) cumulative luminosity, normalized
  type(spectrum_sampler_type), allocatable :: src_spectrum(:)  ! (nsrc) continuous wavelength sampler
  real(kind=wp), allocatable :: src_lumfrac(:,:)        ! (nlam, nsrc) luminosity fraction per bin (for output)
  real(kind=wp), allocatable :: src_Lpacket(:)          ! (nsrc) energy per packet for this source

contains
  !---------------------------------------------------------------
  subroutine setup_sources()
  use mpi
  implicit none
  real(kind=wp), allocatable :: lum(:), lam_f(:), lum_f(:)
  real(kind=wp) :: lsum, lnorm
  integer :: is, il, nsp, ierr
  logical :: is_phys, absolute_is, ok
  logical, allocatable :: derived(:)

  nsrc = par%nsource
  if (nsrc <= 1) then
     use_sources = .false.;  return
  endif

  !--- monochromatic (non-SED) multiple internal sources: luminosity-weighted
  !--- component selection + geometry-based position sampling only, with no
  !--- wavelength/spectrum concept.  The SED-only members (src_spectrum,
  !--- src_lumfrac) are left unallocated and never referenced on this path.
  if (.not. par%use_sed) then
     call setup_sources_mono()
     return
  endif

  use_sources = .true.
  is_phys = trim(par%spectrum_type) /= 'shape'

  allocate(src_lum_cdf(nsrc), src_spectrum(nsrc))
  allocate(src_lumfrac(sed_nlam,nsrc), src_Lpacket(nsrc), lum(nsrc))
  allocate(derived(nsrc));  derived(:) = .false.

  !--- spectrum of each source -> continuous wavelength sampler + luminosity, and
  !--- the source luminosity.  For a physical-type (absolute) file spectrum the
  !--- luminosity is the file integral lnorm [erg/s] when src_lum is unset, or
  !--- src_lum when set (rescales the file).  For a 'shape' file or a Planck
  !--- source the luminosity is src_lum when set, else the equal split of
  !--- par%luminosity (legacy); mixing an unset shape/Planck source with an
  !--- absolute file source is ambiguous and rejected.
  do is = 1, nsrc
     absolute_is = is_phys .and. len_trim(par%src_spectrum(is)) > 0
     if (len_trim(par%src_spectrum(is)) > 0) then
        if (mpar%p_rank == 0) then
           call read_spectrum_file(trim(par%src_spectrum(is)), lam_f, lum_f, nsp)
           call convert_spectrum_units(lam_f, lum_f)
        endif
        call MPI_BCAST(nsp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if (mpar%p_rank /= 0) allocate(lam_f(nsp), lum_f(nsp))
        call MPI_BCAST(lam_f, nsp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(lum_f, nsp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call tabulated_spectrum_sampler(src_spectrum(is), lam_f, lum_f, &
                                        par%lambda_min, par%lambda_max, sed_edge, ok)
        deallocate(lam_f, lum_f)
     else if (par%src_tstar(is) > 0.0_wp) then
        call planck_spectrum_sampler(src_spectrum(is), par%src_tstar(is), &
                                     par%lambda_min, par%lambda_max, sed_edge, ok)
     else
        if (mpar%p_rank == 0) write(*,'(a,i0,a)') &
           'ERROR: source ', is, ' needs src_tstar or src_spectrum.'
        call MPI_FINALIZE(ierr);  stop
     endif
     if (.not. ok) then
        if (mpar%p_rank == 0) write(*,'(a,i0,a)') 'ERROR: source ', is, ' has zero luminosity on the grid.'
        call MPI_FINALIZE(ierr);  stop
     endif
     lnorm = src_spectrum(is)%total

     !--- source luminosity.
     if (absolute_is) then
        if (par%src_lum(is) > 0.0_wp) then
           lum(is) = par%src_lum(is)      ! rescale the absolute file to src_lum
        else
           lum(is) = lnorm                ! derive [erg/s] from the file integral
           derived(is) = .true.
        endif
     else
        if (par%src_lum(is) > 0.0_wp) then
           lum(is) = par%src_lum(is)
        else if (.not. is_phys) then
           lum(is) = par%luminosity/dble(nsrc)   ! legacy equal split
        else
           if (mpar%p_rank == 0) write(*,'(a,i0,a)') 'ERROR: source ', is, &
              ': set src_lum for Planck/shape components when mixing with absolute file spectra.'
           call MPI_FINALIZE(ierr);  stop
        endif
     endif

     !--- luminosity fraction of each bin, an exact integral over the bin.
     do il = 1, sed_nlam
        src_lumfrac(il,is) = band_luminosity_fraction(src_spectrum(is), sed_edge(il), sed_edge(il+1))
     enddo
  enddo
  lsum = sum(lum)
  par%luminosity = lsum          ! total luminosity is the sum of components

  !--- luminosity CDF over sources and packet energy for each source.  Each source
  !--- emits (nphotons * lum(is)/lsum) packets, so a packet carries lsum/nphotons.
  src_lum_cdf(1) = lum(1)
  do is = 2, nsrc
     src_lum_cdf(is) = src_lum_cdf(is-1) + lum(is)
  enddo
  src_lum_cdf(:) = src_lum_cdf(:)/lsum
  src_Lpacket(:) = lsum/dble(par%nphotons)

  if (mpar%p_rank == 0) then
     write(*,'(a)')      '--- Multiple stellar source components (Stage 6) ---'
     write(*,'(a,i0)')   'N sources          : ', nsrc
     write(*,'(2a)')     'spectrum_type      : ', trim(par%spectrum_type)
     do is = 1, nsrc
        write(*,'(a,i2,4a,es11.3,a,es11.3,a)') '  src ', is, '  geom=', trim(par%src_geometry(is)), &
           '  spec=', merge('Planck', 'file  ', par%src_tstar(is) > 0.0_wp), &
           lum(is), '  T/-=', par%src_tstar(is), &
           merge('  (derived)', '           ', derived(is))
     enddo
     write(*,'(a,es12.4)') 'total luminosity   : ', lsum
  endif

  deallocate(lum, derived)
  end subroutine setup_sources

  !---------------------------------------------------------------
  subroutine setup_sources_mono()
  !--- monochromatic multiple internal sources.  Each component has its own
  !--- geometry and luminosity; components are sampled in proportion to their
  !--- luminosity.  Grey dust (par%albedo, par%hgg) is used for every packet.
  use mpi
  implicit none
  real(kind=wp), allocatable :: lum(:)
  real(kind=wp) :: lsum
  integer :: is, ierr

  !--- gen_source_photon does not carry Stokes vectors, so multiple internal
  !--- sources cannot be combined with the Stokes scattering path.
  if (par%use_stokes) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: monochromatic multiple internal sources (nsource>1) do not support par%use_stokes.'
     call MPI_FINALIZE(ierr);  stop
  endif

  allocate(src_lum_cdf(nsrc), src_Lpacket(nsrc), lum(nsrc))

  !--- component luminosity: par%src_lum when set, else the equal split of
  !--- par%luminosity (same convention as the SED shape/Planck path).
  do is = 1, nsrc
     if (par%src_lum(is) > 0.0_wp) then
        lum(is) = par%src_lum(is)
     else
        lum(is) = par%luminosity/dble(nsrc)
     endif
  enddo
  lsum = sum(lum)
  if (.not. (lsum > 0.0_wp)) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: monochromatic multiple internal sources have zero total luminosity.'
     call MPI_FINALIZE(ierr);  stop
  endif
  par%luminosity = lsum          ! total luminosity is the sum of components

  !--- luminosity CDF over components and packet energy for each component.
  src_lum_cdf(1) = lum(1)
  do is = 2, nsrc
     src_lum_cdf(is) = src_lum_cdf(is-1) + lum(is)
  enddo
  src_lum_cdf(:) = src_lum_cdf(:)/lsum
  src_Lpacket(:) = lsum/dble(par%nphotons)

  use_sources = .true.

  if (mpar%p_rank == 0) then
     write(*,'(a)')    '--- Multiple internal sources (monochromatic) ---'
     write(*,'(a,i0)') 'N sources          : ', nsrc
     do is = 1, nsrc
        write(*,'(a,i2,3a,es11.3)') '  src ', is, '  geom=', &
           trim(par%src_geometry(is)), '  lum=', lum(is)
     enddo
     write(*,'(a,es12.4)') 'total luminosity   : ', lsum
  endif

  deallocate(lum)
  end subroutine setup_sources_mono

  !---------------------------------------------------------------
  !--- generate a photon from a luminosity-weighted source component: pick the
  !--- source, sample its position (per geometry), an isotropic direction, and
  !--- a wavelength from that source's spectrum.
  subroutine gen_source_photon(grid, photon)
  implicit none
  type(grid_type),   intent(in)  :: grid
  type(photon_type), intent(inout) :: photon
  real(kind=wp) :: u, sint, cost, phi, rp, rs_max, tanp, bx, by, bz
  integer :: is, lo, hi, mid

  !--- select source by luminosity CDF.
  u = rand_number()
  lo = 1;  hi = nsrc
  do while (lo < hi)
     mid = (lo+hi)/2
     if (u <= src_lum_cdf(mid)) then;  hi = mid;  else;  lo = mid+1;  endif
  enddo
  is = lo

  !--- position from this source's geometry.
  select case (trim(par%src_geometry(is)))
  case ('uniform')
     rp   = rand_number()**(1.0_wp/3.0_wp) * par%rmax
     cost = 2.0_wp*rand_number()-1.0_wp;  sint = sqrt(1.0_wp-cost*cost)
     phi  = twopi*rand_number()
     photon%x = rp*sint*cos(phi);  photon%y = rp*sint*sin(phi);  photon%z = rp*cost
  case ('gaussian')
     photon%x = grid%xrange*rand_number()+grid%xmin
     photon%y = grid%yrange*rand_number()+grid%ymin
     photon%z = par%src_zscale(is)/sqrt(2.0_wp)*rand_gauss()
  case ('exponential', 'sech', 'exp_spiral')
     !--- radially exponential disk (r ~ r*exp(-r/rscale)); vertically
     !--- exponential ('exponential') or sech^2 ('sech'); 'exp_spiral' adds
     !--- log-spiral arms by rejection.  Falls back to a plane-uniform disk when
     !--- src_rscale is not set (backward compatible).
     if (par%src_rscale(is) > 0.0_wp) then
        if (trim(par%src_geometry(is)) == 'exp_spiral' .and. par%spiral_m > 0) then
           tanp = tan(par%spiral_pitch*pi/180.0_wp)
           do
              rp  = par%src_rscale(is) * rand_r1exp(par%rmax/par%src_rscale(is))
              phi = twopi*rand_number()
              if (rp <= 0.0_wp) cycle
              if ((1.0_wp+par%spiral_amp)*rand_number() <= &
                  1.0_wp + par%spiral_amp*sin(par%spiral_m*(log(rp)/tanp - phi))) exit
           enddo
        else
           rp  = par%src_rscale(is) * rand_r1exp(par%rmax/par%src_rscale(is))
           phi = twopi*rand_number()
        endif
        photon%x = rp*cos(phi);  photon%y = rp*sin(phi)
     else
        photon%x = grid%xrange*rand_number()+grid%xmin
        photon%y = grid%yrange*rand_number()+grid%ymin
     endif
     if (trim(par%src_geometry(is)) == 'sech') then
        photon%z = par%src_zscale(is)*rand_sech2(par%zmax/par%src_zscale(is))
     else
        photon%z = par%src_zscale(is)*rand_zexp(par%zmax/par%src_zscale(is))
     endif
  case ('sersic')
     !--- 3-D deprojected Sersic bulge: spherical radius from rand_sersic,
     !--- isotropic direction, oblate flattening by src_axial_ratio.  Reject
     !--- draws outside the box (cylinder rmax / zmax).
     rs_max = sqrt(par%zmax**2 + par%rmax**2)/par%src_reff(is)
     do
        rp   = par%src_reff(is) * rand_sersic(par%src_sersic_index(is), rs_max)
        cost = 2.0_wp*rand_number()-1.0_wp;  sint = sqrt(1.0_wp-cost*cost)
        phi  = twopi*rand_number()
        photon%x = rp*sint*cos(phi)
        photon%y = rp*sint*sin(phi)
        photon%z = rp*cost*par%src_axial_ratio(is)
        if (sqrt(photon%x**2+photon%y**2) <= par%rmax .and. &
            abs(photon%z) <= par%zmax) exit
     enddo
  case ('boxy', 'bar', 'xbar')
     !--- generalized-ellipsoid bulge (boxiness src_boxiness), scaled by
     !--- src_reff and flattened in z by src_axial_ratio.  Sampled in units of
     !--- the bulge scale (hard bound 8); rejected outside the box.
     do
        select case (trim(par%src_geometry(is)))
        case ('boxy')
           call rand_boxy(bx, by, bz, 0.1_wp, 8.0_wp, 8.0_wp, 8.0_wp, &
                          par%src_boxiness(is), par%src_sersic_index(is))
        case ('bar')
           call rand_bar(bx, by, bz, 0.1_wp, 8.0_wp, 8.0_wp, 8.0_wp, par%src_boxiness(is))
        case default   ! 'xbar'
           call rand_xbar(bx, by, bz, 8.0_wp, 8.0_wp, 8.0_wp, par%src_boxiness(is))
        end select
        photon%x = bx*par%src_reff(is)
        photon%y = by*par%src_reff(is)
        photon%z = bz*par%src_reff(is)*par%src_axial_ratio(is)
        if (sqrt(photon%x**2+photon%y**2) <= par%rmax .and. &
            abs(photon%z) <= par%zmax) exit
     enddo
  case default   ! 'point'
     photon%x = par%src_x(is);  photon%y = par%src_y(is);  photon%z = par%src_z(is)
  end select

  photon%icell = floor((photon%x-grid%xmin)/grid%dx)+1
  photon%jcell = floor((photon%y-grid%ymin)/grid%dy)+1
  photon%kcell = floor((photon%z-grid%zmin)/grid%dz)+1

  !--- isotropic direction.
  cost = 2.0_wp*rand_number()-1.0_wp;  sint = sqrt(1.0_wp-cost*cost)
  phi  = twopi*rand_number()
  photon%kx = sint*cos(phi);  photon%ky = sint*sin(phi);  photon%kz = cost

  if (par%use_sed) then
     !--- wavelength from this source's spectrum.
     photon%lambda  = sample_wavelength(src_spectrum(is), rand_number())
     photon%il      = sed_bin_of(photon%lambda)
     photon%s_ext   = sed_sext_at(photon%lambda)
     photon%albedo  = sed_albedo_at(photon%lambda)
     photon%hgg     = sed_hgg_at(photon%lambda)
  else
     !--- monochromatic: grey dust properties, no wavelength.  s_ext is unused
     !--- by the monochromatic raytrace but set to the safe unit value.
     photon%albedo  = par%albedo
     photon%hgg     = par%hgg
     photon%s_ext   = 1.0_wp
  endif

  photon%nscatt  = 0
  photon%inside  = .true.
  photon%wgt     = 1.0_wp
  photon%Lpacket = src_Lpacket(is)
  end subroutine gen_source_photon

  !---------------------------------------------------------------
  !--- quasi-random variant of gen_source_photon: the component (uq(2)),
  !--- direction (uq(4), uq(5)) and wavelength (uq(3), inverse CDF) come from the
  !--- scrambled Sobol point; only a 'point' component has no position draw, so
  !--- every other geometry keeps its Mersenne Twister position sampler
  !--- (rejection / stateful) exactly as in gen_source_photon.
  subroutine gen_source_photon_qmc(grid, photon, uq)
  implicit none
  type(grid_type),   intent(in)  :: grid
  type(photon_type), intent(inout) :: photon
  real(kind=wp),     intent(in)  :: uq(:)
  real(kind=wp) :: sint, cost, phi, rp, rs_max, tanp, bx, by, bz
  integer :: is

  !--- select source by luminosity CDF (quasi-random dimension 2).
  is = cdf_search(src_lum_cdf, uq(2))

  !--- position from this source's geometry (unchanged from gen_source_photon;
  !--- 'point' has no draw, the others keep Mersenne Twister sampling).
  select case (trim(par%src_geometry(is)))
  case ('uniform')
     rp   = rand_number()**(1.0_wp/3.0_wp) * par%rmax
     cost = 2.0_wp*rand_number()-1.0_wp;  sint = sqrt(1.0_wp-cost*cost)
     phi  = twopi*rand_number()
     photon%x = rp*sint*cos(phi);  photon%y = rp*sint*sin(phi);  photon%z = rp*cost
  case ('gaussian')
     photon%x = grid%xrange*rand_number()+grid%xmin
     photon%y = grid%yrange*rand_number()+grid%ymin
     photon%z = par%src_zscale(is)/sqrt(2.0_wp)*rand_gauss()
  case ('exponential', 'sech', 'exp_spiral')
     if (par%src_rscale(is) > 0.0_wp) then
        if (trim(par%src_geometry(is)) == 'exp_spiral' .and. par%spiral_m > 0) then
           tanp = tan(par%spiral_pitch*pi/180.0_wp)
           do
              rp  = par%src_rscale(is) * rand_r1exp(par%rmax/par%src_rscale(is))
              phi = twopi*rand_number()
              if (rp <= 0.0_wp) cycle
              if ((1.0_wp+par%spiral_amp)*rand_number() <= &
                  1.0_wp + par%spiral_amp*sin(par%spiral_m*(log(rp)/tanp - phi))) exit
           enddo
        else
           rp  = par%src_rscale(is) * rand_r1exp(par%rmax/par%src_rscale(is))
           phi = twopi*rand_number()
        endif
        photon%x = rp*cos(phi);  photon%y = rp*sin(phi)
     else
        photon%x = grid%xrange*rand_number()+grid%xmin
        photon%y = grid%yrange*rand_number()+grid%ymin
     endif
     if (trim(par%src_geometry(is)) == 'sech') then
        photon%z = par%src_zscale(is)*rand_sech2(par%zmax/par%src_zscale(is))
     else
        photon%z = par%src_zscale(is)*rand_zexp(par%zmax/par%src_zscale(is))
     endif
  case ('sersic')
     rs_max = sqrt(par%zmax**2 + par%rmax**2)/par%src_reff(is)
     do
        rp   = par%src_reff(is) * rand_sersic(par%src_sersic_index(is), rs_max)
        cost = 2.0_wp*rand_number()-1.0_wp;  sint = sqrt(1.0_wp-cost*cost)
        phi  = twopi*rand_number()
        photon%x = rp*sint*cos(phi)
        photon%y = rp*sint*sin(phi)
        photon%z = rp*cost*par%src_axial_ratio(is)
        if (sqrt(photon%x**2+photon%y**2) <= par%rmax .and. &
            abs(photon%z) <= par%zmax) exit
     enddo
  case ('boxy', 'bar', 'xbar')
     do
        select case (trim(par%src_geometry(is)))
        case ('boxy')
           call rand_boxy(bx, by, bz, 0.1_wp, 8.0_wp, 8.0_wp, 8.0_wp, &
                          par%src_boxiness(is), par%src_sersic_index(is))
        case ('bar')
           call rand_bar(bx, by, bz, 0.1_wp, 8.0_wp, 8.0_wp, 8.0_wp, par%src_boxiness(is))
        case default   ! 'xbar'
           call rand_xbar(bx, by, bz, 8.0_wp, 8.0_wp, 8.0_wp, par%src_boxiness(is))
        end select
        photon%x = bx*par%src_reff(is)
        photon%y = by*par%src_reff(is)
        photon%z = bz*par%src_reff(is)*par%src_axial_ratio(is)
        if (sqrt(photon%x**2+photon%y**2) <= par%rmax .and. &
            abs(photon%z) <= par%zmax) exit
     enddo
  case default   ! 'point'
     photon%x = par%src_x(is);  photon%y = par%src_y(is);  photon%z = par%src_z(is)
  end select

  photon%icell = floor((photon%x-grid%xmin)/grid%dx)+1
  photon%jcell = floor((photon%y-grid%ymin)/grid%dy)+1
  photon%kcell = floor((photon%z-grid%zmin)/grid%dz)+1

  !--- isotropic direction from the quasi-random dimensions 4 (mu) and 5 (phi).
  cost = 2.0_wp*uq(4)-1.0_wp;  sint = sqrt(1.0_wp-cost*cost)
  phi  = twopi*uq(5)
  photon%kx = sint*cos(phi);  photon%ky = sint*sin(phi);  photon%kz = cost

  if (par%use_sed) then
     !--- wavelength from this source's spectrum (inverse CDF, dimension 3).
     photon%lambda  = sample_wavelength(src_spectrum(is), uq(3))
     photon%il      = sed_bin_of(photon%lambda)
     photon%s_ext   = sed_sext_at(photon%lambda)
     photon%albedo  = sed_albedo_at(photon%lambda)
     photon%hgg     = sed_hgg_at(photon%lambda)
  else
     photon%albedo  = par%albedo
     photon%hgg     = par%hgg
     photon%s_ext   = 1.0_wp
  endif

  photon%nscatt  = 0
  photon%inside  = .true.
  photon%wgt     = 1.0_wp
  photon%Lpacket = src_Lpacket(is)
  end subroutine gen_source_photon_qmc

end module sources_mod
