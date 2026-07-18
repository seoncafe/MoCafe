module sources_mod
!--- Multiple stellar source components (MoCafe v2.00, Stage 6).  Each
!--- component has its own spectrum (Planck or 2-column file), luminosity,
!--- geometry and geometry parameters.  Sources are sampled in proportion to
!--- their luminosity; each carries its own wavelength-bin alias table so a
!--- run mixes, e.g., a hot young population and a cool old population with
!--- different spatial distributions.  When par%nsource == 1 this reduces to
!--- the single-source SED path (sed_mod) and is not activated.
  use define
  use random,  only : random_alias_setup, rand_alias_choise, rand_number, rand_gauss, &
                      rand_zexp, rand_sech2, rand_r1exp
  use sed_mod, only : sed_nlam, sed_wave, sed_dwave, sed_sext, sed_albedo, sed_hgg, &
                      planck_shape, read_spectrum_file, interp_clamped, convert_spectrum_units, &
                      cdf_search
  use random_bulge,  only : rand_sersic, rand_boxy, rand_bar, rand_xbar
  implicit none
  private
  public :: setup_sources, gen_source_photon, gen_source_photon_qmc, use_sources

  logical :: use_sources = .false.
  integer :: nsrc = 0
  real(kind=wp), allocatable :: src_lum_cdf(:)          ! (nsrc) cumulative luminosity, normalized
  real(kind=wp), allocatable :: src_pdf(:,:)            ! (nlam, nsrc) spectrum alias PDF for each source
  integer,       allocatable :: src_alias(:,:)          ! (nlam, nsrc)
  real(kind=wp), allocatable :: src_cdf(:,:)            ! (nlam, nsrc) cumulative spectrum (inverse-CDF path)
  real(kind=wp), allocatable :: src_lumfrac(:,:)        ! (nlam, nsrc) luminosity fraction per bin (for output)
  real(kind=wp), allocatable :: src_Lpacket(:)          ! (nsrc) energy per packet for this source

contains
  !---------------------------------------------------------------
  subroutine setup_sources()
  use mpi
  implicit none
  real(kind=wp), allocatable :: lum(:), lam_f(:), lum_f(:), pdf(:)
  real(kind=wp) :: lsum, lnorm
  integer :: is, il, nsp, ierr
  logical :: is_phys, absolute_is
  logical, allocatable :: derived(:)

  nsrc = par%nsource
  if (nsrc <= 1) then
     use_sources = .false.;  return
  endif

  !--- monochromatic (non-SED) multiple internal sources: luminosity-weighted
  !--- component selection + geometry-based position sampling only, with no
  !--- wavelength/spectrum concept.  The SED-only arrays (src_pdf, src_alias,
  !--- src_lumfrac) are left unallocated and never referenced on this path.
  if (.not. par%use_sed) then
     call setup_sources_mono()
     return
  endif

  use_sources = .true.
  is_phys = trim(par%spectrum_type) /= 'shape'

  allocate(src_lum_cdf(nsrc), src_pdf(sed_nlam,nsrc), src_alias(sed_nlam,nsrc), src_cdf(sed_nlam,nsrc))
  allocate(src_lumfrac(sed_nlam,nsrc), src_Lpacket(nsrc), lum(nsrc), pdf(sed_nlam))
  allocate(derived(nsrc));  derived(:) = .false.

  !--- spectrum of each source -> luminosity fraction per bin + alias table, and
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
        do il = 1, sed_nlam
           if (sed_wave(il) < lam_f(1) .or. sed_wave(il) > lam_f(nsp)) then
              pdf(il) = 0.0_wp
           else
              pdf(il) = interp_clamped(log(lam_f), lum_f, log(sed_wave(il)))*sed_dwave(il)
           endif
        enddo
        deallocate(lam_f, lum_f)
     else if (par%src_tstar(is) > 0.0_wp) then
        do il = 1, sed_nlam
           pdf(il) = planck_shape(sed_wave(il), par%src_tstar(is))*sed_dwave(il)
        enddo
     else
        if (mpar%p_rank == 0) write(*,'(a,i0,a)') &
           'ERROR: source ', is, ' needs src_tstar or src_spectrum.'
        call MPI_FINALIZE(ierr);  stop
     endif
     lnorm = sum(pdf)
     if (.not. (lnorm > 0.0_wp)) then
        if (mpar%p_rank == 0) write(*,'(a,i0,a)') 'ERROR: source ', is, ' has zero luminosity on the grid.'
        call MPI_FINALIZE(ierr);  stop
     endif

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

     src_lumfrac(:,is) = pdf(:)/lnorm
     src_pdf(:,is)     = src_lumfrac(:,is)
     call random_alias_setup(src_pdf(:,is), src_alias(:,is))
     !--- cumulative spectrum for the quasi-random inverse-CDF wavelength
     !--- sampler; the alias table above stays the default sampling path.
     src_cdf(1,is) = src_lumfrac(1,is)
     do il = 2, sed_nlam
        src_cdf(il,is) = src_cdf(il-1,is) + src_lumfrac(il,is)
     enddo
     src_cdf(sed_nlam,is) = 1.0_wp
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

  deallocate(lum, pdf, derived)
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
  integer :: is, il, lo, hi, mid

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
     il = rand_alias_choise(src_pdf(:,is), src_alias(:,is))
     photon%il      = il
     photon%lambda  = sed_wave(il)
     photon%s_ext   = sed_sext(il)
     photon%albedo  = sed_albedo(il)
     photon%hgg     = sed_hgg(il)
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
  integer :: is, il

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
     il = cdf_search(src_cdf(:,is), uq(3))
     photon%il      = il
     photon%lambda  = sed_wave(il)
     photon%s_ext   = sed_sext(il)
     photon%albedo  = sed_albedo(il)
     photon%hgg     = sed_hgg(il)
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
