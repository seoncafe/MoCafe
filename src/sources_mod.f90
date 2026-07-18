module sources_mod
!--- Multiple internal source components (monochromatic).  Each component has
!--- its own geometry, geometry parameters and luminosity; components are
!--- sampled in proportion to their luminosity, and every packet carries the
!--- same energy L_tot/nphotons, so the image normalization by par%luminosity
!--- (= sum of the component luminosities) is unchanged.
!---
!--- Activated only when par%nsource > 1; par%nsource <= 1 leaves the single
!--- source_geometry path of gen_photon untouched.  There is no wavelength or
!--- spectrum concept in this version: every packet uses the grey dust
!--- properties par%albedo and par%hgg.
  use define
  use random, only : rand_number, rand_gauss, rand_zexp
  implicit none
  private
  public :: setup_sources, gen_source_photon, gen_source_photon_qmc, use_sources

  logical :: use_sources = .false.
  integer :: nsrc = 0
  real(kind=wp), allocatable :: src_lum_cdf(:)   ! (nsrc) cumulative luminosity, normalized

contains
  !---------------------------------------------------------------
  subroutine setup_sources()
  use mpi
  implicit none
  real(kind=wp), allocatable :: lum(:)
  real(kind=wp) :: lsum
  integer :: is, ierr

  nsrc = par%nsource
  if (nsrc <= 1) then
     use_sources = .false.;  return
  endif
  if (nsrc > MAX_SRC) then
     if (mpar%p_rank == 0) write(*,'(a,i0,a,i0,a)') &
        'ERROR: par%nsource = ', nsrc, ' exceeds MAX_SRC = ', MAX_SRC, '.'
     call MPI_FINALIZE(ierr);  stop
  endif

  !--- gen_source_photon does not set the Stokes vector or its reference frame,
  !--- so multiple internal sources cannot be combined with the Stokes path.
  if (par%use_stokes) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: multiple internal sources (par%nsource > 1) do not support par%use_stokes.'
     call MPI_FINALIZE(ierr);  stop
  endif

  !--- only the geometries that the single-source path of this monochromatic
  !--- version can sample are accepted.
  do is = 1, nsrc
     select case (trim(par%src_geometry(is)))
     case ('point', 'uniform', 'uniform_xy', 'gaussian', 'exponential')
        !--- supported
     case default
        if (mpar%p_rank == 0) write(*,'(a,i0,3a)') 'ERROR: par%src_geometry(', is, ') = ''', &
           trim(par%src_geometry(is)), &
           ''' (use ''point'', ''uniform'', ''uniform_xy'', ''gaussian'' or ''exponential'').'
        call MPI_FINALIZE(ierr);  stop
     end select
  enddo

  allocate(src_lum_cdf(nsrc), lum(nsrc))

  !--- component luminosity: par%src_lum when set, else the equal split of
  !--- par%luminosity.
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
        'ERROR: multiple internal sources have zero total luminosity.'
     call MPI_FINALIZE(ierr);  stop
  endif
  par%luminosity = lsum          ! total luminosity is the sum of the components

  !--- luminosity CDF over the components.  Each component emits
  !--- nphotons*lum(is)/lsum packets, so every packet carries lsum/nphotons.
  src_lum_cdf(1) = lum(1)
  do is = 2, nsrc
     src_lum_cdf(is) = src_lum_cdf(is-1) + lum(is)
  enddo
  src_lum_cdf(:) = src_lum_cdf(:)/lsum

  use_sources = .true.

  if (mpar%p_rank == 0) then
     write(*,'(a)')    '--- Multiple internal sources ---'
     write(*,'(a,i0)') 'N sources                 : ', nsrc
     do is = 1, nsrc
        write(*,'(a,i2,3a,es11.3)') '  src ', is, '  geom=', &
           trim(par%src_geometry(is)), '  lum=', lum(is)
     enddo
     write(*,'(a,es12.4)') 'total luminosity          : ', lsum
  endif

  deallocate(lum)
  end subroutine setup_sources

  !---------------------------------------------------------------
  !--- index of the first CDF entry >= u (binary search).
  function cdf_search(cdf, u) result(idx)
  implicit none
  real(kind=wp), intent(in) :: cdf(:), u
  integer :: idx, lo, hi, mid
  lo = 1;  hi = size(cdf)
  do while (lo < hi)
     mid = (lo+hi)/2
     if (u <= cdf(mid)) then;  hi = mid;  else;  lo = mid+1;  endif
  enddo
  idx = lo
  end function cdf_search

  !---------------------------------------------------------------
  !--- sample the position of component is into the photon.
  subroutine sample_source_position(grid, photon, is)
  implicit none
  type(grid_type),   intent(in)    :: grid
  type(photon_type), intent(inout) :: photon
  integer,           intent(in)    :: is
  real(kind=wp) :: sint, cost, phi, rp

  select case (trim(par%src_geometry(is)))
  case ('uniform')
     if (par%rmax > 0.0_wp) then
        rp   = rand_number()**(1.0_wp/3.0_wp) * par%rmax
        cost = 2.0_wp*rand_number()-1.0_wp
        sint = sqrt(1.0_wp-cost*cost)
        phi  = twopi*rand_number()
        photon%x = rp*sint*cos(phi)
        photon%y = rp*sint*sin(phi)
        photon%z = rp*cost
     else
        photon%x = (2.0_wp*rand_number()-1.0_wp)*grid%xmax
        photon%y = (2.0_wp*rand_number()-1.0_wp)*grid%ymax
        photon%z = (2.0_wp*rand_number()-1.0_wp)*grid%zmax
     endif
  case ('uniform_xy')
     photon%x = grid%xrange*rand_number()+grid%xmin
     photon%y = grid%yrange*rand_number()+grid%ymin
     photon%z = 0.0_wp
  case ('gaussian')
     photon%x = grid%xrange*rand_number()+grid%xmin
     photon%y = grid%yrange*rand_number()+grid%ymin
     photon%z = par%src_zscale(is)/sqrt(2.0_wp)*rand_gauss()
  case ('exponential')
     photon%x = grid%xrange*rand_number()+grid%xmin
     photon%y = grid%yrange*rand_number()+grid%ymin
     photon%z = par%src_zscale(is)*rand_zexp(par%zmax/par%src_zscale(is))
  case default   ! 'point'
     photon%x = par%src_x(is);  photon%y = par%src_y(is);  photon%z = par%src_z(is)
  end select

  photon%icell = floor((photon%x-grid%xmin)/grid%dx)+1
  photon%jcell = floor((photon%y-grid%ymin)/grid%dy)+1
  photon%kcell = floor((photon%z-grid%zmin)/grid%dz)+1
  end subroutine sample_source_position

  !---------------------------------------------------------------
  !--- generate a photon from a luminosity-weighted source component: pick the
  !--- component, sample its position and an isotropic direction, and set the
  !--- grey dust properties.
  subroutine gen_source_photon(grid, photon)
  implicit none
  type(grid_type),   intent(in)    :: grid
  type(photon_type), intent(inout) :: photon
  real(kind=wp) :: sint, cost, phi
  integer :: is

  is = cdf_search(src_lum_cdf, rand_number())
  call sample_source_position(grid, photon, is)

  !--- isotropic direction.
  cost = 2.0_wp*rand_number()-1.0_wp;  sint = sqrt(1.0_wp-cost*cost)
  phi  = twopi*rand_number()
  photon%kx = sint*cos(phi);  photon%ky = sint*sin(phi);  photon%kz = cost

  photon%albedo = par%albedo
  photon%hgg    = par%hgg
  photon%nscatt = 0
  photon%inside = .true.
  photon%wgt    = 1.0_wp
  end subroutine gen_source_photon

  !---------------------------------------------------------------
  !--- quasi-random variant: the component (uq(2)) and the direction (uq(4),
  !--- uq(5)) come from the scrambled Sobol point.  Only a 'point' component has
  !--- no position draw, so every other geometry keeps its Mersenne Twister
  !--- position sampler exactly as in gen_source_photon.
  subroutine gen_source_photon_qmc(grid, photon, uq)
  implicit none
  type(grid_type),   intent(in)    :: grid
  type(photon_type), intent(inout) :: photon
  real(kind=wp),     intent(in)    :: uq(:)
  real(kind=wp) :: sint, cost, phi
  integer :: is

  is = cdf_search(src_lum_cdf, uq(2))
  call sample_source_position(grid, photon, is)

  !--- isotropic direction from the quasi-random dimensions 4 (mu) and 5 (phi).
  cost = 2.0_wp*uq(4)-1.0_wp;  sint = sqrt(1.0_wp-cost*cost)
  phi  = twopi*uq(5)
  photon%kx = sint*cos(phi);  photon%ky = sint*sin(phi);  photon%kz = cost

  photon%albedo = par%albedo
  photon%hgg    = par%hgg
  photon%nscatt = 0
  photon%inside = .true.
  photon%wgt    = 1.0_wp
  end subroutine gen_source_photon_qmc

end module sources_mod
