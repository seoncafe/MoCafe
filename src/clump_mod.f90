module clump_mod
!---------------------------------------------------------------------------
! Clumpy dust medium for MoCafe (par%grid_type = 'clump').
!
! Geometry: N_clumps spherical dust clumps of radius cl_radius(:)
! placed uniformly at random (non-overlapping via RSA) inside a sphere of
! radius sphere_R.  Each clump carries a grey dust opacity cl_rhokap(:)
! (extinction per unit code length).  Outside the clumps the medium is vacuum.
! When all radial profiles are 'constant', every clump entry equals the
! corresponding reference scalar; otherwise the arrays are populated from
! radial profiles.
!
! This is the dust-only slim of LaRT_v2.00/clump_mod.f90: all Lyman-alpha /
! velocity / temperature / Voigt / frequency handling is stripped (see
! AMR_CLUMPS_PLAN.md Part B).  The clump subsystem is purely geometric plus a
! scalar opacity of each clump.
!
! Algorithms:
!   Placement : linked-list RSA, O(N_cl) expected time for f_vol < 35%.
!   Raytrace  : DDA through CSR acceleration grid (Amanatides & Woo 1987).
!   Memory    : MPI-3 shared memory (one copy per node).
!---------------------------------------------------------------------------
  use define
  use memory_mod
  use random
  use mpi
  implicit none
  public

  !--- Physical properties of the clump population
  integer(int64), save :: N_clumps    = 0_int64
  real(kind=wp),  save :: sphere_R    = 0.0_wp   ! outer sphere radius [code units]
  real(kind=wp),  save :: r_min_clump = 0.0_wp   ! inner placement radius [code units]
                                                 ! (= max(0, par%rmin); 0 -> filled sphere)

  !--- Dust properties of each clump (MPI shared memory, dimension N_clumps).
  real(kind=wp), pointer, save :: cl_radius(:)  => null()  ! clump radius [code units]
  real(kind=wp), pointer, save :: cl_radius2(:) => null()  ! cl_radius(icl)**2
  real(kind=wp), pointer, save :: cl_rhokap(:)  => null()  ! grey dust opacity per code unit

  !--- Reference / representative scalars used during initialisation, grid
  !    setup, and FITS-header reporting.  cl_radius_max also sizes the RSA /
  !    CSR acceleration grids so the 27-neighbor overlap search stays complete
  !    when r_cl varies between clumps.
  real(kind=wp), save :: cl_radius_max = 0.0_wp
  real(kind=wp), save :: cl_rhokap_ref = 0.0_wp

  !--- Clump positions [code units] - MPI shared memory.
  real(kind=dp), pointer, save :: cl_x(:) => null()
  real(kind=dp), pointer, save :: cl_y(:) => null()
  real(kind=dp), pointer, save :: cl_z(:) => null()

  !--- Radial-profile handling (geometry only).
  real(kind=wp), save :: base_radius_in  = 0.0_wp     ! par%clump_radius (peak radius)
  real(kind=wp), save :: base_rhokap_in  = 0.0_wp     ! peak rhokap that shapes scale
  logical,       save :: profiles_active = .false.    ! .true. if any profile != 'constant'
  logical,       save :: clumps_from_file = .false.   ! .true. if loaded from clump_input_file
  logical,       save :: has_overlap     = .false.    ! .true. if loaded clumps overlap

  !--- Tabulated radial CDF for inverse-CDF sampling of clump positions.
  integer, parameter :: NPROF = 4001
  real(kind=wp), save :: prof_dr = 0.0_wp
  real(kind=wp), save :: prof_r(NPROF)
  real(kind=wp), save :: prof_shape_number(NPROF)     ! n_cl(r) shape (unnormalized)
  real(kind=wp), save :: prof_shape_radius(NPROF)     ! r_cl(r) shape
  real(kind=wp), save :: prof_shape_density(NPROF)    ! rhokap(r) shape
  real(kind=wp), save :: prof_cdf_pos(NPROF)          ! CDF for r ~ shape_number * r^2

  !--- Optional tabulated profile from clump_profile_file.
  integer, save :: NTAB = 0
  real(kind=wp), allocatable, save :: tab_r(:)
  real(kind=wp), allocatable, save :: tab_radius(:)
  real(kind=wp), allocatable, save :: tab_density(:)
  real(kind=wp), allocatable, save :: tab_number(:)

  !--- CSR acceleration grid (MPI shared memory)
  !    Cell (i,j,k), 0-based -> 1-based index = 1 + i + j*cgx + k*cgx*cgy
  !    cg_start(icell) .. cg_start(icell+1)-1 -> entries in cg_list for that cell
  integer,       save :: cgx = 0, cgy = 0, cgz = 0
  real(kind=wp), save :: cg_xmin, cg_ymin, cg_zmin
  real(kind=wp), save :: cg_dx, cg_dy, cg_dz
  real(kind=wp), save :: cg_inv_dx, cg_inv_dy, cg_inv_dz
  integer(int32), pointer, save :: cg_start(:) => null()   ! size ncells+1
  integer(int32), pointer, save :: cg_list(:)  => null()   ! size total_registrations

contains

  !===========================================================================
  pure integer function cg_cell_idx(i, j, k)
  integer, intent(in) :: i, j, k
  cg_cell_idx = 1 + i + j*cgx + k*cgx*cgy
  end function cg_cell_idx
  !===========================================================================

  !===========================================================================
  ! Radial-profile evaluator. Returns the multiplicative shape factor for a
  ! profile name at radius r (code units, 0 <= r <= R_box).
  !   'constant' : 1;  'powerlaw' : (max(r,0.05 r0)/max(r0,0.05 r0))^(-alpha);
  !   'gaussian' : exp(-(r/r0)^2);  'exponential' : exp(-r/r0);
  !   'file' : linear interpolation from the tabulated column.
  ! axis_id selects the file column: 1=radius, 2=density, 3=number-shape.
  !===========================================================================
  pure real(kind=wp) function profile_shape(shape_name, alpha, r0, r, axis_id) &
       result(f)
  character(len=*), intent(in) :: shape_name
  real(kind=wp),    intent(in) :: alpha, r0, r
  integer,          intent(in) :: axis_id
  real(kind=wp) :: r_eff, r0_eff, r_floor

  select case (trim(shape_name))
  case ('constant', '')
     f = 1.0_wp
  case ('powerlaw', 'power_law')
     if (r0 > 0.0_wp) then
        r_floor = 0.05_wp * r0
        r_eff   = max(r,  r_floor)
        r0_eff  = max(r0, r_floor)
        f = (r_eff / r0_eff) ** (-alpha)
     else
        f = 1.0_wp
     end if
  case ('gaussian')
     if (r0 > 0.0_wp) then
        f = exp(-(r/r0)**2)
     else
        f = 1.0_wp
     end if
  case ('exponential')
     if (r0 > 0.0_wp) then
        f = exp(-r/r0)
     else
        f = 1.0_wp
     end if
  case ('file')
     f = profile_file_interp(r, axis_id)
  case default
     f = 1.0_wp
  end select
  end function profile_shape
  !===========================================================================

  !===========================================================================
  pure real(kind=wp) function profile_file_interp(r, axis_id) result(f)
  real(kind=wp), intent(in) :: r
  integer,       intent(in) :: axis_id
  integer :: lo, hi, mid
  real(kind=wp) :: t

  if (NTAB <= 1) then
     f = 1.0_wp;  return
  end if

  if (r <= tab_r(1)) then
     lo = 1;  hi = 2
  else if (r >= tab_r(NTAB)) then
     lo = NTAB-1;  hi = NTAB
  else
     lo = 1;  hi = NTAB
     do while (hi - lo > 1)
        mid = (lo + hi) / 2
        if (r >= tab_r(mid)) then
           lo = mid
        else
           hi = mid
        end if
     end do
  end if

  t = (r - tab_r(lo)) / max(tab_r(hi) - tab_r(lo), tiny(1.0_wp))
  t = max(0.0_wp, min(1.0_wp, t))

  select case (axis_id)
  case (1)
     if (allocated(tab_radius)) then
        f = (1.0_wp - t) * tab_radius(lo) + t * tab_radius(hi)
     else
        f = 1.0_wp
     end if
  case (2)
     if (allocated(tab_density)) then
        f = (1.0_wp - t) * tab_density(lo) + t * tab_density(hi)
     else
        f = 1.0_wp
     end if
  case (3)
     if (allocated(tab_number)) then
        f = (1.0_wp - t) * tab_number(lo) + t * tab_number(hi)
     else
        f = 1.0_wp
     end if
  case default
     f = 1.0_wp
  end select
  end function profile_file_interp
  !===========================================================================

  !===========================================================================
  pure real(kind=wp) function shape_radius(r) result(f)
  real(kind=wp), intent(in) :: r
  real(kind=wp) :: r0_use
  r0_use = par%clump_radius_r0
  if (r0_use <= 0.0_wp) r0_use = sphere_R
  f = profile_shape(par%clump_radius_profile,  par%clump_radius_alpha, r0_use, r, 1)
  end function shape_radius

  pure real(kind=wp) function shape_density(r) result(f)
  real(kind=wp), intent(in) :: r
  real(kind=wp) :: r0_use
  r0_use = par%clump_density_r0
  if (r0_use <= 0.0_wp) r0_use = sphere_R
  f = profile_shape(par%clump_density_profile, par%clump_density_alpha, r0_use, r, 2)
  end function shape_density

  pure real(kind=wp) function shape_number(r) result(f)
  real(kind=wp), intent(in) :: r
  real(kind=wp) :: r0_use
  r0_use = par%clump_number_r0
  if (r0_use <= 0.0_wp) r0_use = sphere_R
  f = profile_shape(par%clump_number_profile,  par%clump_number_alpha, r0_use, r, 3)
  end function shape_number
  !===========================================================================

  !===========================================================================
  ! Build the radial CDF for inverse-CDF sampling of clump positions.
  ! P(r) ~ shape_number(r) * r^2  on [0, sphere_R].  Also sets cl_radius_max.
  !===========================================================================
  subroutine build_radial_profile_tables(R_box)
  real(kind=wp), intent(in) :: R_box
  integer        :: i
  real(kind=wp)  :: r, integrand, last_integrand, total, rcl_local

  prof_dr = R_box / real(NPROF - 1, wp)
  cl_radius_max = base_radius_in

  prof_r(1)             = 0.0_wp
  prof_shape_number(1)  = shape_number(0.0_wp)
  prof_shape_radius(1)  = shape_radius(0.0_wp)
  prof_shape_density(1) = shape_density(0.0_wp)
  prof_cdf_pos(1)       = 0.0_wp

  if (prof_r(1) < r_min_clump) prof_shape_number(1) = 0.0_wp

  if (prof_shape_number(1) > 0.0_wp) then
     rcl_local = base_radius_in * prof_shape_radius(1)
     if (rcl_local > cl_radius_max) cl_radius_max = rcl_local
  end if

  last_integrand = prof_shape_number(1) * 0.0_wp

  do i = 2, NPROF
     r = real(i-1, wp) * prof_dr
     prof_r(i)             = r
     prof_shape_number(i)  = shape_number(r)
     prof_shape_radius(i)  = shape_radius(r)
     prof_shape_density(i) = shape_density(r)
     if (r < r_min_clump) prof_shape_number(i) = 0.0_wp
     integrand             = prof_shape_number(i) * r * r
     prof_cdf_pos(i) = prof_cdf_pos(i-1) + 0.5_wp * (last_integrand + integrand) * prof_dr
     last_integrand = integrand

     if (prof_shape_number(i) > 0.0_wp) then
        rcl_local = base_radius_in * prof_shape_radius(i)
        if (rcl_local > cl_radius_max) cl_radius_max = rcl_local
     end if
  end do

  total = prof_cdf_pos(NPROF)
  if (total > 0.0_wp) then
     prof_cdf_pos(:) = prof_cdf_pos(:) / total
  else
     do i = 1, NPROF
        r = real(i-1, wp) * prof_dr
        prof_cdf_pos(i) = (r / R_box) ** 3
     end do
  end if
  end subroutine build_radial_profile_tables
  !===========================================================================

  !===========================================================================
  real(kind=wp) function sample_clump_radius() result(r_out)
  real(kind=wp) :: u, t
  integer       :: lo, hi, mid

  u = rand_number()
  if (u <= prof_cdf_pos(1)) then
     r_out = prof_r(1);  return
  end if
  if (u >= prof_cdf_pos(NPROF)) then
     r_out = prof_r(NPROF);  return
  end if

  lo = 1;  hi = NPROF
  do while (hi - lo > 1)
     mid = (lo + hi) / 2
     if (u >= prof_cdf_pos(mid)) then
        lo = mid
     else
        hi = mid
     end if
  end do

  t = (u - prof_cdf_pos(lo)) / max(prof_cdf_pos(hi) - prof_cdf_pos(lo), tiny(1.0_wp))
  t = max(0.0_wp, min(1.0_wp, t))
  r_out = (1.0_wp - t) * prof_r(lo) + t * prof_r(hi)
  end function sample_clump_radius
  !===========================================================================

  !===========================================================================
  pure real(kind=wp) function integrate_table(integrand) result(total)
  real(kind=wp), intent(in) :: integrand(NPROF)
  integer :: i
  total = 0.0_wp
  do i = 2, NPROF
     total = total + 0.5_wp * (integrand(i-1) + integrand(i)) * prof_dr
  end do
  end function integrate_table
  !===========================================================================

  !===========================================================================
  ! f_cov = integral_0^R A_norm * shape_number(r) * pi * r_cl(r)^2 dr
  !===========================================================================
  pure real(kind=wp) function f_cov_LOS_quad(A_norm) result(fcov)
  real(kind=wp), intent(in) :: A_norm
  real(kind=wp) :: rcl, integrand(NPROF)
  integer :: i
  do i = 1, NPROF
     rcl = base_radius_in * prof_shape_radius(i)
     integrand(i) = A_norm * prof_shape_number(i) * pi * rcl * rcl
  end do
  fcov = integrate_table(integrand)
  end function f_cov_LOS_quad
  !===========================================================================

  !===========================================================================
  ! f_vol = (4 pi / V_shell) * integral_0^R A_norm * shape_number(r)*r_cl^3*r^2 dr
  !===========================================================================
  pure real(kind=wp) function f_vol_quad(A_norm, R_box) result(fvol)
  real(kind=wp), intent(in) :: A_norm, R_box
  real(kind=wp) :: rcl, r, integrand(NPROF), V_shell
  integer :: i
  do i = 1, NPROF
     r   = prof_r(i)
     rcl = base_radius_in * prof_shape_radius(i)
     integrand(i) = A_norm * prof_shape_number(i) * rcl**3 * r * r
  end do
  V_shell = max(R_box**3 - r_min_clump**3, tiny(1.0_wp))
  fvol = fourpi * integrate_table(integrand) / V_shell
  end function f_vol_quad
  !===========================================================================

  !===========================================================================
  ! N = integral_0^R A_norm * shape_number(r) * 4 pi r^2 dr
  !===========================================================================
  pure real(kind=wp) function total_count_quad(A_norm) result(Ntot)
  real(kind=wp), intent(in) :: A_norm
  real(kind=wp) :: integrand(NPROF), r
  integer :: i
  do i = 1, NPROF
     r = prof_r(i)
     integrand(i) = A_norm * prof_shape_number(i) * fourpi * r * r
  end do
  Ntot = integrate_table(integrand)
  end function total_count_quad
  !===========================================================================

  !===========================================================================
  ! Geometric factor for the radial-sightline taumax integral through the
  ! clumpy medium (dust analog of LaRT's LOS_geometric_quad with the Voigt
  ! factor removed):
  !   taumax = GF(A_norm) * base_rhokap   (profile mode)
  !===========================================================================
  pure real(kind=wp) function LOS_geometric_quad(A_norm) result(GF)
  real(kind=wp), intent(in) :: A_norm
  real(kind=wp) :: integrand(NPROF)
  integer :: i
  do i = 1, NPROF
     integrand(i) = prof_shape_number(i) * prof_shape_radius(i)**3 * prof_shape_density(i)
  end do
  GF = (fourpi / 3.0_wp) * A_norm * base_radius_in**3 * integrate_table(integrand)
  end function LOS_geometric_quad
  !===========================================================================

  !===========================================================================
  ! Load tabulated radial profile from par%clump_profile_file.
  ! ASCII format (whitespace-separated, '#' = comment), 5 columns:
  !   r [code units]   r_cl_factor   rhokap_factor   T [K, ignored]   n_cl_shape
  ! Column 4 (temperature) is read for LaRT file-format compatibility but
  ! ignored in the dust-only port.
  !===========================================================================
  subroutine load_clump_profile_file(R_box)
  real(kind=wp), intent(in) :: R_box
  integer :: u, ios, nread, count_lines, ierr
  character(len=512) :: linebuf
  real(kind=wp) :: r_v, rc_v, nH_v, T_v, nc_v

  if (len_trim(par%clump_profile_file) == 0) return

  open(newunit=u, file=trim(par%clump_profile_file), status='old', &
       action='read', iostat=ios)
  if (ios /= 0) then
     write(*,*) 'ERROR: cannot open clump_profile_file: ', trim(par%clump_profile_file)
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if

  count_lines = 0
  do
     read(u, '(a)', iostat=ios) linebuf
     if (ios /= 0) exit
     linebuf = adjustl(linebuf)
     if (len_trim(linebuf) == 0) cycle
     if (linebuf(1:1) == '#') cycle
     count_lines = count_lines + 1
  end do
  rewind(u)

  if (count_lines < 2) then
     write(*,*) 'ERROR: clump_profile_file needs at least 2 data rows'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if

  NTAB = count_lines
  if (allocated(tab_r))       deallocate(tab_r)
  if (allocated(tab_radius))  deallocate(tab_radius)
  if (allocated(tab_density)) deallocate(tab_density)
  if (allocated(tab_number))  deallocate(tab_number)
  allocate(tab_r(NTAB), tab_radius(NTAB), tab_density(NTAB), tab_number(NTAB))

  nread = 0
  do
     read(u, '(a)', iostat=ios) linebuf
     if (ios /= 0) exit
     linebuf = adjustl(linebuf)
     if (len_trim(linebuf) == 0) cycle
     if (linebuf(1:1) == '#') cycle
     read(linebuf, *, iostat=ios) r_v, rc_v, nH_v, T_v, nc_v
     if (ios /= 0) cycle
     nread = nread + 1
     if (nread > NTAB) exit
     tab_r(nread)      = r_v
     tab_radius(nread) = rc_v
     tab_density(nread) = nH_v
     tab_number(nread) = nc_v
  end do
  close(u)

  if (nread /= NTAB) NTAB = nread

  if (mpar%p_rank == 0) then
     write(*,'(a,i6,a)') ' Clumps: loaded ', NTAB, ' rows from clump_profile_file'
     write(*,'(a,2es12.4)') ' Clumps: profile r range = ', tab_r(1), tab_r(NTAB)
  end if
  end subroutine load_clump_profile_file
  !===========================================================================

  !===========================================================================
  subroutine init_clumps(R_sphere)
  !---------------------------------------------------------------------------
  ! Derive N_clumps, compute the dust opacity, allocate MPI shared memory,
  ! run RSA placement, and build the CSR acceleration grid.
  !---------------------------------------------------------------------------
  implicit none
  real(kind=wp), intent(in) :: R_sphere
  real(kind=wp) :: A_norm, fvol_unit, fcov_unit, fvol_realized, fcov_realized
  real(kind=wp) :: GF_los, shell_R2_factor
  integer       :: ierr

  sphere_R       = R_sphere
  base_radius_in = par%clump_radius
  if (base_radius_in <= 0.0_wp) then
     if (mpar%p_rank == 0) write(*,*) 'ERROR: clump_radius must be > 0'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if

  r_min_clump = max(0.0_wp, par%rmin)
  if (r_min_clump >= R_sphere) then
     if (mpar%p_rank == 0) write(*,'(a,2es12.4)') &
        'ERROR: par%rmin must be < par%rmax for clump placement; rmin, rmax = ', &
        r_min_clump, R_sphere
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if

  !--- pre-built clump file: load and skip internal generation.
  if (len_trim(par%clump_input_file) > 0) then
     call read_clumps_info(trim(par%clump_input_file), R_sphere)
     call check_has_overlap()
     return
  end if

  profiles_active = (trim(par%clump_radius_profile)  /= 'constant') .or. &
                    (trim(par%clump_density_profile) /= 'constant') .or. &
                    (trim(par%clump_number_profile)  /= 'constant')

  if (profiles_active .and. ( &
       trim(par%clump_radius_profile)  == 'file' .or. &
       trim(par%clump_density_profile) == 'file' .or. &
       trim(par%clump_number_profile)  == 'file')) then
     call load_clump_profile_file(R_sphere)
  end if

  if (profiles_active) then
     call build_radial_profile_tables(R_sphere)
  else
     cl_radius_max = base_radius_in
  end if

  !--- derive N_clumps + realized f_vol/f_cov.
  A_norm = 0.0_wp
  if (.not. profiles_active) then
     if (par%clump_N_clumps > 0) then
        N_clumps = int(par%clump_N_clumps, int64)
     else if (par%clump_f_vol > 0.0_wp) then
        N_clumps = nint(par%clump_f_vol * (R_sphere**3 - r_min_clump**3) / cl_radius_max**3, int64)
     else if (par%clump_f_cov > 0.0_wp) then
        N_clumps = nint((4.0_wp/3.0_wp)*par%clump_f_cov * &
                        (R_sphere**2 + R_sphere*r_min_clump + r_min_clump**2) / &
                        cl_radius_max**2, int64)
     else
        if (mpar%p_rank == 0) &
           write(*,*) 'ERROR: specify clump_N_clumps, clump_f_vol, or clump_f_cov'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     fvol_realized = real(N_clumps, wp) * cl_radius_max**3 / &
                     max(R_sphere**3 - r_min_clump**3, tiny(1.0_wp))
     fcov_realized = 0.75_wp * real(N_clumps, wp) * cl_radius_max**2 / &
                     max(R_sphere**2 + R_sphere*r_min_clump + r_min_clump**2, tiny(1.0_wp))
  else
     if (par%clump_N_clumps > 0) then
        N_clumps = int(par%clump_N_clumps, int64)
        A_norm   = real(N_clumps, wp) / max(total_count_quad(1.0_wp), tiny(1.0_wp))
     else if (par%clump_f_vol > 0.0_wp) then
        fvol_unit = f_vol_quad(1.0_wp, R_sphere)
        A_norm    = par%clump_f_vol / max(fvol_unit, tiny(1.0_wp))
        N_clumps  = nint(total_count_quad(A_norm), int64)
     else if (par%clump_f_cov > 0.0_wp) then
        fcov_unit = f_cov_LOS_quad(1.0_wp)
        A_norm    = par%clump_f_cov / max(fcov_unit, tiny(1.0_wp))
        N_clumps  = nint(total_count_quad(A_norm), int64)
     else
        if (mpar%p_rank == 0) &
           write(*,*) 'ERROR: specify clump_N_clumps, clump_f_vol, or clump_f_cov'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     fvol_realized = f_vol_quad(A_norm, R_sphere)
     fcov_realized = f_cov_LOS_quad(A_norm)
  end if
  if (N_clumps <= 0_int64) N_clumps = 1_int64

  !--- clump dust opacity (peak value).  Direct inputs for each clump take
  !    priority over the system-level back-solve from par%taumax / par%tauhomo.
  shell_R2_factor = R_sphere**2 + R_sphere*r_min_clump + r_min_clump**2
  if (par%clump_tau0 > 0.0_wp) then
     ! dimensionless dust optical depth of one clump, center -> surface
     cl_rhokap_ref = par%clump_tau0 / base_radius_in
  else if (par%clump_ndust > 0.0_wp) then
     ! dust density [code-unit^-3 in physical cm^-3]; MoCafe convention
     ! rhokap = density * cext_dust * distance2cm  (cf. grid_create)
     cl_rhokap_ref = par%clump_ndust * par%cext_dust * par%distance2cm
  else if (par%clump_nH > 0.0_wp) then
     if (par%distance2cm <= 0.0_wp) then
        if (mpar%p_rank == 0) write(*,*) 'ERROR: clump_nH requires distance_unit'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     cl_rhokap_ref = par%clump_nH * par%cext_dust * par%distance2cm
  else if (par%taumax > 0.0_wp .or. par%tauhomo > 0.0_wp) then
     !--- back-solve the peak opacity from the system-level radial target.
     if (profiles_active) then
        GF_los = LOS_geometric_quad(A_norm)
     else
        GF_los = real(N_clumps, wp) * cl_radius_max**3 / max(shell_R2_factor, tiny(1.0_wp))
     end if
     if (GF_los <= 0.0_wp) then
        if (mpar%p_rank == 0) write(*,*) &
           'ERROR: cannot back-solve clump opacity (geometric factor is zero)'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     if (par%taumax > 0.0_wp) then
        cl_rhokap_ref = par%taumax / GF_los
     else
        cl_rhokap_ref = par%tauhomo / GF_los
     end if
  else
     if (mpar%p_rank == 0) write(*,*) &
        'ERROR: specify clump_tau0, clump_ndust, clump_nH, taumax, or tauhomo'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if
  base_rhokap_in = cl_rhokap_ref

  if (mpar%p_rank == 0) then
     write(*,'(a,i14)')    ' Clumps: N_clumps  = ', N_clumps
     write(*,'(a,f12.6)')  ' Clumps: f_vol     = ', fvol_realized
     write(*,'(a,f12.5)')  ' Clumps: f_cov     = ', fcov_realized
     write(*,'(a,2f12.5)') ' Clumps: rmin/rmax = ', r_min_clump, R_sphere
     write(*,'(a,es12.4)') ' Clumps: cl_rhokap = ', cl_rhokap_ref
     if (profiles_active) then
        write(*,'(a,a)') ' Clumps: radius profile  = ', trim(par%clump_radius_profile)
        write(*,'(a,a)') ' Clumps: density profile = ', trim(par%clump_density_profile)
        write(*,'(a,a)') ' Clumps: number profile  = ', trim(par%clump_number_profile)
        write(*,'(a,f12.6,a,f12.6)') ' Clumps: r_cl range      = ', &
              base_radius_in*minval(prof_shape_radius), ' to ', cl_radius_max
     end if
  end if

  !--- shared memory for clump positions and dust properties.
  call create_shared_mem(cl_x,        [int(N_clumps)])
  call create_shared_mem(cl_y,        [int(N_clumps)])
  call create_shared_mem(cl_z,        [int(N_clumps)])
  call create_shared_mem(cl_radius,   [int(N_clumps)])
  call create_shared_mem(cl_radius2,  [int(N_clumps)])
  call create_shared_mem(cl_rhokap,   [int(N_clumps)])

  if (mpar%h_rank == 0) then
     cl_radius(:)  = cl_radius_max
     cl_radius2(:) = cl_radius_max * cl_radius_max
     cl_rhokap(:)  = cl_rhokap_ref
  end if
  call MPI_BARRIER(mpar%hostcomm, ierr)

  !--- p_rank=0 generates positions; broadcast to all h_rank=0 (one per node).
  if (mpar%h_rank == 0) then
     if (mpar%p_rank == 0) call generate_clumps()
     call MPI_BCAST(cl_x, int(N_clumps), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_y, int(N_clumps), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_z, int(N_clumps), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     if (profiles_active) then
        call MPI_BCAST(cl_radius,  int(N_clumps), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
        call MPI_BCAST(cl_radius2, int(N_clumps), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
        call MPI_BCAST(cl_rhokap,  int(N_clumps), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     end if
  end if
  call MPI_BARRIER(mpar%hostcomm, ierr)

  call build_clump_csr()
  call MPI_BARRIER(mpar%hostcomm, ierr)

  end subroutine init_clumps
  !===========================================================================

  !===========================================================================
  subroutine generate_clumps()
  !---------------------------------------------------------------------------
  ! RSA with linked-list grid acceleration.  Called only on h_rank=0.
  !---------------------------------------------------------------------------
  implicit none
  integer        :: rg, ncells_rsa
  real(kind=wp)  :: rg_cell, min_sep2_uniform
  real(kind=wp)  :: xc, yc, zc, dx, dy, dz, d2, sep_pair
  real(kind=wp)  :: r_trial, rcl_trial, cos_theta, sin_theta, phi_az
  real(kind=wp)  :: dens_factor
  real(kind=wp)  :: r_min_center, r_max_center, r_min_center2, r_max_center2
  real(kind=wp)  :: cos_cone_opening, cos_theta_min
  integer        :: ig, jg, kg, ig2, jg2, kg2, icell_rsa, jnb
  integer        :: ierr
  integer(int64) :: icl, n_attempts, stall
  logical        :: overlap
  integer, allocatable :: head(:), nxt(:)
  integer(int64), parameter :: MAX_STALL = 5000000_int64  ! RSA jamming guard

  rg         = min(512, max(32, int(real(N_clumps,wp)**(1.0_wp/3.0_wp)) + 1))
  ncells_rsa = rg**3
  rg_cell    = max(2.0_wp * sphere_R / real(rg, wp), 2.0_wp * cl_radius_max)
  rg         = max(2, int((2.0_wp * sphere_R) / rg_cell) + 1)
  rg_cell    = (2.0_wp * sphere_R) / real(rg, wp)
  ncells_rsa = rg**3
  min_sep2_uniform = (2.0_wp * cl_radius_max)**2

  if (par%clump_fully_inside) then
     if (cl_radius_max >= sphere_R) then
        if (mpar%p_rank == 0) write(*,'(a)') &
           'ERROR: par%clump_fully_inside=.true. but max clump radius >= sphere_R; '// &
           'cannot fit any clump entirely inside the medium.'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     if (r_min_clump + 2.0_wp*cl_radius_max > sphere_R) then
        if (mpar%p_rank == 0) write(*,'(a,3es12.4)') &
           'ERROR: par%clump_fully_inside=.true. but rmin + 2*cl_radius_max > rmax; '// &
           'no clump fits inside the shell. rmin, cl_radius_max, rmax = ', &
           r_min_clump, cl_radius_max, sphere_R
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
  end if

  allocate(head(ncells_rsa),   source=-1)
  allocate(nxt(int(N_clumps)), source=-1)

  if (mpar%p_rank == 0) then
     write(*,'(a,i5,a,i14,a)') ' RSA: grid ', rg, '^3, placing ', N_clumps, ' clumps...'
     if (par%clump_fully_inside) then
        write(*,'(a)') ' RSA: clump_fully_inside = .true.  (clumps must fit inside the shell)'
     else
        write(*,'(a)') ' RSA: clump_fully_inside = .false. (only centers inside the shell)'
     end if
     if (r_min_clump > 0.0_wp) &
        write(*,'(a,2f12.5)') ' RSA: shell rmin/rmax    = ', r_min_clump, sphere_R
  end if

  n_attempts = 0_int64
  stall      = 0_int64
  icl        = 0_int64
  rcl_trial  = base_radius_in

  if (par%cone_opening > 0.0_wp .and. par%cone_opening < 90.0_wp) then
     cos_cone_opening = cos(par%cone_opening * deg2rad)
  else
     cos_cone_opening = -1.0_wp
  end if

  do while (icl < N_clumps)
     n_attempts = n_attempts + 1_int64
     stall      = stall + 1_int64
     if (stall > MAX_STALL) then
        if (mpar%p_rank == 0) write(*,'(a,i0,a)') &
           'ERROR: RSA jamming -- failed to place clump ', icl+1_int64, &
           ' after 5e6 attempts.  Lower clump_f_vol/clump_f_cov or clump_radius '// &
           '(RSA is reliable only for volume filling factor < ~0.35).'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if

     if (profiles_active) then
        r_trial   = sample_clump_radius()
        rcl_trial = base_radius_in * shape_radius(r_trial)
        if (cos_cone_opening > 0.0_wp) then
           cos_theta = cos_cone_opening + (1.0_wp - cos_cone_opening) * rand_number()
           if (rand_number() < 0.5_wp) cos_theta = -cos_theta
        else
           cos_theta = 2.0_wp * rand_number() - 1.0_wp
        end if
        sin_theta = sqrt(max(0.0_wp, 1.0_wp - cos_theta*cos_theta))
        phi_az    = twopi * rand_number()
        xc = r_trial * sin_theta * cos(phi_az)
        yc = r_trial * sin_theta * sin(phi_az)
        zc = r_trial * cos_theta
        if (par%clump_fully_inside) then
           if (r_trial + rcl_trial > sphere_R)    cycle
           if (r_trial - rcl_trial < r_min_clump) cycle
           if (cos_cone_opening > 0.0_wp .and. r_trial > rcl_trial) then
              cos_theta_min = abs(cos_theta) * sqrt(1.0_wp - (rcl_trial/r_trial)**2) &
                            - sin_theta * (rcl_trial / r_trial)
              if (cos_theta_min < cos_cone_opening) cycle
           end if
        end if
     else
        if (par%clump_fully_inside) then
           r_max_center = sphere_R    - base_radius_in
           r_min_center = r_min_clump + base_radius_in
        else
           r_max_center = sphere_R
           r_min_center = r_min_clump
        end if
        r_max_center2 = r_max_center * r_max_center
        r_min_center2 = r_min_center * r_min_center
        if (cos_cone_opening > 0.0_wp) then
           do
              r_trial   = (r_min_center**3 + (r_max_center**3 - r_min_center**3) * rand_number())**(1.0_wp/3.0_wp)
              cos_theta = cos_cone_opening + (1.0_wp - cos_cone_opening) * rand_number()
              if (rand_number() < 0.5_wp) cos_theta = -cos_theta
              sin_theta = sqrt(max(0.0_wp, 1.0_wp - cos_theta*cos_theta))
              if (par%clump_fully_inside .and. r_trial > base_radius_in) then
                 cos_theta_min = abs(cos_theta) * sqrt(1.0_wp - (base_radius_in/r_trial)**2) &
                               - sin_theta * (base_radius_in / r_trial)
                 if (cos_theta_min < cos_cone_opening) cycle
              end if
              exit
           end do
           phi_az = twopi * rand_number()
           xc = r_trial * sin_theta * cos(phi_az)
           yc = r_trial * sin_theta * sin(phi_az)
           zc = r_trial * cos_theta
        else
           do
              xc = (2.0_wp * rand_number() - 1.0_wp) * r_max_center
              yc = (2.0_wp * rand_number() - 1.0_wp) * r_max_center
              zc = (2.0_wp * rand_number() - 1.0_wp) * r_max_center
              d2 = xc*xc + yc*yc + zc*zc
              if (d2 <= r_max_center2 .and. d2 >= r_min_center2) exit
           end do
        end if
     end if

     !--- cell in RSA grid (0-based)
     ig = min(rg-1, max(0, int((xc + sphere_R) / rg_cell)))
     jg = min(rg-1, max(0, int((yc + sphere_R) / rg_cell)))
     kg = min(rg-1, max(0, int((zc + sphere_R) / rg_cell)))

     !--- check 27 neighbors for overlap
     overlap = .false.
     outer: do kg2 = max(0,kg-1), min(rg-1,kg+1)
        do jg2 = max(0,jg-1), min(rg-1,jg+1)
           do ig2 = max(0,ig-1), min(rg-1,ig+1)
              icell_rsa = 1 + ig2 + jg2*rg + kg2*rg*rg
              jnb = head(icell_rsa)
              do while (jnb > 0)
                 dx = xc - cl_x(jnb);  dy = yc - cl_y(jnb);  dz = zc - cl_z(jnb)
                 d2 = dx*dx + dy*dy + dz*dz
                 if (profiles_active) then
                    sep_pair = rcl_trial + cl_radius(int(jnb,int64))
                    if (d2 < sep_pair*sep_pair) then
                       overlap = .true.;  exit outer
                    end if
                 else
                    if (d2 < min_sep2_uniform) then
                       overlap = .true.;  exit outer
                    end if
                 end if
                 jnb = nxt(jnb)
              end do
           end do
        end do
     end do outer
     if (overlap) cycle

     !--- accept
     icl   = icl + 1_int64
     stall = 0_int64
     cl_x(icl) = real(xc, dp);  cl_y(icl) = real(yc, dp);  cl_z(icl) = real(zc, dp)

     if (profiles_active) then
        cl_radius(icl)  = rcl_trial
        cl_radius2(icl) = rcl_trial * rcl_trial
        dens_factor     = shape_density(r_trial)
        cl_rhokap(icl)  = base_rhokap_in * dens_factor
     end if

     !--- insert into RSA linked list
     icell_rsa       = 1 + ig + jg*rg + kg*rg*rg
     nxt(int(icl))   = head(icell_rsa)
     head(icell_rsa) = int(icl)

     if (mod(icl, 1000000_int64) == 0_int64 .and. mpar%p_rank == 0) &
        write(*,'(a,i14,a,i14)') '   placed ', icl, ' / ', N_clumps
  end do

  if (mpar%p_rank == 0) &
     write(*,'(a,f6.1,a)') ' RSA done, acceptance rate = ', &
        real(N_clumps,wp)/real(n_attempts,wp)*100.0_wp, '%'

  deallocate(head, nxt)
  end subroutine generate_clumps
  !===========================================================================

  !===========================================================================
  subroutine build_clump_csr()
  !---------------------------------------------------------------------------
  ! Two-pass CSR construction.  All ranks allocate shared memory; h_rank=0
  ! fills it.  Caller handles barrier after return.
  !---------------------------------------------------------------------------
  implicit none
  integer        :: ncells, total_regs
  integer        :: i, j, k, imin, imax, jmin, jmax, kmin, kmax, icell
  integer(int64) :: icl
  integer,       allocatable :: cnt(:)
  integer        :: ierr

  cgx = min(512, max(32, int(real(N_clumps,wp)**(1.0_wp/3.0_wp)) + 1))
  cgy = cgx;  cgz = cgx
  ncells = cgx * cgy * cgz

  cg_xmin = -(sphere_R + cl_radius_max);  cg_ymin = cg_xmin;  cg_zmin = cg_xmin
  cg_dx   = (2.0_wp*(sphere_R + cl_radius_max)) / real(cgx, wp)
  cg_dy   = cg_dx;  cg_dz = cg_dx
  cg_inv_dx = 1.0_wp / cg_dx
  cg_inv_dy = 1.0_wp / cg_dy
  cg_inv_dz = 1.0_wp / cg_dz

  call create_shared_mem(cg_start, [ncells + 1])

  if (mpar%h_rank == 0) then
     allocate(cnt(ncells), source=0)
     do icl = 1_int64, N_clumps
        call clump_cell_range(icl, imin, imax, jmin, jmax, kmin, kmax)
        do k = kmin, kmax
           do j = jmin, jmax
              do i = imin, imax
                 icell = cg_cell_idx(i, j, k)
                 cnt(icell) = cnt(icell) + 1
              end do
           end do
        end do
     end do
     cg_start(1) = 1
     do icell = 1, ncells
        cg_start(icell+1) = cg_start(icell) + cnt(icell)
     end do
  end if
  call MPI_BARRIER(mpar%hostcomm, ierr)

  total_regs = cg_start(ncells+1) - 1
  call create_shared_mem(cg_list, [total_regs])

  if (mpar%h_rank == 0) then
     cnt(:) = 0
     do icl = 1_int64, N_clumps
        call clump_cell_range(icl, imin, imax, jmin, jmax, kmin, kmax)
        do k = kmin, kmax
           do j = jmin, jmax
              do i = imin, imax
                 icell = cg_cell_idx(i, j, k)
                 cg_list(cg_start(icell) + cnt(icell)) = int(icl, int32)
                 cnt(icell) = cnt(icell) + 1
              end do
           end do
        end do
     end do
     deallocate(cnt)
     if (mpar%p_rank == 0) write(*,'(a,i12,a,i5,a)') &
           ' CSR grid: ', total_regs, ' registrations in ', cgx, '^3 cells.'
  end if

  end subroutine build_clump_csr
  !===========================================================================

  !===========================================================================
  pure subroutine clump_cell_range(icl, imin, imax, jmin, jmax, kmin, kmax)
  integer(int64), intent(in)  :: icl
  integer,        intent(out) :: imin, imax, jmin, jmax, kmin, kmax
  real(kind=wp) :: rcl
  rcl = cl_radius(icl)
  imin = max(0, int((cl_x(icl) - cg_xmin - rcl) * cg_inv_dx))
  imax = min(cgx-1, int((cl_x(icl) - cg_xmin + rcl) * cg_inv_dx))
  jmin = max(0, int((cl_y(icl) - cg_ymin - rcl) * cg_inv_dy))
  jmax = min(cgy-1, int((cl_y(icl) - cg_ymin + rcl) * cg_inv_dy))
  kmin = max(0, int((cl_z(icl) - cg_zmin - rcl) * cg_inv_dz))
  kmax = min(cgz-1, int((cl_z(icl) - cg_zmin + rcl) * cg_inv_dz))
  end subroutine clump_cell_range
  !===========================================================================

  !===========================================================================
  pure subroutine ray_sphere_isect(ox, oy, oz, kx, ky, kz, cx, cy, cz, &
                                    icl, t_entry, t_exit, hit)
  !---------------------------------------------------------------------------
  ! Ray-sphere intersection.  Ray P(t) = origin + t*dir; sphere center
  ! (cx,cy,cz), squared-radius cl_radius2(icl).  hit=.true. & t_entry<=t_exit
  ! when t_exit > 0.
  !---------------------------------------------------------------------------
  real(kind=wp),  intent(in)  :: ox, oy, oz, kx, ky, kz, cx, cy, cz
  integer(int64), intent(in)  :: icl
  real(kind=wp),  intent(out) :: t_entry, t_exit
  logical,        intent(out) :: hit
  real(kind=wp) :: rx, ry, rz, b, disc
  rx = ox - cx;  ry = oy - cy;  rz = oz - cz
  b    = rx*kx + ry*ky + rz*kz
  disc = b*b - (rx*rx + ry*ry + rz*rz) + cl_radius2(icl)
  if (disc < 0.0_wp) then
     hit = .false.;  t_entry = 0.0_wp;  t_exit = 0.0_wp
  else
     disc    = sqrt(disc)
     t_entry = -b - disc
     t_exit  = -b + disc
     hit     = (t_exit > 0.0_wp)
  end if
  end subroutine ray_sphere_isect
  !===========================================================================

  !===========================================================================
  subroutine find_next_clump(xp, yp, zp, kx, ky, kz, skip_icl, t_max, &
                              t_entry, t_exit, icl_found, found)
  !---------------------------------------------------------------------------
  ! DDA through CSR grid to find the nearest clump hit along the ray
  !   P(t) = (xp,yp,zp) + t*(kx,ky,kz),  0 < t <= t_max.
  ! Clump skip_icl (> 0) is excluded (the one just exited).
  ! Pattern: Amanatides & Woo (1987).
  !---------------------------------------------------------------------------
  real(kind=wp),  intent(in)  :: xp, yp, zp, kx, ky, kz
  integer(int64), intent(in)  :: skip_icl
  real(kind=wp),  intent(in)  :: t_max
  real(kind=wp),  intent(out) :: t_entry, t_exit
  integer(int64), intent(out) :: icl_found
  logical,        intent(out) :: found

  integer        :: ci, cj, ck, si, sj, sk, icell, ip
  integer(int64) :: icl
  real(kind=wp)  :: tx, ty, tz, delx, dely, delz, d
  real(kind=wp)  :: te, tx2, best_te, best_tx2
  integer(int64) :: best_icl
  logical        :: hit

  found    = .false.
  best_te  = hugest
  best_tx2 = 0.0_wp
  best_icl = 0_int64
  d        = 0.0_wp
  t_entry  = 0.0_wp
  t_exit   = 0.0_wp
  icl_found = 0_int64

  ci = max(0, min(cgx-1, int((xp - cg_xmin) * cg_inv_dx)))
  cj = max(0, min(cgy-1, int((yp - cg_ymin) * cg_inv_dy)))
  ck = max(0, min(cgz-1, int((zp - cg_zmin) * cg_inv_dz)))

  if (kx > 0.0_wp) then
     si   =  1
     tx   = ((cg_xmin + real(ci+1,wp)*cg_dx) - xp) / kx
     delx =  cg_dx / kx
  else if (kx < 0.0_wp) then
     si   = -1
     tx   = ((cg_xmin + real(ci,wp)*cg_dx) - xp) / kx
     delx = -cg_dx / kx
  else
     si = 0;  tx = hugest;  delx = hugest
  end if

  if (ky > 0.0_wp) then
     sj   =  1
     ty   = ((cg_ymin + real(cj+1,wp)*cg_dy) - yp) / ky
     dely =  cg_dy / ky
  else if (ky < 0.0_wp) then
     sj   = -1
     ty   = ((cg_ymin + real(cj,wp)*cg_dy) - yp) / ky
     dely = -cg_dy / ky
  else
     sj = 0;  ty = hugest;  dely = hugest
  end if

  if (kz > 0.0_wp) then
     sk   =  1
     tz   = ((cg_zmin + real(ck+1,wp)*cg_dz) - zp) / kz
     delz =  cg_dz / kz
  else if (kz < 0.0_wp) then
     sk   = -1
     tz   = ((cg_zmin + real(ck,wp)*cg_dz) - zp) / kz
     delz = -cg_dz / kz
  else
     sk = 0;  tz = hugest;  delz = hugest
  end if

  do while(.true.)
     if (d > best_te .or. d > t_max) exit

     icell = cg_cell_idx(ci, cj, ck)
     do ip = cg_start(icell), cg_start(icell+1) - 1
        icl = int(cg_list(ip), int64)
        if (icl == skip_icl) cycle
        call ray_sphere_isect(xp, yp, zp, kx, ky, kz, &
             real(cl_x(icl),wp), real(cl_y(icl),wp), real(cl_z(icl),wp), &
             icl, te, tx2, hit)
        if (hit .and. tx2 > 0.0_wp .and. te < best_te .and. &
            (te > 0.0_wp .or. icl /= skip_icl)) then
           best_te  = te
           best_tx2 = tx2
           best_icl = icl
        end if
     end do

     if (tx <= ty .and. tx <= tz) then
        d  = tx
        ci = ci + si;  if (ci < 0 .or. ci >= cgx) exit
        tx = tx + delx
     else if (ty <= tz) then
        d  = ty
        cj = cj + sj;  if (cj < 0 .or. cj >= cgy) exit
        ty = ty + dely
     else
        d  = tz
        ck = ck + sk;  if (ck < 0 .or. ck >= cgz) exit
        tz = tz + delz
     end if
  end do

  if (best_icl > 0_int64 .and. best_te <= t_max) then
     found     = .true.
     t_entry   = best_te
     t_exit    = min(best_tx2, t_max)
     icl_found = best_icl
  end if

  end subroutine find_next_clump
  !===========================================================================

  !===========================================================================
  pure real(kind=wp) function clump_exit_dist(xp, yp, zp, kx, ky, kz, icl)
  !---------------------------------------------------------------------------
  ! Distance from (xp,yp,zp) along (kx,ky,kz) to the exit of clump icl.
  ! Assumes the photon is currently inside the clump.
  !---------------------------------------------------------------------------
  real(kind=wp),  intent(in) :: xp, yp, zp, kx, ky, kz
  integer(int64), intent(in) :: icl
  real(kind=wp) :: rx, ry, rz, b, disc
  rx = xp - real(cl_x(icl),wp)
  ry = yp - real(cl_y(icl),wp)
  rz = zp - real(cl_z(icl),wp)
  b    = rx*kx + ry*ky + rz*kz
  disc = b*b - (rx*rx + ry*ry + rz*rz) + cl_radius2(icl)
  if (disc < 0.0_wp) disc = 0.0_wp
  clump_exit_dist = max(0.0_wp, -b + sqrt(disc))
  end function clump_exit_dist
  !===========================================================================

  !===========================================================================
  pure real(kind=wp) function sphere_exit_dist(xp, yp, zp, kx, ky, kz)
  !---------------------------------------------------------------------------
  ! Distance to the exit of the outer sphere (radius sphere_R) from (xp,yp,zp).
  !---------------------------------------------------------------------------
  real(kind=wp), intent(in) :: xp, yp, zp, kx, ky, kz
  real(kind=wp) :: b, disc
  b    = xp*kx + yp*ky + zp*kz
  disc = b*b - (xp*xp + yp*yp + zp*zp) + sphere_R*sphere_R
  if (disc < 0.0_wp) disc = 0.0_wp
  sphere_exit_dist = max(0.0_wp, -b + sqrt(disc))
  end function sphere_exit_dist
  !===========================================================================

  !===========================================================================
  subroutine check_has_overlap()
  !---------------------------------------------------------------------------
  ! Scan clump pairs via the CSR grid for any overlap.  Sets has_overlap on
  ! rank 0 and broadcasts.  RSA-generated populations never overlap; a loaded
  ! population may.  The dust raytrace assumes non-overlap, so a warning is
  ! emitted if overlap is found.
  !---------------------------------------------------------------------------
  implicit none
  integer        :: imin, imax, jmin, jmax, kmin, kmax
  integer        :: i, j, k, icell, ip
  integer(int64) :: icl_i, icl_j
  real(kind=dp)  :: dx, dy, dz, dist2, r_sum
  integer        :: ierr

  has_overlap = .false.
  if (mpar%h_rank == 0) then
     outer: do icl_i = 1_int64, N_clumps
        call clump_cell_range(icl_i, imin, imax, jmin, jmax, kmin, kmax)
        do k = kmin, kmax
           do j = jmin, jmax
              do i = imin, imax
                 icell = cg_cell_idx(i, j, k)
                 do ip = cg_start(icell), cg_start(icell+1) - 1
                    icl_j = int(cg_list(ip), int64)
                    if (icl_j <= icl_i) cycle
                    dx    = cl_x(icl_j) - cl_x(icl_i)
                    dy    = cl_y(icl_j) - cl_y(icl_i)
                    dz    = cl_z(icl_j) - cl_z(icl_i)
                    dist2 = dx*dx + dy*dy + dz*dz
                    r_sum = real(cl_radius(icl_i) + cl_radius(icl_j), dp)
                    if (dist2 < r_sum*r_sum) then
                       has_overlap = .true.
                       exit outer
                    end if
                 end do
              end do
           end do
        end do
     end do outer
  end if
  call MPI_BCAST(has_overlap, 1, MPI_LOGICAL, 0, mpar%SAME_HRANK_COMM, ierr)
  if (mpar%p_rank == 0 .and. has_overlap) write(*,'(a)') &
     ' Clumps: WARNING -- overlapping clumps detected; dust raytrace assumes '// &
     'non-overlap (overlap-aware path is not ported).'
  end subroutine check_has_overlap
  !===========================================================================

  !===========================================================================
  subroutine active_set_at_point(xp, yp, zp, active, n_active)
  !---------------------------------------------------------------------------
  ! Return all clump indices containing point (xp,yp,zp).  active(1:n_active)
  ! are the 1-based clump indices.  Used to set photon%icell_clump at birth.
  !---------------------------------------------------------------------------
  implicit none
  real(kind=wp),  intent(in)  :: xp, yp, zp
  integer(int64), intent(out) :: active(:)
  integer,        intent(out) :: n_active

  integer        :: ci, cj, ck, i, j, k, icell, ip
  integer(int64) :: icl
  real(kind=dp)  :: rx, ry, rz

  n_active = 0
  if (cgx <= 0) return
  ci = max(0, min(cgx-1, int((xp - cg_xmin) * cg_inv_dx)))
  cj = max(0, min(cgy-1, int((yp - cg_ymin) * cg_inv_dy)))
  ck = max(0, min(cgz-1, int((zp - cg_zmin) * cg_inv_dz)))

  do k = max(0, ck-1), min(cgz-1, ck+1)
     do j = max(0, cj-1), min(cgy-1, cj+1)
        do i = max(0, ci-1), min(cgx-1, ci+1)
           icell = cg_cell_idx(i, j, k)
           do ip = cg_start(icell), cg_start(icell+1) - 1
              icl = int(cg_list(ip), int64)
              rx  = real(xp, dp) - cl_x(icl)
              ry  = real(yp, dp) - cl_y(icl)
              rz  = real(zp, dp) - cl_z(icl)
              if (rx*rx + ry*ry + rz*rz <= real(cl_radius2(icl), dp)) then
                 if (n_active == 0 .or. all(active(1:n_active) /= icl)) then
                    if (n_active < size(active)) then
                       n_active = n_active + 1
                       active(n_active) = icl
                    end if
                 end if
              end if
           end do
        end do
     end do
  end do
  end subroutine active_set_at_point
  !===========================================================================

  !===========================================================================
  subroutine destroy_clumps()
  implicit none
  ! destroy_shared_mem_all() frees the windows collectively; then nullify.
  call destroy_shared_mem_all()
  nullify(cl_x, cl_y, cl_z)
  nullify(cl_radius, cl_radius2, cl_rhokap)
  nullify(cg_start, cg_list)
  N_clumps = 0_int64
  end subroutine destroy_clumps
  !===========================================================================

  !===========================================================================
  subroutine compute_clump_scalars(tauhomo_out, taumax_out)
  !---------------------------------------------------------------------------
  ! System-level dust optical-depth scalars from the clump arrays.
  !   tauhomo = Sum_i rhokap_i * r_i^3 / (R^2 + R*r0 + r0^2)
  !             (uniform-smear of clump opacity over the shell, traversed
  !              radially)
  !   taumax  = Sum_i rhokap_i * r_i^3 / (3 * max(d_i^2, r_i^2))
  !             (small-angle expected radial tau from the origin; the cap on
  !              d^2 avoids the singularity for a clump straddling the origin)
  ! Uniform-isotropic limit -> both reduce to N r_cl^3 rhokap / (R^2+R r0+r0^2)
  !  = (4/3) f_cov tau_per_clump.  Computed on p_rank=0, then broadcast.
  !---------------------------------------------------------------------------
  implicit none
  real(kind=wp), intent(out) :: tauhomo_out, taumax_out
  integer        :: ierr
  integer(int64) :: i
  real(kind=wp)  :: di2, shell_R2_factor, w_homo, w_max
  real(kind=wp)  :: vals(2)

  shell_R2_factor = sphere_R**2 + sphere_R*r_min_clump + r_min_clump**2
  if (shell_R2_factor <= 0.0_wp) shell_R2_factor = sphere_R**2

  vals(:) = 0.0_wp
  if (mpar%p_rank == 0 .and. N_clumps > 0_int64) then
     do i = 1_int64, N_clumps
        di2     = max(cl_x(i)**2 + cl_y(i)**2 + cl_z(i)**2, cl_radius(i)**2)
        w_homo  = cl_radius(i)**3 / shell_R2_factor
        w_max   = cl_radius(i)**3 / (3.0_wp * di2)
        vals(1) = vals(1) + cl_rhokap(i) * w_homo
        vals(2) = vals(2) + cl_rhokap(i) * w_max
     end do
  end if
  call MPI_BCAST(vals, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  tauhomo_out = vals(1)
  taumax_out  = vals(2)
  end subroutine compute_clump_scalars
  !===========================================================================

  !===========================================================================
  subroutine write_clumps_info(fname)
  !---------------------------------------------------------------------------
  ! Write clump positions and dust opacity to a FITS/HDF5 binary table.
  ! Called only on p_rank == 0 after init_clumps().
  !   HDU 2 (BinTable): columns X, Y, Z [code units]; R_CLUMP, RHOKAP written
  !   only when non-constant (otherwise captured by header keywords).
  !---------------------------------------------------------------------------
  use define
  use iofile_mod
  implicit none
  character(len=*), intent(in) :: fname

  type(io_file_type) :: iofh
  integer :: status, bitpix
  real(kind=real64), allocatable :: tmp(:)
  integer(int64) :: ncl, i
  real(kind=wp)  :: f_vol_actual, f_cov_actual
  real(kind=wp)  :: rcl_min, rcl_max, rcl_mean, tau_mean
  real(kind=wp)  :: kap_min_w, kap_max_w, kap_mean_w
  logical        :: write_radius, write_rhokap
  real(kind=wp), parameter :: const_tol = 1.0e-3_wp

  status = 0
  ncl    = N_clumps
  bitpix = -32

  if (clumps_from_file .or. profiles_active) then
     f_vol_actual = 0.0_wp
     f_cov_actual = 0.0_wp
     do i = 1_int64, ncl
        f_vol_actual = f_vol_actual + cl_radius(i)**3
        f_cov_actual = f_cov_actual + cl_radius(i)**2
     end do
     f_vol_actual = f_vol_actual / max(sphere_R**3 - r_min_clump**3, tiny(1.0_wp))
     f_cov_actual = 0.75_wp * f_cov_actual / &
                    max(sphere_R**2 + sphere_R*r_min_clump + r_min_clump**2, tiny(1.0_wp))
  else
     f_vol_actual = real(ncl,wp) * cl_radius_max**3 / &
                    max(sphere_R**3 - r_min_clump**3, tiny(1.0_wp))
     f_cov_actual = 0.75_wp * real(ncl,wp) * cl_radius_max**2 / &
                    max(sphere_R**2 + sphere_R*r_min_clump + r_min_clump**2, tiny(1.0_wp))
  end if

  rcl_min  = huge(1.0_wp);  rcl_max = 0.0_wp;  rcl_mean = 0.0_wp;  tau_mean = 0.0_wp
  do i = 1_int64, ncl
     rcl_min  = min(rcl_min,  cl_radius(i))
     rcl_max  = max(rcl_max,  cl_radius(i))
     rcl_mean = rcl_mean + cl_radius(i)
     tau_mean = tau_mean + cl_rhokap(i) * cl_radius(i)
  end do
  if (ncl > 0_int64) then
     rcl_mean = rcl_mean / real(ncl, wp)
     tau_mean = tau_mean / real(ncl, wp)
  end if

  if (mpar%p_rank == 0) then
     write(*,'(a,es12.4)') ' Clumps: realized f_vol      = ', f_vol_actual
     write(*,'(a,es12.4)') ' Clumps: realized f_cov      = ', f_cov_actual
     write(*,'(a,3es12.4)')' Clumps: r_cl min/mean/max   = ', rcl_min, rcl_mean, rcl_max
     write(*,'(a,es12.4)') ' Clumps: mean tau_per_clump  = ', tau_mean
  end if

  call io_open_new(iofh, trim(fname), status)
  if (status /= 0) then
     write(*,*) 'WARNING: write_clumps_info: cannot open ', trim(fname)
     return
  end if

  allocate(tmp(ncl))
  tmp = cl_x(1:ncl)
  call io_write_table_column(iofh, 'X', tmp, status, bitpix)
  tmp = cl_y(1:ncl)
  call io_write_table_column(iofh, 'Y', tmp, status, bitpix)
  tmp = cl_z(1:ncl)
  call io_write_table_column(iofh, 'Z', tmp, status, bitpix)

  kap_min_w  = minval(cl_rhokap(1:ncl))
  kap_max_w  = maxval(cl_rhokap(1:ncl))
  kap_mean_w = sum(cl_rhokap(1:ncl)) / max(real(ncl,wp), 1.0_wp)

  write_radius = (rcl_max  - rcl_min ) > const_tol * max(abs(rcl_mean ), tiny(1.0_wp))
  write_rhokap = (kap_max_w - kap_min_w) > const_tol * max(abs(kap_mean_w), tiny(1.0_wp))

  if (write_radius) then
     tmp = cl_radius(1:ncl)
     call io_write_table_column(iofh, 'R_CLUMP', tmp, status, bitpix)
  end if
  if (write_rhokap) then
     tmp = cl_rhokap(1:ncl)
     call io_write_table_column(iofh, 'RHOKAP', tmp, status, bitpix)
  end if
  deallocate(tmp)

  call io_put_keyword(iofh, 'N_CLUMPS', ncl,                'number of clumps (realized)',         status)
  call io_put_keyword(iofh, 'SPHERE_R', sphere_R,           'outer sphere radius [code units]',     status)
  call io_put_keyword(iofh, 'RMIN',     r_min_clump,        'inner placement radius [code units]',  status)
  call io_put_keyword(iofh, 'CL_RAD',   cl_radius_max,      'clump radius (max) [code units]',      status)
  call io_put_keyword(iofh, 'F_VOL',    f_vol_actual,       'volume filling factor (realized)',     status)
  call io_put_keyword(iofh, 'F_COV',    f_cov_actual,       'covering factor (realized)',           status)
  call io_put_keyword(iofh, 'TAU0',     par%clump_tau0,     'dust tau of one clump (center->surf)', status)
  call io_put_keyword(iofh, 'RHOKAP',   cl_rhokap_ref,      'dust opacity/code-unit (reference)',   status)
  call io_put_keyword(iofh, 'RMAX',     par%rmax,           'outer sphere radius input (par%rmax)', status)
  call io_put_keyword(iofh, 'IN_FCOV',  par%clump_f_cov,    'covering factor (input)',              status)
  call io_put_keyword(iofh, 'IN_FVOL',  par%clump_f_vol,    'volume filling factor (input)',        status)
  call io_put_keyword(iofh, 'IN_NCL',   par%clump_N_clumps, 'N_clumps (input)',                     status)
  call io_put_keyword(iofh, 'IN_NDUST', par%clump_ndust,    'clump dust density input',             status)
  call io_put_keyword(iofh, 'IN_NH',    par%clump_nH,       'clump nH density input [cm^-3]',       status)
  call io_put_keyword(iofh, 'DISTUNIT', par%distance_unit,  'Distance Unit',                        status)
  call io_put_keyword(iofh, 'DIST_CM',  par%distance2cm,    'Distance Unit (cm)',                   status)

  call io_close(iofh, status)
  if (mpar%p_rank == 0) write(*,'(2a)') ' Clumps saved to ', trim(fname)
  end subroutine write_clumps_info
  !===========================================================================

  !===========================================================================
  subroutine read_perclump_or_keyword(iofh, colnames, keyname, scratch, ncl, dst, found_name)
  !---------------------------------------------------------------------------
  ! Try each comma-separated column name in `colnames`; the first present is
  ! read into `dst`.  If none is present, fall back to header keyword
  ! `keyname` and fill `dst` uniformly.  Optional found_name returns the used
  ! column name (or empty for the keyword fallback).
  !---------------------------------------------------------------------------
  use iofile_mod
  implicit none
  type(io_file_type), intent(in)             :: iofh
  character(len=*),   intent(in)             :: colnames, keyname
  real(real32),       intent(inout)          :: scratch(:)
  integer(int64),     intent(in)             :: ncl
  real(kind=wp),      intent(out)            :: dst(:)
  character(len=*),   intent(out),  optional :: found_name

  integer            :: status, colnum, ierr, ks, ke, lstr
  real(kind=wp)      :: keyval
  character(len=80)  :: cmt
  character(len=64)  :: candidate
  character(len=64)  :: tried
  logical            :: found

  found = .false.
  tried = ''
  if (present(found_name)) found_name = ''
  lstr  = len_trim(colnames)
  ks    = 1
  do while (ks <= lstr .and. .not. found)
     ke = index(colnames(ks:lstr), ',')
     if (ke == 0) then
        candidate = adjustl(colnames(ks:lstr))
        ks = lstr + 1
     else
        candidate = adjustl(colnames(ks:ks+ke-2))
        ks = ks + ke
     end if
     if (len_trim(candidate) == 0) cycle
     tried = trim(tried)//' '//trim(candidate)
     status = 0
     call io_get_column_number(iofh, trim(candidate), colnum, status)
     if (status == 0) then
        call io_read_table_column(iofh, colnum, scratch, status)
        dst(1:ncl) = real(scratch(1:ncl), wp)
        if (present(found_name)) found_name = trim(candidate)
        found = .true.
     end if
  end do

  if (.not. found) then
     status = 0
     call io_get_keyword(iofh, trim(keyname), keyval, status, cmt)
     if (status /= 0) then
        write(*,'(5a)') 'ERROR: clump_input_file is missing every column in {', &
             trim(adjustl(tried)), ' } and header keyword ', trim(keyname), &
             '; cannot fall back to a constant value'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     dst(1:ncl) = keyval
  end if
  end subroutine read_perclump_or_keyword
  !===========================================================================

  !===========================================================================
  subroutine read_clumps_info(fname, R_sphere)
  !---------------------------------------------------------------------------
  ! Read a clump population from a FITS/HDF5 file written by write_clumps_info
  ! (or the python make_clumps.py generator).  Reads X, Y, Z (mandatory) and
  ! R_CLUMP, RHOKAP (for each clump, or via header keywords CL_RAD / RHOKAP).  A
  ! DENSITY/DENS column is converted to opacity via cext_dust * distance2cm.
  ! Builds the CSR grid.  Sets clumps_from_file = .true.
  !---------------------------------------------------------------------------
  use define
  use iofile_mod
  use mpi
  implicit none
  character(len=*), intent(in) :: fname
  real(kind=wp),    intent(in) :: R_sphere

  type(io_file_type) :: iofh
  integer        :: status, ierr, colnum
  integer(int64) :: ncl, i
  real(kind=real32), allocatable :: tmp(:)
  character(len=80)  :: cmt
  character(len=64)  :: rhokap_colname
  character(len=128) :: distunit_file
  real(kind=wp)      :: distcm_file
  logical            :: rhokap_is_density, distance_unit_ok, user_set_distunit
  character(len=:), allocatable :: resolved_fname

  status   = 0
  sphere_R = R_sphere
  base_radius_in   = par%clump_radius
  profiles_active  = .false.
  clumps_from_file = .true.
  r_min_clump      = max(0.0_wp, par%rmin)

  resolved_fname = io_resolve_filename(fname)

  if (mpar%p_rank == 0) then
     write(*,'(2a)') ' Clumps: reading from ', resolved_fname
     call io_open_old(iofh, resolved_fname, status)
     if (status /= 0) then
        write(*,'(2a)') 'ERROR: cannot open clump_input_file: ', resolved_fname
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     call io_move_to_next_section(iofh, status)
     if (status /= 0) then
        write(*,'(a)') 'ERROR: clump_input_file: cannot reach binary-table HDU 2'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     call io_get_keyword(iofh, 'N_CLUMPS', ncl, status, cmt)
     if (status /= 0 .or. ncl <= 0_int64) then
        write(*,'(a)') 'ERROR: clump_input_file: N_CLUMPS keyword missing or <= 0'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
  end if
  call MPI_BCAST(ncl, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
  N_clumps = ncl

  call create_shared_mem(cl_x,       [int(N_clumps)])
  call create_shared_mem(cl_y,       [int(N_clumps)])
  call create_shared_mem(cl_z,       [int(N_clumps)])
  call create_shared_mem(cl_radius,  [int(N_clumps)])
  call create_shared_mem(cl_radius2, [int(N_clumps)])
  call create_shared_mem(cl_rhokap,  [int(N_clumps)])

  if (mpar%p_rank == 0) then
     allocate(tmp(ncl))

     call io_get_column_number(iofh, 'X', colnum, status)
     call io_read_table_column(iofh, colnum, tmp, status)
     cl_x = real(tmp, dp)
     call io_get_column_number(iofh, 'Y', colnum, status)
     call io_read_table_column(iofh, colnum, tmp, status)
     cl_y = real(tmp, dp)
     call io_get_column_number(iofh, 'Z', colnum, status)
     call io_read_table_column(iofh, colnum, tmp, status)
     cl_z = real(tmp, dp)

     call read_perclump_or_keyword(iofh, 'R_CLUMP',             'CL_RAD', tmp, ncl, cl_radius)
     call read_perclump_or_keyword(iofh, 'RHOKAP,DENSITY,DENS', 'RHOKAP', tmp, ncl, cl_rhokap, &
                                   found_name=rhokap_colname)

     status        = 0
     distunit_file = ''
     call io_get_keyword(iofh, 'DISTUNIT', distunit_file, status, cmt)
     if (status /= 0) distunit_file = ''
     status        = 0
     distcm_file   = 0.0_wp
     call io_get_keyword(iofh, 'DIST_CM', distcm_file, status, cmt)
     if (status /= 0) distcm_file = 0.0_wp

     deallocate(tmp)
     call io_close(iofh, status)

     user_set_distunit = (len_trim(par%distance_unit) > 0 .and. par%distance2cm > 1.0_wp)
     if (.not. user_set_distunit) then
        if (len_trim(distunit_file) > 0) par%distance_unit = trim(distunit_file)
        if (distcm_file > 0.0_wp)        par%distance2cm   = distcm_file
     end if

     do i = 1_int64, ncl
        cl_radius2(i) = cl_radius(i) * cl_radius(i)
     end do

     !--- DENSITY/DENS column -> dust opacity per code unit (MoCafe convention).
     rhokap_is_density = (trim(rhokap_colname) == 'DENSITY' .or. trim(rhokap_colname) == 'DENS')
     distance_unit_ok  = (len_trim(par%distance_unit) > 0 .or. par%distance2cm > 1.0_wp)
     if (rhokap_is_density) then
        if (distance_unit_ok) then
           do i = 1_int64, ncl
              cl_rhokap(i) = cl_rhokap(i) * par%cext_dust * par%distance2cm
           end do
           write(*,'(3a)') ' Clumps: ', trim(rhokap_colname), &
                ' column converted to dust opacity via cext_dust * distance2cm.'
        else
           write(*,'(3a)')  ' Clumps: WARNING -- column ', trim(rhokap_colname), &
                ' found but distance_unit not set; treating values as RHOKAP.'
        end if
     end if
  end if

  call MPI_BCAST(par%distance_unit, len(par%distance_unit), MPI_CHARACTER,        0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(par%distance2cm,   1,                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  call MPI_BARRIER(mpar%hostcomm, ierr)
  if (mpar%h_rank == 0) then
     call MPI_BCAST(cl_x,       int(ncl), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_y,       int(ncl), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_z,       int(ncl), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_radius,  int(ncl), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_radius2, int(ncl), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_rhokap,  int(ncl), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
  end if
  call MPI_BARRIER(mpar%hostcomm, ierr)

  cl_radius_max  = maxval(cl_radius(1:ncl))
  cl_rhokap_ref  = cl_rhokap(1)
  base_rhokap_in = cl_rhokap_ref

  call rescale_loaded_clumps_to_target()

  if (mpar%p_rank == 0) then
     write(*,'(a,i14)')    ' Clumps: N_clumps  = ', N_clumps
     write(*,'(a,es12.4)') ' Clumps: cl_rhokap = ', cl_rhokap_ref
     write(*,'(a,2es12.4)')' Clumps: r_cl min/max  = ', minval(cl_radius(1:ncl)), cl_radius_max
  end if

  call build_clump_csr()
  call MPI_BARRIER(mpar%hostcomm, ierr)

  end subroutine read_clumps_info
  !===========================================================================

  !===========================================================================
  subroutine rescale_loaded_clumps_to_target()
  !---------------------------------------------------------------------------
  ! When par%taumax or par%tauhomo (priority order) is supplied with a loaded
  ! clump file, multiply every cl_rhokap(i) by one factor so the realized
  ! scalar matches the target.  Otherwise leave cl_rhokap unchanged.
  !---------------------------------------------------------------------------
  implicit none
  integer           :: ierr
  integer(int64)    :: i
  real(kind=wp)     :: realized, target_val, alpha
  real(kind=wp)     :: tauhomo_r, taumax_r
  character(len=12) :: which_target

  if (par%taumax <= 0.0_wp .and. par%tauhomo <= 0.0_wp) return

  call compute_clump_scalars(tauhomo_r, taumax_r)

  if (par%taumax > 0.0_wp) then
     realized = taumax_r;   target_val = par%taumax;   which_target = 'taumax'
  else
     realized = tauhomo_r;  target_val = par%tauhomo;  which_target = 'tauhomo'
  end if

  if (realized <= 0.0_wp) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        ' Clumps: WARNING -- realized tau from loaded clumps is zero; cannot rescale.'
     return
  end if
  alpha = target_val / realized

  if (mpar%h_rank == 0) then
     do i = 1_int64, N_clumps
        cl_rhokap(i) = cl_rhokap(i) * alpha
     end do
  end if
  call MPI_BARRIER(mpar%hostcomm, ierr)

  cl_rhokap_ref  = cl_rhokap_ref * alpha
  base_rhokap_in = cl_rhokap_ref

  if (mpar%p_rank == 0) write(*,'(3a,es12.4,a,es12.4)') &
     ' Clumps: rescaled to par%', trim(which_target), ' = ', target_val, ', alpha = ', alpha
  end subroutine rescale_loaded_clumps_to_target
  !===========================================================================

end module clump_mod
