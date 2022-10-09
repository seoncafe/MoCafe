module define
  use, intrinsic :: iso_fortran_env, only: real32, real64, real128, int32, int64
public

!! kind parameter for single and double precision (6 and 15 significant decimal digits)
!! precision of 32-, 64-, and 128-bit reals
!  integer, parameter :: sp = selected_real_kind(6, 37)
!  integer, parameter :: dp = selected_real_kind(15, 307)
!  integer, parameter :: qp = selected_real_kind(33, 4931)
!! precision of 32-, 64- bit integers
!  integer, parameter :: i4b = selected_int_kind(r=9)
!  integer, parameter :: i8b = selected_int_kind(r=18)
  integer, parameter :: sp  = real32
  integer, parameter :: dp  = real64
  integer, parameter :: qp  = real128
  integer, parameter :: i4b = int32
  integer, parameter :: i8b = int64

! working precision
! Note that (1) iwp (integer working precision) applies only to seed of KISS random number generator.
!           (2) nphotons is always 64-bit integer.
!           (3) other integer variables are always 32-bit integers.
!  integer, parameter :: wp  = sp
!  integer, parameter :: iwp = i4b
  integer, parameter :: wp  = dp
  integer, parameter :: iwp = i8b

! numerical values
  real(kind=wp), parameter :: pi      = 3.141592653589793238462643383279502884197_wp
  real(kind=wp), parameter :: twopi   = 6.283185307179586476925286766559005768394_wp
  real(kind=wp), parameter :: fourpi  = 12.56637061435917295385057353311801153679_wp
  real(kind=wp), parameter :: halfpi  = pi/2.0_wp
  real(kind=wp), parameter :: sqrttwo = sqrt(2.0_wp)

! conversion factor (radian to degree, degree to radian)
  real(kind=wp), parameter :: rad2deg = 180.0_wp/pi
  real(kind=wp), parameter :: deg2rad = pi/180.0_wp

! speed of light in km/s
  real(kind=wp), parameter :: speedc  = 2.99792458e5_wp

! h_planck = Planck's constant in m^2 kg/s
! massH    = Hydrogen mass in kg
  real(kind=wp), parameter :: h_planck = 6.62607004e-34_wp
  real(kind=wp), parameter :: massH    = 1.6737236e-27_wp
  real(kind=wp), parameter :: g0       = h_planck/massH

! distances
  real(kind=wp), parameter :: pc2cm  = 3.0856776e18_wp
  real(kind=wp), parameter :: kpc2cm = pc2cm * 1e3_wp
  real(kind=wp), parameter :: au2cm  = 1.4960e13_wp
  real(kind=wp), parameter :: ang2m  = 1.0e-10_wp
  real(kind=wp), parameter :: ang2km = 1.0e-13_wp

! dust mass per hydrogen nuclei (g/H)
  real(kind=wp), parameter :: h2dust = 1.87e-26_wp

! kappa_v : dust extinction cross-section (cm^2/g) at V-band
  real(kind=wp), parameter :: kappa_v = 25992.6_wp

! tinest = the smallest positive number
! eps    = the least positive number
!          that added to 1 returns a number that is greater than 1
! hugest = the hugest positive number
  real(kind=wp), parameter :: tinest = tiny(0.0_wp)
  real(kind=wp), parameter :: eps    = epsilon(0.0_wp)
  real(kind=sp), parameter :: eps_sp = epsilon(0.0_sp)
  real(kind=dp), parameter :: eps_dp = epsilon(0.0_dp)
  real(kind=wp), parameter :: hugest = huge(1.0_wp)

! define NaN
  real(real64),  parameter :: nan64 = transfer(-2251799813685248_int64, 1._real64)
  real(real32),  parameter :: nan32 = transfer(-4194304_int32, 1._real32)

! photon type
  type photon_type
     real(kind=wp) :: x,y,z
     real(kind=wp) :: kx,ky,kz
     real(kind=wp) :: nx,ny,nz
     real(kind=wp) :: mx,my,mz
     integer       :: icell,jcell,kcell
     integer       :: nscatt
     real(kind=wp) :: wgt
     logical       :: inside
     real(kind=wp) :: albedo, hgg
     ! Stokes parameters
     real(kind=wp) :: I,Q,U,V
  end type photon_type

! cylindrical grid type
! nr, np, nz          - number of cells
! rface, pface, zface - locations of cell faces
! opacity             - opacity per unit length
  type grid_type
     integer :: nx,ny,nz
     integer :: i0 = 0, j0 = 0, k0 = 0
     real(kind=wp) :: rmax
     real(kind=wp) :: xmin,xmax,xrange
     real(kind=wp) :: ymin,ymax,yrange
     real(kind=wp) :: zmin,zmax,zrange
     real(kind=wp) :: dx,dy,dz
     real(kind=wp), pointer :: xface(:)      => null()
     real(kind=wp), pointer :: yface(:)      => null()
     real(kind=wp), pointer :: zface(:)      => null()
     real(kind=wp), pointer :: rhokap(:,:,:) => null()
  end type grid_type

! input parameters
  type params_type
     integer(kind=int64) :: nphotons   = 1e5
     integer             :: nprint     = 1e7
     real(kind=wp)       :: no_photons = 1e5
     real(kind=wp)       :: no_print   = 1e7
     integer       :: iseed        = 0
     real(kind=wp) :: luminosity   = 1.0
     real(kind=wp) :: tauhomo      = -999.0
     real(kind=wp) :: taumax       = -999.0
     real(kind=wp) :: lambda0      = 6563.0_wp
     !--- grid geometry parameters
     logical       :: xyz_symmetry    = .false.
     logical       :: xy_periodic     = .false.
     logical       :: z_symmetry      = .false.
     logical       :: sightline_tau   = .false.
     integer       :: nx   = 1
     integer       :: ny   = 1
     integer       :: nz   = 11
     real(kind=wp) :: xmax = 1.0_wp
     real(kind=wp) :: ymax = 1.0_wp
     real(kind=wp) :: zmax = 1.0_wp
     !--- to simulate spherical geometry
     integer       :: nr   = -999
     real(kind=wp) :: rmax = -999.0
     real(kind=wp) :: rmin = 0.0
     !real(kind=wp) :: rsource = 0.0_wp
     !--- continuum parameters
     real(kind=wp) :: source_zscale        = nan64
     real(kind=wp) :: density_rscale       = nan64
     real(kind=wp) :: distance2cm          = kpc2cm
     character(len=128) :: distance_unit   = 'kpc'
     character(len=128) :: source_geometry = 'point'
     character(len=128) :: base_name       = ''
     character(len=128) :: out_file        = ''
     integer            :: out_bitpix      = -32
     !--- location of a point source.
     real(kind=wp) :: xs_point = 0.0_wp
     real(kind=wp) :: ys_point = 0.0_wp
     real(kind=wp) :: zs_point = 0.0_wp
     !--- shared memory & master-slave algorithm
     integer       :: num_send_at_once   = 10000
     logical       :: use_shared_memory  = .false.
     logical       :: use_master_slave   = .true.
     logical       :: use_reduced_wgt    = .true.
     !--- Dust-related parameters
     real(kind=wp) :: hgg           = 0.6761
     real(kind=wp) :: albedo        = 0.3253
     real(kind=wp) :: cext_dust     = 1.6059e-21
     logical       :: use_stokes    = .true.
     character(len=128) :: scatt_mat_file = ''
     !--- external density file
     character(len=128) :: density_file  = ''
     integer            :: reduce_factor = 1
     integer            :: centering     = 0
     !--- parameters for PEELING-OFF
     real(kind=wp) :: distance          = nan64
     real(kind=wp) :: inclination_angle = nan64
     real(kind=wp) :: position_angle    = nan64
     real(kind=wp) :: phase_angle       = nan64
     real(kind=wp) :: alpha = nan64
     real(kind=wp) :: beta  = nan64
     real(kind=wp) :: gamma = nan64
     real(kind=wp) :: obsx  = nan64
     real(kind=wp) :: obsy  = nan64
     real(kind=wp) :: obsz  = nan64
     real(kind=wp) :: dxim  = nan64
     real(kind=wp) :: dyim  = nan64
     integer :: nxim = 129
     integer :: nyim = 129
     !--- number of scatterings, these are not input parameters.
     real(kind=wp) :: nscatt_tot  = 0.0_wp
     real(kind=wp) :: exetime     = 0.0_wp
  end type params_type

  !--- observer type
  type observer_type
     real(kind=wp) :: x,y,z
     real(kind=wp) :: inclination_angle = 0.0_wp
     real(kind=wp) :: position_angle    = 0.0_wp
     real(kind=wp) :: phase_angle       = 0.0_wp
     real(kind=wp) :: distance          = -999.9_wp
     real(kind=wp) :: rmatrix(3,3)
     integer       :: nxim,nyim
     real(kind=wp) :: dxim,dyim
     real(kind=wp) :: steradian_pix
     real(kind=wp), pointer :: tau(:,:)   => null()
     real(kind=wp), pointer :: scatt(:,:) => null()
     real(kind=wp), pointer :: direc(:,:) => null()
     real(kind=wp), pointer :: I(:,:)     => null()
     real(kind=wp), pointer :: Q(:,:)     => null()
     real(kind=wp), pointer :: U(:,:)     => null()
     real(kind=wp), pointer :: V(:,:)     => null()
  end type observer_type

! scattering related global variables
  type scattering_matrix_type
     integer :: nPDF
     real(kind=wp), pointer :: coss(:)      => null()
     real(kind=wp), pointer :: S11(:)       => null()
     real(kind=wp), pointer :: S12(:)       => null()
     real(kind=wp), pointer :: S33(:)       => null()
     real(kind=wp), pointer :: S34(:)       => null()
     real(kind=wp), pointer :: phase_PDF(:) => null()
     integer,       pointer :: alias(:)     => null()
  end type scattering_matrix_type

  type params_mpi
     !--- p_rank    = process (thread) rank
     !--- h_rank    = host rank (rank defined for each process of a node)
     !--- nproc     = number of total processes
     !--- hostcomm  = communicator for host
     !--- SAME_HRANK_COMM  = communicator for processes with the same h_rank
     !--- SAME_HRANK_NPROC = number of processes with the same h_rank
     integer :: p_rank
     integer :: h_rank
     integer :: nproc
     integer :: hostcomm
     integer :: SAME_HRANK_COMM
     integer :: SAME_HRANK_NPROC
  end type params_mpi

  ! globar parameters
  type(params_type) :: par
  type(params_mpi)  :: mpar
  type(observer_type) :: observer
  type(scattering_matrix_type) :: scatt_mat

  !------------------------------------------------------------
  ! define procedure pointers
  procedure(raytrace_tau), pointer :: raytrace_to_tau => null()
  abstract interface
     subroutine raytrace_tau(photon,grid,tau_in)
     import
     type(photon_type), intent(inout) :: photon
     type(grid_type),   intent(inout) :: grid
     real(kind=wp),     intent(in)    :: tau_in
     end subroutine raytrace_tau
  end interface

  procedure(raytrace_edge), pointer :: raytrace_to_edge => null()
  abstract interface
     subroutine raytrace_edge(photon,grid,tau)
     import
     type(photon_type), intent(in)  :: photon
     type(grid_type),   intent(in)  :: grid
     real(kind=wp),     intent(out) :: tau
   end subroutine raytrace_edge
  end interface

  procedure(peeling_direct), pointer :: peeling_direct_photon => null()
  abstract interface
     subroutine peeling_direct(photon,grid,mode)
     import
     type(photon_type), intent(in) :: photon
     type(grid_type),   intent(in) :: grid
     integer, intent(in), optional :: mode
   end subroutine peeling_direct
  end interface

  procedure(peeling_scattered), pointer :: peeling_scattered_photon => null()
  abstract interface
     subroutine peeling_scattered(photon,grid)
     import
     type(photon_type), intent(in) :: photon
     type(grid_type),   intent(in) :: grid
   end subroutine peeling_scattered
  end interface

  procedure(scatter), pointer :: scattering => null()
  abstract interface
     subroutine scatter(photon,grid)
     import
     type(photon_type), intent(inout) :: photon
     type(grid_type),   intent(inout) :: grid
     end subroutine scatter
  end interface

  procedure(run_sim), pointer :: run_simulation => null()
  abstract interface
     subroutine run_sim(grid)
     import
     type(grid_type),   intent(inout) :: grid
     end subroutine run_sim
  end interface

  procedure(write_out), pointer :: write_output => null()
  abstract interface
     subroutine write_out(filename,grid)
     import
     character(len=*),  intent(in) :: filename
     type(grid_type),   intent(in) :: grid
     end subroutine write_out
  end interface
  !------------------------------------------------------------
end module define
