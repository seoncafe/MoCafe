module define
public

! kind parameter for single and double precision (6 and 15 significant decimal digits)
! precision of 32-, 64-, and 128-bit reals
  integer, parameter :: sp = selected_real_kind(6, 37)
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: qp = selected_real_kind(33, 4931)
! precision of 32-, 64- bit integers
  integer, parameter :: i4b = selected_int_kind(r=9)
  integer, parameter :: i8b = selected_int_kind(r=18)

! working precision
! Note that (1) nphotons is always 64-bit integer.
!           (2) other integer variables are always 32-bit integers.
!  integer, parameter :: wp  = sp
!  integer, parameter :: iwp = i4b
  integer, parameter :: wp  = dp
  integer, parameter :: iwp = i8b

! tinest = the smallest positive number
! eps    = the least positive number
!          that added to 1 returns a number that is greater than 1
! hugest = the hugest positive number
  real(kind=wp), parameter :: tinest = tiny(0.0_wp)
  real(kind=wp), parameter :: eps    = epsilon(0.0_wp)
  real(kind=sp), parameter :: eps_sp = epsilon(0.0_sp)
  real(kind=dp), parameter :: eps_dp = epsilon(0.0_dp)
  real(kind=wp), parameter :: hugest = huge(1.0_wp)

! numerical values
  real(kind=wp), parameter :: pi     = 3.141592653589793238462643383279502884197_wp
  real(kind=wp), parameter :: twopi  = 6.283185307179586476925286766559005768394_wp
  real(kind=wp), parameter :: fourpi = 12.56637061435917295385057353311801153679_wp
  real(kind=wp), parameter :: halfpi = pi/2.0_wp

! conversion factor (radian to degree, degree to radian)
  real(kind=wp), parameter :: rad2deg = 180.0_wp/pi
  real(kind=wp), parameter :: deg2rad = pi/180.0_wp

! distances
  real(kind=wp), parameter :: pc2cm  = 3.0856776e18_wp
  real(kind=wp), parameter :: kpc2cm = pc2cm * 1e3_wp
  real(kind=wp), parameter :: au2cm  = 1.4960e13_wp
  real(kind=wp) :: distance2cm

! dust mass per hydrogen nuclei (g/H)
  real(kind=wp), parameter :: h2dust = 1.87e-26_wp

! kappa_v : dust extinction cross-section (cm^2/g) at V-band
  real(kind=wp), parameter :: kappa_v = 25992.6_wp

! photon type
  type photon_type
     real(kind=wp) :: x
     real(kind=wp) :: y
     real(kind=wp) :: z
     real(kind=wp) :: vx
     real(kind=wp) :: vy
     real(kind=wp) :: vz
     integer :: rcell
     integer :: pcell
     integer :: zcell
     real(kind=wp) :: wgt
     real(kind=wp) :: lscale
  end type photon_type

! cylindrical grid type
! nr, np, nz          - number of cells
! rface, pface, zface - locations of cell faces
! opacity             - opacity per unit length
  type grid_type
     integer :: nr
     integer :: np
     integer :: nz
     real(kind=wp) :: rmin,rmax,rrange
     real(kind=wp) :: pmin,pmax,prange
     real(kind=wp) :: zmin,zmax,zrange
     real(kind=wp) :: r_alpha
     real(kind=wp) :: z_alpha
     real(kind=wp) :: zmin_eq,zmax_eq,dz_eq
     real(kind=wp) :: rmin_eq,rmax_eq,dr_eq
     real(kind=wp) :: dp
     real(kind=wp), pointer :: dr(:)
     real(kind=wp), pointer :: dz(:)
     real(kind=wp), pointer :: rface(:)
     real(kind=wp), pointer :: pface(:)
     real(kind=wp), pointer :: zface(:)
     real(kind=wp), pointer :: sinp(:)
     real(kind=wp), pointer :: cosp(:)
     real(kind=wp), pointer :: opacity(:,:,:)
  end type grid_type

! observer type
  type observer_type
     sequence
     real(kind=wp) :: x
     real(kind=wp) :: y
     real(kind=wp) :: z
     real(kind=wp) :: rmatrix(3,3)
  end type observer_type

! output type
  type output_type
     integer :: nx
     integer :: ny
     real(kind=wp) :: dx
     real(kind=wp) :: dy
     real(kind=wp) :: steradian_pix
     real(kind=wp), pointer :: scatt(:,:)
     real(kind=wp), pointer :: direc(:,:)
     real(kind=wp), pointer :: tot(:,:)
     real(kind=wp), pointer :: scatt_sig(:,:)
     real(kind=wp), pointer :: direc_sig(:,:)
     real(kind=wp), pointer :: tot_sig(:,:)
  end type output_type

! dust distribution type
  type dust_type
     character(len=128) :: name  = ''
     real(kind=wp) :: tau_faceon = 0.0_wp    ! central face-on optical depth
     real(kind=wp) :: rscale     = 6.0_wp
     real(kind=wp) :: zscale     = 0.2_wp
     real(kind=wp) :: rmax       = 18.0_wp
     real(kind=wp) :: zmax       = 18.0_wp
     real(kind=wp) :: rcenter    = 6.0_wp    ! ring-like geometry
  end type dust_type

! source distribution type
  type source_type
     character(len=128) :: diskname = 'exponential'
     real(kind=wp) :: rscale = 6.0          ! disk-type
     real(kind=wp) :: zscale = 0.2          ! disk-type
     real(kind=wp) :: rmax   = 18.0         ! disk-type
     real(kind=wp) :: zmax   = 18.0         ! disk-type
     character(len=128) :: bulgename = 'sersic'
     real(kind=wp) :: sersic_index = 4.0    ! bulge-type
     real(kind=wp) :: Reff         = 2.3    ! bulge-type
     real(kind=wp) :: axial_ratio  = 0.5    ! bulge-type
     real(kind=wp) :: BulgeToDisk  = 0.0    ! bulge-type
  end type source_type

!- input parameters
  integer, parameter :: MAX_dust_geometry   = 10
  type params_type
     integer(kind=i8b) :: nphotons  = 1e6
     integer(kind=i8b) :: nprint    = 1e7
     real(kind=wp) :: hgg           = 0.5
     real(kind=wp) :: albedo        = 0.5
     real(kind=wp) :: luminosity    = 1.0
     real(kind=wp) :: tau_crit      = 2.0_wp
     real(kind=wp) :: dist_crit     = 1.0_wp
     real(kind=wp) :: wgt_min       = 1e-6_wp
     real(kind=wp) :: p_survival    = 0.5_wp
     real(kind=wp) :: inclination_angle = 90.0_wp
     real(kind=wp) :: position_angle    = 0.0_wp
     real(kind=wp) :: phase_angle       = 0.0_wp
     real(kind=wp) :: distance          = 9500.0_wp
     real(kind=wp) :: xshift = 0.0_wp
     real(kind=wp) :: yshift = 0.0_wp
     integer       :: nr   = 40
     integer       :: np   = 1
     integer       :: nz   = 350
     real(kind=wp) :: rmax = 18.0
     real(kind=wp) :: zmax = 18.0
     real(kind=wp) :: r_alpha = -999.0_wp
     real(kind=wp) :: z_alpha = -999.0_wp
     integer       :: nxim = 100
     integer       :: nyim = 100
     real(kind=wp) :: dxim = 0.000416667
     real(kind=wp) :: dyim = 0.000416667
     character(len=128) :: distance_unit = 'kpc'
     character(len=128) :: out_file      = 'out.fits.gz'
     character(len=128) :: psf_file      = ''
     integer            :: ndust   = 1
     type(dust_type)    :: dust(MAX_dust_geometry)
     type(source_type)  :: source
     logical            :: left_right_fold = .true.
     integer            :: nthreads    = 8
     integer            :: output_mode = 1
     integer            :: iseed = 0
     !--- Let's keep these here for now.
     real(kind=wp) :: gg1, gg2, gg3
  end type params_type

!------------------------
! Global Parameters
  type(params_type) :: par
end module define
