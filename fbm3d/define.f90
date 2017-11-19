module define
public
! kind parameter for single and double precision (6 and 15 significant decimal digits)
! precision of 32-, 64-, and 128-bit reals
  integer, parameter :: sp = selected_real_kind(6, 37)
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: qp = selected_real_kind(33, 4931)
  integer, parameter :: wp = sp

  integer, parameter :: i4b = selected_int_kind(r=9)
  integer, parameter :: i8b = selected_int_kind(r=18)

! numerical values
  real(kind=wp), parameter :: pi     = 3.141592653589793238462643383279502884197_wp
  real(kind=wp), parameter :: twopi  = 6.283185307179586476925286766559005768394_wp
  real(kind=wp), parameter :: fourpi = 12.56637061435917295385057353311801153679_wp

! conversion factor
  real(kind=wp), parameter :: rad2deg = 180.0_wp/pi
  real(kind=wp), parameter :: deg2rad = pi/180.0_wp

! tinest = the smallest positive number
! eps    = the least positive number
!          that added to 1 returns a number that is greater than 1
  real(kind=wp), parameter :: tinest = tiny(0.0_wp)
  real(kind=wp), parameter :: eps    = epsilon(0.0_wp)
  real(kind=wp), parameter :: hugest = huge(1.0_wp)

end module define
