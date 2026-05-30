  implicit none
  integer, parameter :: niso=1000
  real  :: x(niso),y(niso)
  integer :: i
  real, parameter :: ximax = 100

  call isothermal(x,y,ximax)
  do i=1,niso
    print*, x(i), y(i)
  enddo

  stop
  end
!---------------------------------------------
! isothermal gas sphere - Bonnor-Ebert sphere
!
! Kwang-Il Seon, KASI, 2012
!------------
! how to use:
! call isothermal(x,y,ximax)
!   Then the subroutine gives x(1:niso) and y(1:niso).
!   x array is radius normalized to lscale defined by the following equations.
!   y array is density normalized to rho_c (central density of cloud)
!------------
! ximax   = 7.0, in most cases
! radius  = x * lscale
!           lscale = sqrt(c_s^2/(4*pi*G*rho_c))
!           c_s    = 7.67855e3 (T/K) cm s^-1
!           G      = 6.67384e-8 cm^3 g^-1 s^-2
!           rho_c  = central density
! density = y * rho_c
!---------------------------------------------
  subroutine isothermal(x,y,ximax)

  use define, only : wp
  implicit none
  integer,parameter :: niso=1000
  real(kind=wp) :: x(niso), u(niso), y(niso), ximax

  integer :: i
  real(kind=wp) :: dx, xmax

  xmax   = ximax
  dx     = xmax/(dble(niso)-1)
  do i=1,niso
    x(i) = dx*(dble(i)-1)
  enddo

  u(1:2) = 0.0
  do i=3,niso
    u(i) = 2.0*u(i-1)-u(i-2)+dx**2*(-2.0/(x(i-1)*dx)*(u(i-1)-u(i-2)) + exp(-u(i-1)) )
  enddo
  y(:) = exp(-u(:))

  end
