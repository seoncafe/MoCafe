module scatter
use define, only : wp
use random, only : rand_number
implicit none
public

contains
!==================================================================
function rand_hg(g) result(hgdev)
!-----------------------------------------------------------
! Randoum Number Generator for Henyey-Greenstein function
! g  = asymmetry phase factor
!-----------------------------------------------------------
  implicit none
  real(kind=wp), intent(in) :: g
  real(kind=wp) :: hgdev
  ! local variables
  real(kind=wp) :: x,g2

  x = rand_number()
  if (g == 0.0_wp) then
     hgdev = 2.0_wp*x-1.0_wp
  else
     g2    = g*g
     hgdev = ((1.0_wp+g2)-((1.0_wp-g2)/(1.0_wp-g+2.0_wp*g*x))**2)/(2.0_wp*g)
  endif
  return
end function rand_hg

!-----------------------------------------------------------
function rand_draine(g,alpha) result(draine)
!-----------------------------------------------------------
! Randoum Number Generator for Draine's Phase Function
! Henyey-Greenstein function is obatined for alpha = 0.
! g  = asymmetry phase factor
! 08/09/2006 Kwang-il Seon
!-----------------------------------------------------------
  implicit none
  real(kind=wp), intent(in) :: g, alpha
  real(kind=wp) :: draine
  ! local variables
  real(kind=wp) :: x,y,g2

  do while(.true.)
     x = rand_number()
     if (g == 0.0_wp) then
        draine = 2.0_wp*x-1.0_wp
     else
        g2 = g*g
        draine = ((1.0_wp+g2)-((1.0_wp-g2)/(1.0_wp-g+2.0_wp*g*x))**2)/(2.0_wp*g)
     endif
     if (alpha /= 0.0_wp) then
        x = rand_number()
        y = (1.0_wp+alpha*draine**2)/(2.0_wp+alpha)
        if (x <= y) exit
     endif
  enddo
  return
end function rand_draine

subroutine random_sph(radius,x,y,z)
!-----------------------------------------------------
! random locations in the sphere with a radius r.
!-----------------------------------------------------
  implicit none
  real(kind=wp), intent(in)  :: radius
  real(kind=wp), intent(out) :: x,y,z

! local variables
  real(kind=wp) :: r,cost,sint,phi
  real(kind=wp), parameter :: twopi = 6.283185307179586476925286766559005768394_wp

  if (radius > 0.0_wp) then
    r    = radius*(rand_number())**(1.0_wp/3.0_wp)
    cost = 2.0_wp*rand_number()-1.0_wp
    sint = sqrt(1.0_wp-cost*cost)
    phi  = twopi*rand_number()
    x    = r*sint*cos(phi)
    y    = r*sint*sin(phi)
    z    = r*cost
  else
    x    = 0.0_wp
    y    = 0.0_wp
    z    = 0.0_wp
  endif
  return
end subroutine random_sph

subroutine scatt(photon)
  use define
  use random
  implicit none
!
! Calculate new direction cosines after scattering
! Written by Kwang-il Seon
!
  type(photon_type), intent(inout) :: photon

! local variables
  real(kind=wp) :: cost,sint,phi,cosp,sinp
  real(kind=wp) :: vx,vy,vz
  real(kind=wp) :: vx1,vy1,vz1,r

! New cos(theta)
  cost = rand_hg(par%hgg)
  sint = 1.0_wp - cost**2
  if (sint > 0.0_wp) then
    sint = sqrt(sint)
  else
    sint = 0.0_wp
  endif

! New phi
  phi  = twopi * rand_number()
  cosp = cos(phi)
  sinp = sin(phi)

! New direction
  vx1 = photon%vx
  vy1 = photon%vy
  vz1 = photon%vz
  if (abs(vz1) >= 0.99999999999_wp) then
     vx = sint*cosp
     vy = sint*sinp
     vz = cost
  else
     r = sqrt(1.0_wp - vz1*vz1)
     vx = cost*vx1 + sint*(vz1*vx1*sinp + vy1*cosp)/r
     vy = cost*vy1 + sint*(vz1*vy1*sinp - vx1*cosp)/r
     vz = cost*vz1 - sint*sinp*r
  endif
  r  = sqrt(vx*vx + vy*vy + vz*vz)
  photon%vx = vx/r
  photon%vy = vy/r
  photon%vz = vz/r

  return
end subroutine scatt

end module scatter
