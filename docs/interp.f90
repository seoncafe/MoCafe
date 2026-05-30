  subroutine interp(nx,x,y,xnew,ynew)
  use define, only : wp
  implicit none
  integer, intent(in) :: nx
  real(kind=wp), intent(in)  :: x(nx),y(nx),xnew
  real(kind=wp), intent(out) :: ynew
!---------------
! local variable
  integer :: i
!
  call locate(x,nx,xnew,i)
  if (i == 0) then
     ynew = y(1)
  else if (i == nx) then
     ynew = y(nx)
  else
     ynew = y(i) + (y(i+1)-y(i))*(xnew-x(i))/(x(i+1)-x(i))
  endif
  return
  end
!=====================================================================
  SUBROUTINE locate(xx,n,x,j)
  use define, only : wp
  IMPLICIT NONE
  INTEGER :: j,n
  REAL(kind=wp) :: x,xx(n)
! Given an array xx(1:n), and given a value x, returns a value j such that
! x is between xx(j) and xx(j+1). xx(1:n) must be monotonic, either increasing
! or decreasing. j=0 or j=n is returned to indicate that x is out of range.
  INTEGER :: jl,jm,ju
   jl=0   !Initialize lower
   ju=n+1 !and upper limits.
10 if(ju-jl.gt.1)then !If we are not yet done,
     jm=(ju+jl)/2     !compute a midpoint,
     if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
       jl=jm          !and replace either the lower limit
     else
       ju=jm          !or the upper limit, as appropriate.
     endif
   goto 10            !Repeat until
   endif              !the test condition 10 is satised.
   if(x.eq.xx(1))then !Then set the output
     j=1
   else if(x.eq.xx(n))then
     j=n-1
   else
     j=jl
   endif
   return
  END
