module mathlib
  use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64
  implicit none
  private sp,dp,wp, PI, TWOPI
  integer, parameter  :: sp = real32
  integer, parameter  :: dp = real64
  integer, parameter  :: wp = dp
  real(wp), parameter :: PI = 3.141592653589793238462643383279502884197_wp
  real(wp), parameter :: TWOPI = PI * 2.0_wp
contains
  function arctan(y,x) result(angle)
     real(kind=wp), intent(in) :: x,y
     real(kind=wp) :: angle
     angle = atan2(y,x)
     if (angle < 0.0_wp) angle = angle + twopi
  end function arctan
!=====================================================================
  subroutine interp(nx,x,y,xnew,ynew)
  !use define, only : wp
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
  end subroutine interp
!=====================================================================
  SUBROUTINE locate(xx,n,x,j)
  !use define, only : wp
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
  END SUBROUTINE locate

!--------
!  These routines are from Numerical Recipes.
!    2015-11-16
!-------- Gammln
  FUNCTION gammln(xx)
  IMPLICIT NONE
  REAL(WP), INTENT(IN) :: xx
  REAL(DP) :: gammln
  REAL(DP) :: tmp,x
  REAL(DP) :: ser,y
  REAL(DP) :: stp = 2.5066282746310005_dp
  REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
                                    -86.50532032941677_dp,&
                                     24.01409824083091_dp,&
                                     -1.231739572450155_dp,&
                                      0.1208650973866179e-2_dp,&
                                     -0.5395239384953e-5_dp/)
  INTEGER :: i
  x   = xx
  tmp = x+5.5_dp
  tmp = (x+0.5_dp)*log(tmp)-tmp
  ser = 1.000000000190015_dp
  y   = x
  do i=1, size(coef)
     y   = y+1.0_dp
     ser = ser+coef(i)/y
  end do
  gammln = tmp + log(stp*ser/x)
  END FUNCTION gammln
!-------- Gcf
  FUNCTION gcf(a,x,gln)
  IMPLICIT NONE
  REAL(WP), INTENT(IN) :: a,x
  REAL(WP), OPTIONAL, INTENT(OUT) :: gln
  REAL(DP) :: gcf
  INTEGER,  PARAMETER :: ITMAX=100
  REAL(DP), PARAMETER :: EPS=epsilon(x), FPMIN=tiny(x)/EPS
  INTEGER :: i
  REAL(DP) :: an,b,c,d,del,h
  if (x == 0.0) then
     gcf = 1.0
     RETURN
  end if
  b = x+1.0_dp-a
  c = 1.0_dp/FPMIN
  d = 1.0_dp/b
  h = d
  do i=1,ITMAX
     an = -i*(i-a)
     b  = b+2.0_dp
     d  = an*d+b
     if (abs(d) < FPMIN) d = FPMIN
     c  = b+an/c
     if (abs(c) < FPMIN) c = FPMIN
     d   = 1.0_dp/d
     del = d*c
     h   = h*del
     if (abs(del-1.0_dp) <= EPS) exit
  end do
  if (i > ITMAX) then
     print*,'a too large, ITMAX too small in gcf'
     stop
  endif
  if (present(gln)) then
     gln = gammln(a)
     gcf = exp(-x+a*log(x)-gln)*h
  else
     gcf = exp(-x+a*log(x)-gammln(a))*h
  end if
  END FUNCTION gcf

!-------- Gser
  FUNCTION gser(a,x,gln)
  IMPLICIT NONE
  REAL(WP), INTENT(IN) :: a,x
  REAL(WP), OPTIONAL, INTENT(OUT) :: gln
  REAL(DP) :: gser
  INTEGER,  PARAMETER :: ITMAX=100
  REAL(DP), PARAMETER :: EPS=epsilon(x)
  INTEGER  :: n
  REAL(DP) :: ap,del,summ
!
  if (x == 0.0) then
     gser = 0.0
     RETURN
  end if
  ap   = a
  summ = 1.0_dp/a
  del  = summ
  do n=1, ITMAX
     ap   = ap+1.0_dp
     del  = del*x/ap
     summ = summ+del
     if (abs(del) < abs(summ)*EPS) exit
  end do
  if (n > ITMAX) then
     print*,'a too large, ITMAX too small in gser'
     stop
  endif
  if (present(gln)) then
     gln  = gammln(a)
     gser = summ*exp(-x+a*log(x)-gln)
  else
     gser = summ*exp(-x+a*log(x)-gammln(a))
  end if
  END FUNCTION gser

!-------- GammaP
  FUNCTION gammp(a,x)
  IMPLICIT NONE
  REAL(WP), INTENT(IN) :: a,x
  REAL(DP) :: gammp
  if (x < a+1.0_dp) then
     gammp = gser(a,x)
  else
     gammp = 1.0_dp-gcf(a,x)
  end if
  END FUNCTION gammp

!-------- GammaQ
  FUNCTION gammq(a,x)
  IMPLICIT NONE
  REAL(WP), INTENT(IN) :: a,x
  REAL(DP) :: gammq
  if (x < a+1.0_dp) then
     gammq = 1.0_dp-gser(a,x)
  else
     gammq = gcf(a,x)
  end if
  END FUNCTION gammq

  FUNCTION erfc(x)
  implicit none
  real(wp), intent(in) :: x
  real(dp) :: erfc
  erfc = merge(1.0_dp+gammp(0.5_dp,x**2),gammq(0.5_dp,x**2), x < 0.0)
  END FUNCTION erfc

!-------- convolution_2D
   subroutine convolution_2D(img,psf,output)
!--- Convolution of two 2D images "img" and "psf" using FFTW (www.fftw.org).
!--- 2015/11/21, Kwang-Il Seon
!--- 2016/01/30, Bug-fixed, center locations for FFTW were wrong.
   implicit none
   real(kind=wp), intent(inout) :: img(:,:)
   real(kind=wp), intent(in)    :: psf(:,:)
   !real(kind=wp), allocatable, intent(in)  :: img(:,:),psf(:,:)
   real(kind=wp), intent(out), optional :: output(:,:)

   ! local varables
   real(kind=wp),    allocatable :: img_tmp(:,:), psf_tmp(:,:), con_tmp(:,:)
   complex(kind=wp), allocatable :: img_fft(:,:), psf_fft(:,:), con_fft(:,:)
   ! plan should always be 8-byte integer
   integer(8) :: plan
   integer    :: nx,ny,nx_img,ny_img,nx_psf,ny_psf
   integer    :: i1,i2,j1,j2,dx,dy,ia,ib,ja,jb
   integer, parameter :: FFTW_ESTIMATE = 64
   logical :: watch_out

   nx_img = size(img,1)
   ny_img = size(img,2)
   nx_psf = size(psf,1)
   ny_psf = size(psf,2)

   if (present(output)) then
      if (size(output,1) /= nx_img .or. size(output,2) /= ny_img) then
         write(*,*) 'size of optional output should be the same as the input img.'
         write(*,*) 'input img will be overwritten.'
         watch_out = .true.
      else
         watch_out = .false.
      endif
   endif

   nx = nx_img
   ny = ny_img
   if (nx < nx_psf) nx = nx_psf
   if (ny < ny_psf) ny = ny_psf

   if (.not. allocated(img_tmp)) allocate(img_tmp(nx,ny))
   if (.not. allocated(psf_tmp)) allocate(psf_tmp(nx,ny))
   if (.not. allocated(con_tmp)) allocate(con_tmp(nx,ny))

   !--- Zero Padding
   !---    do not change the order. PSF first and then IMG next.
   !---    i1, i2 for img will be used later.
   dx = (nx/2+1) - (nx_psf/2+1)
   dy = (ny/2+1) - (ny_psf/2+1)
   i1 = 1      + dx
   i2 = nx_psf + dx
   j1 = 1      + dy
   j2 = ny_psf + dy
   psf_tmp = 0.0
   !--- Normalize PSF
   !psf_tmp(i1:i2,j1:j2) = psf(:,:)/sum(psf)
   where (psf > 0.0_wp)
      psf_tmp(i1:i2,j1:j2) = psf(:,:)/sum(psf)
   endwhere

   dx = (nx/2+1) - (nx_img/2+1)
   dy = (ny/2+1) - (ny_img/2+1)
   i1 = 1      + dx
   i2 = nx_img + dx
   j1 = 1      + dy
   j2 = ny_img + dy
   img_tmp = 0.0
   !img_tmp(i1:i2,j1:j2) = img(:,:)
   ! dividing by sqrt(nx*ny) twice rather than dividing by nx*ny might be better to prevent numerical overflow.
   img_tmp(i1:i2,j1:j2) = img(:,:)/sqrt(real(nx*ny))

   if (.not. allocated(img_fft)) allocate(img_fft(nx/2+1,ny))
   if (.not. allocated(psf_fft)) allocate(psf_fft(nx/2+1,ny))
   if (.not. allocated(con_fft)) allocate(con_fft(nx/2+1,ny))

   ! FFTW manual says:
   ! You should normally pass the same input/output arrays that were used when creating the plan. This is always safe.
   if (wp == sp) then
      call sfftw_plan_dft_r2c_2d(plan,nx,ny,img_tmp,img_fft,FFTW_ESTIMATE)
      call sfftw_execute_dft_r2c(plan,img_tmp,img_fft)
      call sfftw_destroy_plan(plan)

      call sfftw_plan_dft_r2c_2d(plan,nx,ny,psf_tmp,psf_fft,FFTW_ESTIMATE)
      call sfftw_execute_dft_r2c(plan,psf_tmp,psf_fft)
      call sfftw_destroy_plan(plan)
   else
      call dfftw_plan_dft_r2c_2d(plan,nx,ny,img_tmp,img_fft,FFTW_ESTIMATE)
      call dfftw_execute_dft_r2c(plan,img_tmp,img_fft)
      call dfftw_destroy_plan(plan)

      call dfftw_plan_dft_r2c_2d(plan,nx,ny,psf_tmp,psf_fft,FFTW_ESTIMATE)
      call dfftw_execute_dft_r2c(plan,psf_tmp,psf_fft)
      call dfftw_destroy_plan(plan)
   endif

   ! Convolution in k-space
   con_fft = psf_fft * img_fft

   if (wp == sp) then
      ! bug-fixed 2015-12-06
      call sfftw_plan_dft_c2r_2d(plan,nx,ny,con_fft,con_tmp,FFTW_ESTIMATE)
      call sfftw_execute_dft_c2r(plan,con_fft,con_tmp)
      call sfftw_destroy_plan(plan)
   else
      ! bug-fixed 2015-12-06
      call dfftw_plan_dft_c2r_2d(plan,nx,ny,con_fft,con_tmp,FFTW_ESTIMATE)
      call dfftw_execute_dft_c2r(plan,con_fft,con_tmp)
      call dfftw_destroy_plan(plan)
   endif

   ! FFTW manual says:
   !   FFTW computes an unnormalized DFT.
   !   Thus, computing a forward followed by a backward transform (or vice versa) results in the original array scaled by n.
   !   normalize PSF...
   !con_tmp = con_tmp/real(nx*ny)
   con_tmp = con_tmp/sqrt(real(nx*ny))

   ! Swap the quadrants
   !   center for FFTW is defined as nx/2+1, ny/2+1.
   ib = nx/2    + 1
   ia = nx - ib + 1
   jb = ny/2    + 1
   ja = ny - jb + 1
   img_tmp(1:ia,:)    = con_tmp(ib:nx,:)
   img_tmp(ia+1:nx,:) = con_tmp(1:ib-1,:)
   con_tmp(:,1:ja)    = img_tmp(:,jb:ny)
   con_tmp(:,ja+1:ny) = img_tmp(:,1:jb-1)

   ! Return convolved image with the size of the input IMG.
   if (present(output)) then
      !if (.not. allocated(output)) allocate(output(nx_img,ny_img))
      if (.not.(watch_out)) then
         output(:,:) = con_tmp(i1:i2,j1:j2)
      else
         img(:,:)    = con_tmp(i1:i2,j1:j2)
      endif
   else
      img(:,:)       = con_tmp(i1:i2,j1:j2)
   endif

   if (allocated(img_tmp)) deallocate(img_tmp)
   if (allocated(psf_tmp)) deallocate(psf_tmp)
   if (allocated(con_tmp)) deallocate(con_tmp)
   if (allocated(img_fft)) deallocate(img_fft)
   if (allocated(psf_fft)) deallocate(psf_fft)
   if (allocated(con_fft)) deallocate(con_fft)

   end subroutine convolution_2D
!------------------------------------------------
   subroutine loop_divide(nloop,nthreads,myid,mystart,myend)
     implicit none
     integer, intent(in)  :: nloop,  nthreads, myid
     integer, intent(out) :: mystart, myend
     integer :: remainder

     remainder = mod(nloop, nthreads)
     mystart   = myid * (nloop/nthreads) + 1
     if (remainder > myid) then
        mystart = mystart + myid
        myend   = mystart + (nloop/nthreads)
     else
        mystart = mystart + remainder
        myend   = mystart + (nloop/nthreads) - 1
     endif
   end subroutine loop_divide

end module mathlib
