module random_sersic
!--- 3-D deprojected Sersic profile sampler for a galaxy bulge (MoCafe v2.00).
!--- The deprojected 3-D cumulative luminosity of a 2-D Sersic profile
!--- I(R) = I_0 exp(-b (R/Re)^(1/m)) is computed by sersic_cumulative_3D
!--- (Abel transform; Seon et al.), ported from LaRT_v2.00/random_sersic.f90.
!--- The radius is drawn with the alias method (rand_alias_linear) instead of
!--- the inverse-CDF locate used in the reference: the bin probabilities are
!--- built from the cumulative and fed to random_alias_setup once (cached), then
!--- rand_alias_linear returns a continuous radius, matching the alias usage in
!--- external_radiation.f90 / scattering_car.f90.
  use define,  only : wp, dp, pi, halfpi
  use random,  only : rand_number, random_alias_setup, rand_alias_linear
  implicit none
  private
  public :: rand_sersic, sersic_cumulative_3D
contains
  !=============================
  ! Deprojected 3-D cumulative luminosity from a 2-D Sersic profile.
  !   m = Sersic index (>= 0.5); radius in units of Re; profile normalized later.
  subroutine sersic_cumulative_3D(m, radius, profile, rmax_in)
    real(wp), intent(in)    :: m
    real(wp), intent(inout) :: radius(:), profile(:)
    real(wp), optional, intent(in) :: rmax_in

    integer,  parameter :: nx   = 20001
    real(dp), parameter :: xmax = 10000
    real(dp), save, allocatable :: sersic_x(:), sersic_fx(:)
    logical,  save :: sersic_fx_saved = .false.
    real(dp), save :: sersic_dlnx
    !$OMP THREADPRIVATE(sersic_x,sersic_fx,sersic_fx_saved,sersic_dlnx)

    real(dp) :: m2, fr, sum1, sum2, norm, rmin, rmax, dlnr, r, b
    integer  :: nr, nr1, i, j
    real(dp), parameter :: bcoeff(4) = (/ 4d0/405d0, 46d0/25515d0, &
                                          131d0/1148175d0, -2194697d0/30690717750d0 /)

    !--- b(m): Ciotti & Bertin (1999, A&A, 352, 447)
    b = 0.0d0
    do i=2,1,-1
       b = (b + bcoeff(i))/m
    enddo
    b  = b + 2d0*m - 1d0/3d0
    m2 = 2d0*m

    nr   = size(radius)
    nr1  = nr - 1
    rmax = 2.5d0*((14.995674d0+4.0947738d0*m-0.052804581d0*m*m)/b)**m
    rmin = -0.27566682d0+0.21713972d0*m+0.037967891d0*m*m
    if (present(rmax_in)) then
       if (rmax_in > 0.0_wp) rmax = rmax_in
    endif
    if (rmin < 0.0d0) then
       rmin = rmax/1d4
    else
       rmin = 0.005d0*(rmin/b)**m
    endif
    dlnr = log(rmax/rmin)/(nr1-1)

    !--- x-integral of the Abel transform (cached, m-independent).
    if (.not.(sersic_fx_saved)) then
       sersic_dlnx = log(xmax)/dble(nx-1)
       if (.not. allocated(sersic_x))   allocate(sersic_x(nx))
       if (.not. allocated(sersic_fx))  allocate(sersic_fx(nx))
       do i=1,nx
          sersic_x(i) = exp((i-1)*sersic_dlnx)
          if (sersic_x(i) == 1.d0) then
             sersic_fx(i) = halfpi
          else
             sersic_fx(i) = -sqrt(1d0-1d0/sersic_x(i)**2) + sersic_x(i)*atan(1d0/sqrt(sersic_x(i)**2-1d0))
          endif
       enddo
       sersic_fx_saved = .true.
    endif

    radius(1)  = 0.0_wp
    profile(1) = 0.0_wp
    do j=1,nr1
       r    = rmin * exp((j-1)*dlnr)
       sum1 = gammp(m2+1d0, b*r**(1d0/m))
       sum2 = 0.0d0
       do i=1,nx
          fr = exp(-b*(r*sersic_x(i))**(1d0/m)) * (r*sersic_x(i))**(1d0/m)
          if (i == 1 .or. i == nx) then
             sum2 = sum2 + sersic_x(i)*fr*sersic_fx(i)*0.5d0
          else
             sum2 = sum2 + sersic_x(i)*fr*sersic_fx(i)
          endif
       enddo
       norm = (2d0/pi) * (b**(m2+1))/m * exp(-gammln(m2+1d0))
       sum2 = norm * r**2 * sum2*sersic_dlnx
       radius(j+1)  = r
       profile(j+1) = sum1 + sum2
    enddo
  end subroutine sersic_cumulative_3D

  !=============================
  ! Draw a spherical radius (in units of Re) from the 3-D Sersic profile with
  ! the alias method.  rmax is the truncation radius in Re units.
  function rand_sersic(m, rmax) result(rr)
    real(wp), intent(in) :: m, rmax
    real(wp) :: rr

    integer, parameter :: nr = 200
    real(wp), save :: rad(nr), prob(nr), palias(nr-1)
    integer,  save :: alias(nr-1)
    real(wp), save :: m_sav = -1.0_wp, rmax_sav = -1.0_wp
    logical,  save :: initialized = .false.
    !$OMP THREADPRIVATE(rad,prob,palias,alias,m_sav,rmax_sav,initialized)

    real(wp) :: cum(nr), psum
    integer  :: i

    if (.not.initialized .or. m /= m_sav .or. rmax /= rmax_sav) then
       call sersic_cumulative_3D(m, rad, cum, rmax_in=rmax)
       cum(:) = cum(:) / cum(nr)                     ! normalized cumulative
       !--- node PDF from the cumulative (finite differences)
       prob(1) = (cum(2)-cum(1))/(rad(2)-rad(1))
       do i=2,nr-1
          prob(i) = (cum(i+1)-cum(i-1))/(rad(i+1)-rad(i-1))
       enddo
       prob(nr) = (cum(nr)-cum(nr-1))/(rad(nr)-rad(nr-1))
       prob(:)  = max(prob(:), 0.0_wp)
       psum = 0.5_wp*sum((prob(2:nr)+prob(1:nr-1))*(rad(2:nr)-rad(1:nr-1)))
       prob(:) = prob(:)/psum
       !--- selection probability of each bin (trapezoid area), then alias table
       do i=1,nr-1
          palias(i) = 0.5_wp*(prob(i)+prob(i+1))*(rad(i+1)-rad(i))
       enddo
       palias(:) = palias(:)/sum(palias)
       call random_alias_setup(palias, alias)
       m_sav = m;  rmax_sav = rmax;  initialized = .true.
    endif

    rr = rand_alias_linear(palias, alias, rad, prob)
  end function rand_sersic

  !=============================
  ! Incomplete gamma P(a,x) and log-gamma (Numerical Recipes; via LaRT).
  function gammln(xx) result(g)
    real(wp), intent(in) :: xx
    real(dp) :: g, tmp, x, ser, y
    real(dp), parameter :: stp = 2.5066282746310005_dp
    real(dp), parameter :: coef(6) = (/76.18009172947146_dp, -86.50532032941677_dp, &
         24.01409824083091_dp, -1.231739572450155_dp, 0.1208650973866179e-2_dp, &
         -0.5395239384953e-5_dp/)
    integer :: i
    x = xx;  tmp = x+5.5_dp;  tmp = (x+0.5_dp)*log(tmp)-tmp
    ser = 1.000000000190015_dp;  y = x
    do i=1,size(coef)
       y = y+1.0_dp;  ser = ser+coef(i)/y
    enddo
    g = tmp + log(stp*ser/x)
  end function gammln

  function gser(a,x) result(g)
    real(wp), intent(in) :: a, x
    real(dp) :: g, ap, del, summ
    integer, parameter :: ITMAX=100
    real(dp), parameter :: EPS=epsilon(x)
    integer :: n
    if (x == 0.0_wp) then; g = 0.0_dp; return; endif
    ap = a;  summ = 1.0_dp/a;  del = summ
    do n=1,ITMAX
       ap = ap+1.0_dp;  del = del*x/ap;  summ = summ+del
       if (abs(del) < abs(summ)*EPS) exit
    enddo
    g = summ*exp(-x+a*log(x)-gammln(a))
  end function gser

  function gcf(a,x) result(g)
    real(wp), intent(in) :: a, x
    real(dp) :: g, an, b, c, d, del, h
    integer, parameter :: ITMAX=100
    real(dp), parameter :: EPS=epsilon(x), FPMIN=tiny(x)/EPS
    integer :: i
    if (x == 0.0_wp) then; g = 1.0_dp; return; endif
    b = x+1.0_dp-a;  c = 1.0_dp/FPMIN;  d = 1.0_dp/b;  h = d
    do i=1,ITMAX
       an = -i*(i-a);  b = b+2.0_dp;  d = an*d+b
       if (abs(d) < FPMIN) d = FPMIN
       c = b+an/c
       if (abs(c) < FPMIN) c = FPMIN
       d = 1.0_dp/d;  del = d*c;  h = h*del
       if (abs(del-1.0_dp) <= EPS) exit
    enddo
    g = exp(-x+a*log(x)-gammln(a))*h
  end function gcf

  function gammp(a,x) result(g)
    real(wp), intent(in) :: a, x
    real(dp) :: g
    if (x < a+1.0_dp) then
       g = gser(a,x)
    else
       g = 1.0_dp - gcf(a,x)
    endif
  end function gammp

end module random_sersic
