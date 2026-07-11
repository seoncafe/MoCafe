module random_bulge
!--- Boxy / bar / x-bar bulge samplers for a galaxy (MoCafe v2.00), ported from
!--- ~/MoCafe/Galaxy_v1.7/rand_bulge.f90.  Each returns a position (in units of
!--- the bulge scale) drawn from a generalized-ellipsoid density with boxiness
!--- exponent n: r = (|x|^n + |y|^n + |z|^n)^(1/n).
!---   boxy : power-law density  f ~ r^(-pow_idx)  (cored at rcen)
!---   bar  : Gaussian density   f ~ exp(-r^2)     (cored at rcen)
!---   xbar : peanut/x-shaped     f ~ (1 - (|x|^n+|y|^n-z^2))^2
!--- These are Gibbs samplers: each draw is conditioned on the previous one
!--- (module-saved state), so consecutive samples are correlated -- adequate in
!--- aggregate for a bulge but statistically inferior to the alias-sampled
!--- Sersic bulge (random_sersic).
  use define, only : wp
  use random, only : rand_number
  implicit none
  private
  public :: rand_boxy, rand_bar, rand_xbar
contains
  !---------------------------------------------------------------
  subroutine rand_boxy(x, y, z, rcen, xmax, ymax, zmax, n, pow_idx)
    real(wp), intent(in)  :: rcen, xmax, ymax, zmax, n, pow_idx
    real(wp), intent(out) :: x, y, z
    logical,  save :: initialized = .false.
    real(wp), save :: x_sav = 0.0_wp, y_sav = 0.0_wp, z_sav = 0.0_wp
    !$OMP THREADPRIVATE(initialized,x_sav,y_sav,z_sav)
    real(wp) :: xn, yn, zn, r, f, f0, fn, ninv

    if (.not. initialized) then
       x_sav = (2.0_wp*rand_number()-1.0_wp)*xmax
       y_sav = (2.0_wp*rand_number()-1.0_wp)*ymax
       z_sav = (2.0_wp*rand_number()-1.0_wp)*zmax
       initialized = .true.
    endif
    x = x_sav;  y = y_sav;  z = z_sav
    xn = abs(x)**n;  yn = abs(y)**n;  zn = abs(z)**n
    ninv = 1.0_wp/n
    r = (xn+yn+zn)**ninv
    if (r <= rcen) then
       f = rcen**(-pow_idx)
    else
       f = r**(-pow_idx)
    endif
    f0 = f*rand_number()
    fn = f0**(-n/pow_idx)
    x  = (fn-(yn+zn))**ninv;  x = (2.0_wp*rand_number()-1.0_wp)*minval([x,xmax]);  xn = abs(x)**n
    y  = (fn-(xn+zn))**ninv;  y = (2.0_wp*rand_number()-1.0_wp)*minval([y,ymax]);  yn = abs(y)**n
    z  = (fn-(xn+yn))**ninv;  z = (2.0_wp*rand_number()-1.0_wp)*minval([z,zmax])
    x_sav = x;  y_sav = y;  z_sav = z
  end subroutine rand_boxy

  !---------------------------------------------------------------
  subroutine rand_bar(x, y, z, rcen, xmax, ymax, zmax, n)
    real(wp), intent(in)  :: rcen, xmax, ymax, zmax, n
    real(wp), intent(out) :: x, y, z
    logical,  save :: initialized = .false.
    real(wp), save :: x_sav = 0.0_wp, y_sav = 0.0_wp, z_sav = 0.0_wp
    !$OMP THREADPRIVATE(initialized,x_sav,y_sav,z_sav)
    real(wp) :: xn, yn, zn, r, f, f0, fn, ninv

    if (.not. initialized) then
       x_sav = (2.0_wp*rand_number()-1.0_wp)*xmax
       y_sav = (2.0_wp*rand_number()-1.0_wp)*ymax
       z_sav = (2.0_wp*rand_number()-1.0_wp)*zmax
       initialized = .true.
    endif
    x = x_sav;  y = y_sav;  z = z_sav
    xn = abs(x)**n;  yn = abs(y)**n;  zn = abs(z)**n
    ninv = 1.0_wp/n
    r = (xn+yn+zn)**ninv
    if (r <= rcen) then
       f = exp(-rcen**2)
    else
       f = exp(-r**2)
    endif
    f0 = f*rand_number()
    fn = (-log(f0))**(n/2.0_wp)
    x  = (fn-(yn+zn))**ninv;  x = (2.0_wp*rand_number()-1.0_wp)*minval([x,xmax]);  xn = abs(x)**n
    y  = (fn-(xn+zn))**ninv;  y = (2.0_wp*rand_number()-1.0_wp)*minval([y,ymax]);  yn = abs(y)**n
    z  = (fn-(xn+yn))**ninv;  z = (2.0_wp*rand_number()-1.0_wp)*minval([z,zmax])
    x_sav = x;  y_sav = y;  z_sav = z
  end subroutine rand_bar

  !---------------------------------------------------------------
  subroutine rand_xbar(x, y, z, xmax, ymax, zmax, n)
    real(wp), intent(in)  :: xmax, ymax, zmax, n
    real(wp), intent(out) :: x, y, z
    logical,  save :: initialized = .false.
    real(wp), save :: x_sav = 0.0_wp, y_sav = 0.0_wp, z_sav = 0.0_wp
    !$OMP THREADPRIVATE(initialized,x_sav,y_sav,z_sav)
    real(wp) :: xn, yn, z2, f, f0, f1, ninv

    if (.not. initialized) then
       x_sav = (2.0_wp*rand_number()-1.0_wp)*xmax
       y_sav = (2.0_wp*rand_number()-1.0_wp)*ymax
       z_sav = (2.0_wp*rand_number()-1.0_wp)*zmax
       initialized = .true.
    endif
    x = x_sav;  y = y_sav;  z = z_sav
    ninv = 1.0_wp/n
    xn = abs(x)**n;  yn = abs(y)**n;  z2 = abs(z)**2
    f  = (1.0_wp-(xn+yn-z2))**2
    f0 = f*rand_number()
    f1 = 1.0_wp-sqrt(f0)
    x  = (f1-(yn-z2))**ninv;  x = (2.0_wp*rand_number()-1.0_wp)*minval([x,xmax]);  xn = abs(x)**n
    y  = (f1-(xn-z2))**ninv;  y = (2.0_wp*rand_number()-1.0_wp)*minval([y,ymax]);  yn = abs(y)**n
    z  = sqrt(max(xn+yn-f1, 0.0_wp));  z = (2.0_wp*rand_number()-1.0_wp)*minval([z,zmax])
    x_sav = x;  y_sav = y;  z_sav = z
  end subroutine rand_xbar
end module random_bulge
