module spectrum_sampler_mod
!--- Continuous wavelength sampling from a radiation spectrum.
!---
!--- A sampler holds a non-negative luminosity density L_lambda represented as a
!--- piecewise-linear function of wavelength on an ascending node set, together
!--- with its analytically integrated (trapezoidal) cumulative distribution.
!--- Inverting that cumulative distribution is exact for the stored density, so a
!--- photon receives a continuous wavelength instead of the center of a
!--- wavelength bin: the spectral structure inside a bin, and the variation of
!--- the dust cross sections across it, both survive.
!---
!--- Node sets always contain the transfer grid's bin edges, so the density is
!--- resolved at least as finely as that grid and the luminosity of every bin is
!--- an exact difference of the cumulative distribution rather than a midpoint
!--- estimate.
!---
!--- A blackbody spectrum is refined adaptively to a relative tolerance before it
!--- enters the same representation.
  use define, only : wp
  implicit none
  private

  integer,       parameter :: MAX_PLANCK_NODES = 200000
  integer,       parameter :: MAX_REFINE_DEPTH = 50
  real(kind=wp), parameter :: hc_over_k = 1.43877687750393e4_wp  ! [um K]

  type, public :: spectrum_sampler_type
     real(kind=wp), allocatable :: lam(:)    ! ascending nodes [um]
     real(kind=wp), allocatable :: dens(:)   ! L_lambda at the nodes (>= 0)
     real(kind=wp), allocatable :: cdf(:)    ! cumulative luminosity, normalized
     real(kind=wp) :: total = 0.0_wp         ! integral of L_lambda over the support
  end type spectrum_sampler_type

  public :: tabulated_spectrum_sampler, planck_spectrum_sampler
  public :: sample_wavelength, spectrum_cdf, band_luminosity_fraction
  public :: planck_lambda_shape, sampler_is_ready
  public :: sample_power_law_bin, logarithmic_mean

contains
  !---------------------------------------------------------------
  !--- Planck function B_lambda (arbitrary normalization), lambda in um.
  pure function planck_lambda_shape(lam_um, T) result(b)
  implicit none
  real(kind=wp), intent(in) :: lam_um, T
  real(kind=wp) :: b, x
  if (.not. (lam_um > 0.0_wp) .or. .not. (T > 0.0_wp)) then
     b = 0.0_wp;  return
  endif
  x = hc_over_k/(lam_um*T)
  if (x > 700.0_wp) then
     b = 0.0_wp
  else
     b = 1.0_wp/(lam_um**5 * (exp(x) - 1.0_wp))
  endif
  end function planck_lambda_shape

  !---------------------------------------------------------------
  pure function sampler_is_ready(sampler) result(ready)
  implicit none
  type(spectrum_sampler_type), intent(in) :: sampler
  logical :: ready
  ready = allocated(sampler%cdf)
  if (ready) ready = size(sampler%cdf) >= 2 .and. sampler%total > 0.0_wp
  end function sampler_is_ready

  !---------------------------------------------------------------
  !--- Sampler for a tabulated luminosity density (lam_tab, dens_tab [ascending
  !--- in lambda]), restricted to [lam_lo, lam_hi] and to the range the table
  !--- covers.  Between its own nodes the table is read linearly in ln(lambda),
  !--- the model the rest of the code applies to spectra; knots(:) are extra
  !--- ascending node positions forced into the representation (the transfer bin
  !--- edges).  ok = .false. means the table and the requested band do not
  !--- overlap, or the enclosed luminosity vanishes.
  subroutine tabulated_spectrum_sampler(sampler, lam_tab, dens_tab, lam_lo, lam_hi, knots, ok)
  implicit none
  type(spectrum_sampler_type), intent(out) :: sampler
  real(kind=wp), intent(in)  :: lam_tab(:), dens_tab(:), lam_lo, lam_hi, knots(:)
  logical,       intent(out) :: ok
  real(kind=wp), allocatable :: node(:), dens(:)
  real(kind=wp) :: lo, hi
  integer :: i, n, nnode

  ok = .false.
  n  = size(lam_tab)
  if (n < 2 .or. size(dens_tab) /= n)   return
  if (any(dens_tab < 0.0_wp))           return
  do i = 1, n-1
     if (.not. (lam_tab(i+1) > lam_tab(i))) return
  enddo

  lo = max(lam_lo, lam_tab(1))
  hi = min(lam_hi, lam_tab(n))
  if (.not. (hi > lo)) return

  call merge_ascending_nodes(lam_tab, knots, lo, hi, node, nnode)
  allocate(dens(nnode))
  do i = 1, nnode
     dens(i) = log_lambda_interp(lam_tab, dens_tab, node(i))
  enddo

  call piecewise_linear_sampler(sampler, node(1:nnode), dens(1:nnode), ok)
  deallocate(node, dens)
  end subroutine tabulated_spectrum_sampler

  !---------------------------------------------------------------
  !--- Sampler for a blackbody of temperature T over [lam_lo, lam_hi].  The
  !--- interval between consecutive mandatory nodes (the band ends and the
  !--- transfer bin edges in knots) is bisected until the trapezoidal and
  !--- Simpson estimates of its luminosity agree to the relative tolerance rtol,
  !--- so the piecewise-linear density carries the blackbody to that accuracy.
  subroutine planck_spectrum_sampler(sampler, T, lam_lo, lam_hi, knots, ok, rtol)
  implicit none
  type(spectrum_sampler_type), intent(out) :: sampler
  real(kind=wp), intent(in)  :: T, lam_lo, lam_hi, knots(:)
  logical,       intent(out) :: ok
  real(kind=wp), intent(in), optional :: rtol
  real(kind=wp), allocatable :: mandatory(:), node(:), dens(:)
  real(kind=wp), allocatable :: empty(:)
  real(kind=wp) :: tol, scale, span
  integer :: i, nmand, nnode

  ok = .false.
  if (.not. (T > 0.0_wp) .or. .not. (lam_hi > lam_lo) .or. .not. (lam_lo > 0.0_wp)) return
  tol = 1.0e-8_wp
  if (present(rtol)) tol = rtol
  if (.not. (tol > 0.0_wp)) return

  !--- mandatory nodes: the band ends plus every knot strictly inside the band.
  allocate(empty(0))
  call merge_ascending_nodes(knots, empty, lam_lo, lam_hi, mandatory, nmand)
  deallocate(empty)
  if (nmand < 2) then
     deallocate(mandatory);  return
  endif

  scale = planck_band_estimate(T, lam_lo, lam_hi)
  span  = lam_hi - lam_lo
  if (.not. (scale > 0.0_wp)) then
     deallocate(mandatory);  return
  endif

  allocate(node(MAX_PLANCK_NODES), dens(MAX_PLANCK_NODES))
  nnode      = 1
  node(1)    = mandatory(1)
  dens(1)    = planck_lambda_shape(mandatory(1), T)
  do i = 1, nmand-1
     call refine_planck_interval(T, mandatory(i), mandatory(i+1), &
          planck_lambda_shape(mandatory(i),   T), &
          planck_lambda_shape(mandatory(i+1), T), &
          tol, scale, span, 0, node, dens, nnode, ok)
     if (.not. ok) then
        deallocate(mandatory, node, dens);  return
     endif
  enddo

  call piecewise_linear_sampler(sampler, node(1:nnode), dens(1:nnode), ok)
  deallocate(mandatory, node, dens)
  end subroutine planck_spectrum_sampler

  !---------------------------------------------------------------
  !--- Bisect [l0, l1] until the piecewise-linear luminosity of the interval is
  !--- converged, then append its right end.  The left end is already stored.
  recursive subroutine refine_planck_interval(T, l0, l1, f0, f1, tol, scale, span, &
                                              depth, node, dens, nnode, ok)
  implicit none
  real(kind=wp), intent(in)    :: T, l0, l1, f0, f1, tol, scale, span
  integer,       intent(in)    :: depth
  real(kind=wp), intent(inout) :: node(:), dens(:)
  integer,       intent(inout) :: nnode
  logical,       intent(inout) :: ok
  real(kind=wp) :: lm, fm, coarse, fine, err, allowance

  lm     = 0.5_wp*(l0 + l1)
  fm     = planck_lambda_shape(lm, T)
  coarse = 0.5_wp*(f0 + f1)*(l1 - l0)
  fine   = 0.25_wp*(f0 + 2.0_wp*fm + f1)*(l1 - l0)
  err    = abs(fine - coarse)/3.0_wp
  allowance = tol*(abs(fine) + scale*(l1 - l0)/span)

  if (err <= allowance .or. depth >= MAX_REFINE_DEPTH) then
     nnode = nnode + 1
     if (nnode > size(node)) then
        ok = .false.;  return
     endif
     node(nnode) = l1
     dens(nnode) = f1
     ok = .true.
  else
     call refine_planck_interval(T, l0, lm, f0, fm, tol, scale, span, depth+1, node, dens, nnode, ok)
     if (.not. ok) return
     call refine_planck_interval(T, lm, l1, fm, f1, tol, scale, span, depth+1, node, dens, nnode, ok)
  endif
  end subroutine refine_planck_interval

  !---------------------------------------------------------------
  !--- Coarse trapezoidal luminosity of a blackbody band, the scale against which
  !--- the refinement tolerance is measured.
  pure function planck_band_estimate(T, lam_lo, lam_hi) result(total)
  implicit none
  real(kind=wp), intent(in) :: T, lam_lo, lam_hi
  real(kind=wp) :: total, l0, l1
  integer, parameter :: n = 256
  integer :: i
  total = 0.0_wp
  do i = 1, n
     l0 = lam_lo + (lam_hi - lam_lo)*real(i-1, wp)/real(n, wp)
     l1 = lam_lo + (lam_hi - lam_lo)*real(i,   wp)/real(n, wp)
     total = total + 0.5_wp*(planck_lambda_shape(l0, T) + planck_lambda_shape(l1, T))*(l1 - l0)
  enddo
  end function planck_band_estimate

  !---------------------------------------------------------------
  !--- Store a piecewise-linear density and integrate it into a normalized
  !--- cumulative distribution.
  subroutine piecewise_linear_sampler(sampler, lam, dens, ok)
  implicit none
  type(spectrum_sampler_type), intent(out) :: sampler
  real(kind=wp), intent(in)  :: lam(:), dens(:)
  logical,       intent(out) :: ok
  integer :: i, n

  ok = .false.
  n  = size(lam)
  if (n < 2 .or. size(dens) /= n) return
  if (any(dens < 0.0_wp))         return
  do i = 1, n-1
     if (.not. (lam(i+1) > lam(i))) return
  enddo

  allocate(sampler%lam(n), sampler%dens(n), sampler%cdf(n))
  sampler%lam(:)  = lam(:)
  sampler%dens(:) = dens(:)
  sampler%cdf(1)  = 0.0_wp
  do i = 1, n-1
     sampler%cdf(i+1) = sampler%cdf(i) + 0.5_wp*(dens(i) + dens(i+1))*(lam(i+1) - lam(i))
  enddo
  sampler%total = sampler%cdf(n)
  if (.not. (sampler%total > 0.0_wp)) then
     deallocate(sampler%lam, sampler%dens, sampler%cdf)
     sampler%total = 0.0_wp
     return
  endif
  sampler%cdf(:) = sampler%cdf(:)/sampler%total
  sampler%cdf(n) = 1.0_wp
  ok = .true.
  end subroutine piecewise_linear_sampler

  !---------------------------------------------------------------
  !--- Draw a wavelength by inverting the cumulative distribution at u in (0,1).
  !--- Inside the bracketing interval the density is linear, so the cumulative
  !--- luminosity is quadratic in the offset and the inversion is closed-form.
  !--- The root is taken in the form 2A/(y0 + sqrt(y0^2 + 2 m A)), which stays
  !--- accurate when the slope m is small and the usual formula cancels.
  function sample_wavelength(sampler, u) result(lam)
  implicit none
  type(spectrum_sampler_type), intent(in) :: sampler
  real(kind=wp), intent(in) :: u
  real(kind=wp) :: lam, target, dlam, y0, y1, slope, area, part, disc, offset, ymax
  integer :: lo, hi, mid, n

  n = size(sampler%lam)
  target = min(max(u, 0.0_wp), 1.0_wp)
  if (target <= 0.0_wp) then
     lam = sampler%lam(1);  return
  else if (target >= 1.0_wp) then
     lam = sampler%lam(n);  return
  endif

  !--- bracket: largest lo with cdf(lo) <= target.
  lo = 1;  hi = n
  do while (hi - lo > 1)
     mid = (lo + hi)/2
     if (sampler%cdf(mid) <= target) then;  lo = mid;  else;  hi = mid;  endif
  enddo

  area = sampler%cdf(hi) - sampler%cdf(lo)
  if (.not. (area > 0.0_wp)) then
     lam = sampler%lam(hi);  return
  endif
  dlam = sampler%lam(hi) - sampler%lam(lo)
  y0   = sampler%dens(lo)
  y1   = sampler%dens(hi)
  !--- luminosity still to be covered inside this interval.
  part  = (target - sampler%cdf(lo))/area * (0.5_wp*(y0 + y1)*dlam)
  slope = (y1 - y0)/dlam
  ymax  = max(abs(y0), abs(y1), tiny(1.0_wp))
  if (abs(slope*dlam) <= 64.0_wp*epsilon(1.0_wp)*ymax) then
     if (y0 > 0.0_wp) then
        offset = part/y0
     else
        offset = 0.0_wp
     endif
  else
     disc   = sqrt(max(0.0_wp, y0*y0 + 2.0_wp*slope*part))
     offset = 2.0_wp*part/(y0 + disc)
  endif
  lam = sampler%lam(lo) + min(max(offset, 0.0_wp), dlam)
  end function sample_wavelength

  !---------------------------------------------------------------
  !--- Logarithmic mean of two positive values, (b-a)/ln(b/a), continued to
  !--- a = b.  A spectrum read as a power law across a bin has exactly this
  !--- mean of its two end values, weighted per logarithmic wavelength.
  pure function logarithmic_mean(a, b) result(m)
  implicit none
  real(kind=wp), intent(in) :: a, b
  real(kind=wp) :: m, r
  if (a <= 0.0_wp .or. b <= 0.0_wp) then
     m = 0.5_wp*(max(a, 0.0_wp) + max(b, 0.0_wp))
     return
  endif
  r = b/a
  if (abs(r - 1.0_wp) < 1.0e-6_wp) then
     !--- series in (r-1): the ratio (r-1)/ln(r) loses all digits at r = 1.
     m = a*(1.0_wp + 0.5_wp*(r - 1.0_wp) - (r - 1.0_wp)**2/12.0_wp)
  else
     m = a*(r - 1.0_wp)/log(r)
  endif
  end function logarithmic_mean

  !---------------------------------------------------------------
  !--- Draw a wavelength inside one bin whose spectrum is read as a power law:
  !--- lamj = lambda*j_lambda runs from lamj0 at lam0 to lamj1 at lam1, straight
  !--- in log-log.  This is the reading the tabulated spectra already get, and
  !--- across a bin spanning a decade of emissivity it is far closer than a
  !--- straight line in wavelength.  With x = ln(lambda/lam0)/ln(lam1/lam0) the
  !--- energy up to x is (q**x - 1)/(q - 1), q = lamj1/lamj0, which inverts in
  !--- closed form.
  pure function sample_power_law_bin(lam0, lam1, lamj0, lamj1, u) result(lam)
  implicit none
  real(kind=wp), intent(in) :: lam0, lam1, lamj0, lamj1, u
  real(kind=wp) :: lam, t, q, lnq, uu, dln

  uu  = min(max(u, 0.0_wp), 1.0_wp)
  dln = log(lam1/lam0)
  if (.not. (dln > 0.0_wp)) then
     lam = lam0;  return
  endif
  if (lamj0 <= 0.0_wp .or. lamj1 <= 0.0_wp) then
     !--- one end carries no emission: the power law is undefined there, so the
     !--- bin is sampled flat per logarithmic wavelength.
     lam = lam0*exp(uu*dln);  return
  endif

  q = lamj1/lamj0
  if (abs(q - 1.0_wp) < 1.0e-8_wp) then
     t = uu                                   ! uniform in ln(lambda)
  else
     lnq = log(q)
     t   = log(1.0_wp + uu*(q - 1.0_wp))/lnq
  endif
  lam = lam0*exp(min(max(t, 0.0_wp), 1.0_wp)*dln)
  end function sample_power_law_bin


  !---------------------------------------------------------------
  !--- Cumulative luminosity fraction below lambda.
  function spectrum_cdf(sampler, lam) result(f)
  implicit none
  type(spectrum_sampler_type), intent(in) :: sampler
  real(kind=wp), intent(in) :: lam
  real(kind=wp) :: f, dlam, off, y0, slope, whole, part
  integer :: lo, hi, mid, n

  n = size(sampler%lam)
  if (lam <= sampler%lam(1)) then
     f = 0.0_wp;  return
  else if (lam >= sampler%lam(n)) then
     f = 1.0_wp;  return
  endif
  lo = 1;  hi = n
  do while (hi - lo > 1)
     mid = (lo + hi)/2
     if (sampler%lam(mid) <= lam) then;  lo = mid;  else;  hi = mid;  endif
  enddo
  dlam  = sampler%lam(hi) - sampler%lam(lo)
  off   = lam - sampler%lam(lo)
  y0    = sampler%dens(lo)
  slope = (sampler%dens(hi) - y0)/dlam
  whole = 0.5_wp*(y0 + sampler%dens(hi))*dlam
  part  = y0*off + 0.5_wp*slope*off*off
  if (whole > 0.0_wp) then
     f = sampler%cdf(lo) + (sampler%cdf(hi) - sampler%cdf(lo))*part/whole
  else
     f = sampler%cdf(lo)
  endif
  end function spectrum_cdf

  !---------------------------------------------------------------
  !--- Luminosity fraction between two wavelengths.  Exact when both are nodes,
  !--- which is why the transfer bin edges are forced into every node set.
  function band_luminosity_fraction(sampler, lam1, lam2) result(f)
  implicit none
  type(spectrum_sampler_type), intent(in) :: sampler
  real(kind=wp), intent(in) :: lam1, lam2
  real(kind=wp) :: f
  f = spectrum_cdf(sampler, lam2) - spectrum_cdf(sampler, lam1)
  if (f < 0.0_wp) f = 0.0_wp
  end function band_luminosity_fraction

  !---------------------------------------------------------------
  !--- Merge two ascending node lists, keep what lies in [lo, hi], and add both
  !--- ends.  Values closer than a few rounding units are treated as one node,
  !--- since the cumulative distribution needs strictly increasing wavelengths.
  subroutine merge_ascending_nodes(a, b, lo, hi, out, nout)
  implicit none
  real(kind=wp), intent(in) :: a(:), b(:), lo, hi
  real(kind=wp), allocatable, intent(out) :: out(:)
  integer, intent(out) :: nout
  real(kind=wp) :: next
  integer :: ia, ib, na, nb

  na = size(a);  nb = size(b)
  allocate(out(na + nb + 2))
  nout = 1;  out(1) = lo
  ia = 1;  ib = 1
  do
     do while (ia <= na)
        if (a(ia) > lo) exit
        ia = ia + 1
     enddo
     do while (ib <= nb)
        if (b(ib) > lo) exit
        ib = ib + 1
     enddo
     if (ia > na .and. ib > nb) exit
     if (ib > nb) then
        next = a(ia);  ia = ia + 1
     else if (ia > na) then
        next = b(ib);  ib = ib + 1
     else if (a(ia) <= b(ib)) then
        next = a(ia);  ia = ia + 1
     else
        next = b(ib);  ib = ib + 1
     endif
     if (next >= hi) exit
     if (next > out(nout)*(1.0_wp + 8.0_wp*epsilon(1.0_wp))) then
        nout = nout + 1
        out(nout) = next
     endif
  enddo
  if (hi > out(nout)*(1.0_wp + 8.0_wp*epsilon(1.0_wp))) then
     nout = nout + 1
     out(nout) = hi
  else
     out(nout) = hi
  endif
  end subroutine merge_ascending_nodes

  !---------------------------------------------------------------
  !--- Read a tabulated spectrum linearly in ln(lambda), clamped at the ends.
  pure function log_lambda_interp(x, y, xnew) result(ynew)
  implicit none
  real(kind=wp), intent(in) :: x(:), y(:), xnew
  real(kind=wp) :: ynew, t, lx
  integer :: n, lo, hi, mid
  n = size(x)
  if (xnew <= x(1)) then
     ynew = y(1);  return
  else if (xnew >= x(n)) then
     ynew = y(n);  return
  endif
  lo = 1;  hi = n
  do while (hi - lo > 1)
     mid = (lo + hi)/2
     if (x(mid) <= xnew) then;  lo = mid;  else;  hi = mid;  endif
  enddo
  lx = log(xnew)
  t  = (lx - log(x(lo)))/(log(x(hi)) - log(x(lo)))
  ynew = y(lo) + t*(y(hi) - y(lo))
  if (ynew < 0.0_wp) ynew = 0.0_wp
  end function log_lambda_interp

end module spectrum_sampler_mod
