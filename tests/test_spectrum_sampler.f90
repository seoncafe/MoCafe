program test_spectrum_sampler
!--- Checks of the continuous wavelength sampler:
!---   1. a piecewise-linear spectrum is reproduced exactly by the inversion
!---   2. the inversion and the forward cumulative distribution are inverses
!---   3. a blackbody is carried to the requested tolerance
!---   4. the closed-form bin index brackets the sampled wavelength
!---   5. the size of the bias the bin-center scheme carried
  use define, only : wp
  use spectrum_sampler_mod
  implicit none

  type(spectrum_sampler_type) :: s
  real(kind=wp), allocatable  :: edge(:), lam_t(:), dens_t(:), frac(:), cnt(:)
  real(kind=wp) :: lmin, lmax, dln, u, lam, f, worst, err, T
  real(kind=wp) :: exact, mid, rel, relmax, relsum, meanl, meanl_ref
  integer :: nlam, i, il, n, nbad, ntest
  logical :: ok

  write(*,'(a)') '================ continuous wavelength sampler ================'

  !--------------------------------------------------------------
  !--- 1. piecewise-linear spectrum: the inversion must be exact.
  !--------------------------------------------------------------
  lmin = 0.1_wp;  lmax = 3.0_wp;  nlam = 20
  allocate(edge(nlam+1), frac(nlam), cnt(nlam))
  dln = log(lmax/lmin)/nlam
  do i = 1, nlam+1
     edge(i) = lmin*exp((i-1)*dln)
  enddo
  edge(nlam+1) = lmax

  allocate(lam_t(4), dens_t(4))
  lam_t  = [0.10_wp, 0.30_wp, 1.00_wp, 3.00_wp]
  dens_t = [1.00_wp, 5.00_wp, 2.00_wp, 0.50_wp]
  call tabulated_spectrum_sampler(s, lam_t, dens_t, lmin, lmax, edge, ok)
  write(*,'(a,l2,a,i7)') ' [1] tabulated sampler built: ', ok, '   nodes =', size(s%lam)

  do i = 1, nlam
     frac(i) = band_luminosity_fraction(s, edge(i), edge(i+1))
  enddo
  write(*,'(a,es12.5)') '     sum of bin fractions - 1        : ', sum(frac) - 1.0_wp

  !--- stratified quantiles: bin counts must match N*fraction to within one draw.
  n = 10000000
  cnt(:) = 0.0_wp
  do i = 1, n
     u   = (real(i,wp) - 0.5_wp)/real(n,wp)
     lam = sample_wavelength(s, u)
     il  = floor(log(lam/lmin)/dln) + 1
     if (il < 1)    il = 1
     if (il > nlam) il = nlam
     cnt(il) = cnt(il) + 1.0_wp
  enddo
  worst = 0.0_wp
  do i = 1, nlam
     worst = max(worst, abs(cnt(i) - frac(i)*real(n,wp)))
  enddo
  write(*,'(a,i10)')    '     stratified draws                : ', n
  write(*,'(a,f12.4)')  '     worst |count - N*fraction|      : ', worst
  write(*,'(a)')        '     (quantile sampling: must be <= 1 draw)'

  !--------------------------------------------------------------
  !--- 2. inversion and forward cumulative distribution are inverses.
  !--------------------------------------------------------------
  worst = 0.0_wp;  ntest = 200001
  do i = 1, ntest
     u   = real(i,wp)/real(ntest+1,wp)
     lam = sample_wavelength(s, u)
     f   = spectrum_cdf(s, lam)
     worst = max(worst, abs(f - u))
  enddo
  write(*,'(a,es12.5)') ' [2] max |cdf(sample(u)) - u|        : ', worst

  !--------------------------------------------------------------
  !--- 3. blackbody: band fractions against a fine numerical integral.
  !--------------------------------------------------------------
  deallocate(edge, frac, cnt, lam_t, dens_t)
  T = 30000.0_wp
  lmin = 0.0912_wp;  lmax = 10.0_wp;  nlam = 100
  allocate(edge(nlam+1), frac(nlam))
  dln = log(lmax/lmin)/nlam
  do i = 1, nlam+1
     edge(i) = lmin*exp((i-1)*dln)
  enddo
  edge(nlam+1) = lmax
  call planck_spectrum_sampler(s, T, lmin, lmax, edge, ok, rtol=1.0e-10_wp)
  write(*,'(a,l2,a,i7)') ' [3] Planck sampler built    : ', ok, '   nodes =', size(s%lam)

  relmax = 0.0_wp;  relsum = 0.0_wp
  do i = 1, nlam
     frac(i) = band_luminosity_fraction(s, edge(i), edge(i+1))
     exact   = planck_band_integral(T, edge(i), edge(i+1))/planck_band_integral(T, lmin, lmax)
     rel     = abs(frac(i) - exact)/max(exact, 1.0e-300_wp)
     if (exact > 1.0e-8_wp) then
        relmax = max(relmax, rel)
        relsum = relsum + rel
     endif
  enddo
  write(*,'(a,es12.5)') '     max rel. error of bin fraction  : ', relmax
  write(*,'(a,es12.5)') '     mean rel. error of bin fraction : ', relsum/real(nlam,wp)

  !--- mean wavelength from stratified draws against the analytic mean.
  n = 2000000;  meanl = 0.0_wp
  do i = 1, n
     u = (real(i,wp) - 0.5_wp)/real(n,wp)
     meanl = meanl + sample_wavelength(s, u)
  enddo
  meanl = meanl/real(n,wp)
  meanl_ref = planck_band_moment(T, lmin, lmax)/planck_band_integral(T, lmin, lmax)
  write(*,'(a,2es14.6)') '     mean lambda: sampled, exact     : ', meanl, meanl_ref
  write(*,'(a,es12.5)')  '     rel. error                      : ', abs(meanl-meanl_ref)/meanl_ref

  !--------------------------------------------------------------
  !--- 4. the closed-form bin index must bracket the wavelength.
  !--------------------------------------------------------------
  nbad = 0;  n = 2000000
  do i = 1, n
     u   = (real(i,wp) - 0.5_wp)/real(n,wp)
     lam = sample_wavelength(s, u)
     il  = floor(log(lam/lmin)/dln) + 1
     if (il < 1)    il = 1
     if (il > nlam) il = nlam
     if (lam < edge(il) .or. lam > edge(il+1)) nbad = nbad + 1
  enddo
  write(*,'(a,i10)') ' [4] samples outside their own bin   : ', nbad
  il = bin_index(lmin, lmin, dln, nlam);  write(*,'(a,i6)') '     bin of lambda_min               : ', il
  il = bin_index(lmax, lmin, dln, nlam);  write(*,'(a,i6)') '     bin of lambda_max (must be nlam): ', il

  !--------------------------------------------------------------
  !--- 5. the bias of the bin-center scheme, on this same blackbody.
  !--------------------------------------------------------------
  relmax = 0.0_wp;  relsum = 0.0_wp;  ntest = 0
  exact = planck_band_integral(T, lmin, lmax)
  do i = 1, nlam
     mid = planck_lambda_shape(sqrt(edge(i)*edge(i+1)), T)*(edge(i+1) - edge(i))
     err = planck_band_integral(T, edge(i), edge(i+1))
     if (err/exact > 1.0e-6_wp) then
        rel    = abs(mid - err)/err
        relmax = max(relmax, rel)
        relsum = relsum + rel
        ntest  = ntest + 1
     endif
  enddo
  write(*,'(a)')        ' [5] bin-center luminosity vs exact band integral'
  write(*,'(a,i5,a,i5)')'     bins carrying weight            : ', ntest, ' of ', nlam
  write(*,'(a,f10.4,a)')'     worst relative error            : ', 100.0_wp*relmax, ' %'
  write(*,'(a,f10.4,a)')'     mean  relative error            : ', 100.0_wp*relsum/real(ntest,wp), ' %'

  deallocate(edge, frac)
  write(*,'(a)') '=============================================================='

contains

  !--- Planck band luminosity by fine Simpson integration, the reference the
  !--- sampler is measured against.
  function planck_band_integral(T, l0, l1) result(v)
  implicit none
  real(kind=wp), intent(in) :: T, l0, l1
  real(kind=wp) :: v, h, x
  integer, parameter :: m = 20000
  integer :: k
  h = (l1 - l0)/real(m,wp)
  v = planck_lambda_shape(l0, T) + planck_lambda_shape(l1, T)
  do k = 1, m-1
     x = l0 + real(k,wp)*h
     if (mod(k,2) == 1) then
        v = v + 4.0_wp*planck_lambda_shape(x, T)
     else
        v = v + 2.0_wp*planck_lambda_shape(x, T)
     endif
  enddo
  v = v*h/3.0_wp
  end function planck_band_integral

  function planck_band_moment(T, l0, l1) result(v)
  implicit none
  real(kind=wp), intent(in) :: T, l0, l1
  real(kind=wp) :: v, h, x
  integer, parameter :: m = 20000
  integer :: k
  h = (l1 - l0)/real(m,wp)
  v = l0*planck_lambda_shape(l0, T) + l1*planck_lambda_shape(l1, T)
  do k = 1, m-1
     x = l0 + real(k,wp)*h
     if (mod(k,2) == 1) then
        v = v + 4.0_wp*x*planck_lambda_shape(x, T)
     else
        v = v + 2.0_wp*x*planck_lambda_shape(x, T)
     endif
  enddo
  v = v*h/3.0_wp
  end function planck_band_moment

  pure function bin_index(lam, lmin_in, dln_in, nlam_in) result(il)
  implicit none
  real(kind=wp), intent(in) :: lam, lmin_in, dln_in
  integer,       intent(in) :: nlam_in
  integer :: il
  il = floor(log(lam/lmin_in)/dln_in) + 1
  if (il < 1)       il = 1
  if (il > nlam_in) il = nlam_in
  end function bin_index

end program test_spectrum_sampler
