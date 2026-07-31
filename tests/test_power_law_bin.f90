program test_power_law_bin
!--- Checks of the within-bin power-law draw used by the dust emission:
!---   1. the draw never leaves its bin, and the ends map to the ends
!---   2. the sampled distribution matches the analytic power-law CDF
!---   3. the logarithmic mean reproduces the analytic bin integral
!---   4. the degenerate cases (equal ends, a zero end) behave
  use define, only : wp
  use spectrum_sampler_mod
  implicit none

  real(kind=wp) :: l0, l1, a, b, u, lam, t, q, want, got, worst, dln
  real(kind=wp) :: exact, viaLog
  integer :: i, n, nbad

  write(*,'(a)') '============ within-bin power-law draw ============'

  l0 = 1.0_wp;  l1 = 1.1051709180756477_wp     ! dln = 0.1
  dln = log(l1/l0)

  !--- 1 & 2: a steep bin (the end values differ by 30x).
  a = 1.0_wp;  b = 30.0_wp
  q = b/a
  n = 2000001;  worst = 0.0_wp;  nbad = 0
  do i = 1, n
     u   = real(i,wp)/real(n+1,wp)
     lam = sample_power_law_bin(l0, l1, a, b, u)
     if (lam < l0 .or. lam > l1) nbad = nbad + 1
     !--- analytic CDF: with x = ln(lam/l0)/dln, F = (q**x - 1)/(q - 1)
     t    = log(lam/l0)/dln
     want = u
     got  = (q**t - 1.0_wp)/(q - 1.0_wp)
     worst = max(worst, abs(got - want))
  enddo
  write(*,'(a,i10)')    ' [1] draws outside their bin     : ', nbad
  write(*,'(a,2es14.6)')'     lam(u=0), lam(u=1)          : ', &
       sample_power_law_bin(l0,l1,a,b,0.0_wp), sample_power_law_bin(l0,l1,a,b,1.0_wp)
  write(*,'(a,2es14.6)')'     bin ends                    : ', l0, l1
  write(*,'(a,es14.6)') ' [2] max |CDF(sample(u)) - u|    : ', worst

  !--- 3: the bin integral.  For lam*j = a*(b/a)**x the energy is
  !---    Int j dlam = Int (lam*j) dln(lam) = dln * (b-a)/ln(b/a).
  exact  = dln*(b - a)/log(b/a)
  viaLog = dln*logarithmic_mean(a, b)
  write(*,'(a,2es16.8)') ' [3] bin integral: analytic, code: ', exact, viaLog
  write(*,'(a,es14.6)')  '     relative difference         : ', abs(viaLog-exact)/exact

  !--- 4: degenerate ends.
  write(*,'(a)') ' [4] degenerate ends'
  write(*,'(a,2es14.6)') '     a=b=2: logmean vs 2        : ', logarithmic_mean(2.0_wp,2.0_wp), 2.0_wp
  write(*,'(a,es14.6)')  '     a=b: draw at u=0.5 (mid ln): ', sample_power_law_bin(l0,l1,2.0_wp,2.0_wp,0.5_wp)
  write(*,'(a,es14.6)')  '     expected l0*exp(0.5*dln)   : ', l0*exp(0.5_wp*dln)
  write(*,'(a,es14.6)')  '     a=0: draw at u=0.5         : ', sample_power_law_bin(l0,l1,0.0_wp,3.0_wp,0.5_wp)
  write(*,'(a,es14.6)')  '     ratio b/a near 1 (1+1e-9)  : ', logarithmic_mean(1.0_wp,1.0_wp+1.0e-9_wp)

  !--- a shallow bin, where the series branch of the logarithmic mean is used.
  a = 1.0_wp;  b = 1.0_wp + 1.0e-7_wp
  exact = dln*(b - a)/log(b/a)
  viaLog = dln*logarithmic_mean(a, b)
  write(*,'(a,es14.6)') '     shallow bin rel. difference: ', abs(viaLog-exact)/exact

  write(*,'(a)') '=================================================='
end program test_power_law_bin
