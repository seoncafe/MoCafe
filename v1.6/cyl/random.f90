module random
!-----
! Random Number Generator Module
! Marsaglis KISS (Keep It Simple and Stupid) 64 bits RNGs
!
! Author: Kwang-Il Seon
!
! Note:
!    - gfortran random_number uses KISS algorithm while ifort uses Lecuyer algorithm.
!    - intrinsic random_number routine is not thread-safe with openmp.
!
! Usage:
!    call init_random_seed()
!    call random_number(harvest) or harvest = rand_number()
!         where harvest is a number (32- or 64-bit) or an array.
!
! History:
!   v2b(2016/02/20)
!        added random_mvnorm and make_mvcovmat2
!   v2a(2016/02/16)
!        added rand_t, rand_gamma, rand_scaled_inv_chi2
!   v2 (2016/01/15)
!        modified to use the same calling method as the intrinsic procedure defined in fortran.
!   2015/12/11,
!        To make "save" variables thread safe, they should be declared in THREADPRIVATE directive.
!        http://sc.tamu.edu/IBM.Tutorial/docs/Compilers/xlf_8.1/html/lr225.HTM
!   2012/05/04, Initial version
!-----
  use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64
  implicit none
  private
  integer, parameter  :: sp = real32
  integer, parameter  :: dp = real64
  integer, parameter  :: wp = dp
  real(wp), parameter :: PI = 3.141592653589793238462643383279502884197_wp

  public init_random_seed
  public random_seed, random_number, random_gauss
  public rand_number, rand_gauss, rand_exp, rand_rexp, rand_zexp, rand_sech2
  public rand_permutation, rand_cyclic_permutation, rand_index, rand_binomial, rand_multinomial
  public rand_t, rand_gamma, rand_scaled_inv_chi2
  public make_mvcovmat2, random_mvnorm, make_mvcovmat

  integer(int64), save :: seed_kiss(4)   = [ 1234567890987654321_int64, 362436362436362436_int64,&
                                                1066149217761810_int64, 123456123456123456_int64 ]
  integer(int32), save :: seed_kiss32(4) = [ 123456789, 362436069, 521288629, 916191069 ]
  !$OMP THREADPRIVATE(seed_kiss,seed_kiss32)

  interface random_seed
     module procedure random_kiss_seed
  end interface random_seed

  interface random_number
     module procedure random_kiss32,random_kiss64,random_kiss_v32,random_kiss_v64,&
                      random_kiss_vv32,random_kiss_vv64,random_kiss_vvv32,random_kiss_vvv64
  end interface random_number

  interface random_gauss
     module procedure random_gauss1,random_gauss2,random_gauss3
  end interface random_gauss

  interface rand_number
     module procedure kiss64
  end interface rand_number

  interface rand_gauss
     module procedure rand_gauss1
  end interface rand_gauss

  interface rand_gamma
     module procedure gamdev, gamdev_r
  end interface rand_gamma

  interface rand_exp
     module procedure rand_exp32,rand_exp64
  end interface rand_exp

  interface rand_index
     module procedure rand_index1, rand_indexn
  end interface rand_index
contains
!==================================================================
  subroutine random_kiss32(harvest)
     implicit none
     real(sp), intent(out) :: harvest
     harvest = kiss32()
  end subroutine random_kiss32
  subroutine random_kiss64(harvest)
     implicit none
     real(dp), intent(out) :: harvest
     harvest = kiss64()
  end subroutine random_kiss64
  subroutine random_kiss_v32(harvest)
     implicit none
     real(sp), intent(out) :: harvest(:)
     integer :: i
     do i=1, size(harvest)
        harvest(i) = kiss32()
     enddo
  end subroutine random_kiss_v32
  subroutine random_kiss_v64(harvest)
     implicit none
     real(dp), intent(out) :: harvest(:)
     integer :: i
     do i=1, size(harvest)
        harvest(i) = kiss64()
     enddo
  end subroutine random_kiss_v64
  subroutine random_kiss_vv32(harvest)
     implicit none
     real(sp), intent(out) :: harvest(:,:)
     real(sp) :: vv(size(harvest))
     call random_kiss_v32(vv)
     harvest = reshape(vv, [size(harvest,1), size(harvest,2)])
  end subroutine random_kiss_vv32
  subroutine random_kiss_vv64(harvest)
     implicit none
     real(dp), intent(out) :: harvest(:,:)
     real(dp) :: vv(size(harvest))
     call random_kiss_v64(vv)
     harvest = reshape(vv, [size(harvest,1), size(harvest,2)])
  end subroutine random_kiss_vv64
  subroutine random_kiss_vvv32(harvest)
     implicit none
     real(sp), intent(out) :: harvest(:,:,:)
     real(sp) :: vv(size(harvest))
     call random_kiss_v32(vv)
     harvest = reshape(vv, [size(harvest,1), size(harvest,2), size(harvest,3)])
  end subroutine random_kiss_vvv32
  subroutine random_kiss_vvv64(harvest)
     implicit none
     real(dp), intent(out) :: harvest(:,:,:)
     real(dp) :: vv(size(harvest))
     call random_kiss_v64(vv)
     harvest = reshape(vv, [size(harvest,1), size(harvest,2), size(harvest,3)])
  end subroutine random_kiss_vvv64
  !----------------------------------
  function kiss32() result(u)
  !------------------------------------------------------------------
  ! Pseudo-random uniform 32-bit integer random number generator
  ! The Marsaglia KISS (Keep It Simple Stupid) random number generator
  ! returns a sequence of pseudo-random 32-bit integers. It combines:
  ! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
  ! (2) A 3-shift shift-register generator, period 2^32-1,
  ! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497 > 2^59
  !  The overall period is > 2^123. Consider this a mininum standard for
  ! generating random integers.  Integers range from -huge(0) to +huge(0)
  ! (this is the range [-2^31-1, 2^31-1] for 32-bit integers.
  !
  ! http://computer-programming-forum.com/49-fortran/b89977aa62f72ee8.htm
  !-----------------
  ! Kwang-Il Seon, 2013 April 06
  !    Modified to return 0.0 < u < 1.0 (u1 <= u <= u2)
  !------------------------------------------------------------------
    implicit none
    integer :: m,k,n
    integer :: kiss
    real(sp) :: u

    real(sp), parameter :: x1 = 129.0, x2 = 257.0
    real(sp), parameter :: a  = 1.0/(2.0*real(huge(0),sp)+1.0+x1+x2)
    real(sp), parameter :: b  = (real(huge(0),sp)+1.0+x1)/(2.0*real(huge(0),sp)+1.0+x1+x2)

    m(k, n) = ieor(k, ishft(k, n))
    seed_kiss32(1) = 69069 * seed_kiss32(1) + 1327217885
    seed_kiss32(2) = m(m(m(seed_kiss32(2), 13), -17), 5)
    seed_kiss32(3) = 18000 * iand(seed_kiss32(3), 65535) + ishft(seed_kiss32(3), - 16)
    seed_kiss32(4) = 30903 * iand(seed_kiss32(4), 65535) + ishft(seed_kiss32(4), - 16)
    kiss = seed_kiss32(1) + seed_kiss32(2) + ishft(seed_kiss32(3), 16) + seed_kiss32(4)
    u    = b + a*real(kiss,sp)
  end function kiss32
  !----------------------------------
  function kiss64() result(u)
  !--------------------------------------------------------------------
  ! The 64-bit KISS (Keep It Simple Stupid) random number generator.
  ! Components:
  !  (1) Xorshift (XSH), period 2^64-1,
  !  (2) Multiply-with-carry (MWC), period (2^121+2^63-1)
  !  (3) Congruential generator (CNG), period 2^64.
  ! Overall period: (2^250+2^192+2^64-2^186-2^129)/6 ~= 2^(247.42) or 10^(74.48)
  !
  ! https://groups.google.com/forum/#!msg/comp.lang.fortran/qFv18ql_WlU/IK8KGZZFJx4J
  ! or http://fortranwiki.org/fortran/show/kiss64
  !-----------------
  ! Kwang-Il Seon, 2014 October 08, 2015-06-03
  !    Modified to return 0.0 < u < 1.0.
  !--------------------------------------------------------------------
    implicit none
    integer(int64) :: kiss,t
    integer(int64) :: m,s,x,k
    real(int64)    :: u

    real(dp), parameter :: x1 = 1025.0_dp, x2 = 2049.0_dp
    real(dp), parameter :: a  = 1.0_dp/(2.0_dp*real(huge(0_int64),dp)+1.0_dp+x1+x2)
    real(dp), parameter :: b  = (real(huge(0_int64),dp)+1.0_dp+x1)/(2.0_dp*real(huge(0_int64),dp)+1.0_dp+x1+x2)
    ! statement functions
    m(x,k) = ieor(x,ishft(x,k))
    s(x)   = ishft(x,-63)

    t = ishft(seed_kiss(1),58) + seed_kiss(4)
    if (s(seed_kiss(1)) == s(t)) then
      seed_kiss(4) = ishft(seed_kiss(1),-6) + s(seed_kiss(1))
    else
      seed_kiss(4) = ishft(seed_kiss(1),-6) + 1 - s(seed_kiss(1)+t)
    endif
    seed_kiss(1) = t + seed_kiss(1)
    seed_kiss(2) = m(m(m(seed_kiss(2),13_int64), -17_int64), 43_int64)
    seed_kiss(3) = 6906969069_int64 * seed_kiss(3) + 1234567_int64
    kiss         = seed_kiss(1) + seed_kiss(2) + seed_kiss(3)
    u            = b + a*real(kiss,dp)
  end function kiss64
!----------------------------------
  subroutine random_kiss_seed_(seed_size,put,get)
    ! size is an intrinsic function so that we need to define an intermediate routine
    ! which use the argument name "seed_size" instead of "size"
    ! random_kiss_seed has the optional argument "size"
    implicit none
    integer, optional, intent(out) :: seed_size
    integer, optional, intent(in)  :: put(:)
    integer, optional, intent(out) :: get(:)
    ! local variable
    integer :: n

    if (present(seed_size)) then
       seed_size = 4
    endif
    if (present(put)) then
       n              = minval([4,size(put)])
       seed_kiss(1:n) = put(1:n)
    endif
    if (present(get)) then
       n        = minval([4,size(get)])
       get(1:n) = seed_kiss(1:n)
    endif
  end subroutine random_kiss_seed_
  subroutine random_kiss_seed(size,put,get)
    implicit none
    integer, optional, intent(out) :: size
    integer, optional, intent(in)  :: put(:)
    integer, optional, intent(out) :: get(:)
    call random_kiss_seed_(seed_size=size,put=put,get=get)
  end subroutine random_kiss_seed
!----------
  subroutine init_random_seed()
#ifdef _OPENMP
    use omp_lib
#endif
    implicit none
    integer, allocatable :: seed(:)
    integer        :: i, n, un, istat, dt(8), pid
    integer(int64) :: t
    integer        :: getpid

    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
    ! return an integer uniform random deviate between [0 and huge(0)] from /dev/urandom if successful.
    open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(t)
       if (t == 0) then
          call date_and_time(values=dt)
          t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
       end if
#ifdef _OPENMP
       ! OPENMP does not launches multiple instances.
       pid = getpid() + omp_get_thread_num()
#else
       pid = getpid()
#endif
       t   = ieor(t, int(pid, kind(t)))
       do i = 1, n
          seed(i) = shr3(t)
       end do
    end if
    call random_seed(put=seed)
    contains
      ! This is sufficient for seeding a better RNG.
      function shr3(iseed) result(iran)
        implicit none
        integer(int64), intent(inout) :: iseed
        integer(int64) :: iran
        integer(int64) :: iseed0
        integer(int64) :: m,k,n
        m(k, n) = ieor(k, ishft(k, n))
        iseed0  = iseed
        iseed   = m(m(m(iseed, 13_int64), -17_int64), 43_int64)
        iran    = iseed0 + iseed
      end function shr3
  end subroutine init_random_seed
  !-------------------------------------
  !   Generate a random normal deviate using the polar method.
  !   Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for generating
  !              normal variables', Siam Rev., vol.6, 260-264, 1964.
  !           or Numerical Recipes
  !-------------------------------------
  function rand_gauss1() result(gasdev)
    implicit none
    real(dp) :: gasdev
    real(dp) :: rsq,v1,v2
    real(dp), save :: gset
    logical,  save :: gaus_stored = .false.
    real(dp), parameter :: one = 1.0_dp
    !$OMP THREADPRIVATE(gset,gaus_stored)

    if (gaus_stored) then
       gasdev = gset
       gaus_stored = .false.
    else
       do
          v1  = 2.0_dp*rand_number() - 1.0_dp
          v2  = 2.0_dp*rand_number() - 1.0_dp
          rsq = v1*v1 + v2*v2
          if (rsq > 0.0 .and. rsq < one) exit
       enddo
       rsq    = sqrt(-2.0_dp*log(rsq)/rsq)
       gset   = v1*rsq
       gasdev = v2*rsq
       gaus_stored = .true.
    endif
  end function rand_gauss1
  subroutine random_gauss1(harvest)
    implicit none
    real(wp), intent(out) :: harvest
    real(wp) :: rsq,v1,v2
    real(wp), save :: g
    logical,  save :: gaus_stored=.false.
    !$OMP THREADPRIVATE(g,gaus_stored)
    if (gaus_stored) then
       harvest     = g
       gaus_stored = .false.
    else
       do
          v1  = 2.0_sp*rand_number()-1.0_sp
          v2  = 2.0_sp*rand_number()-1.0_sp
          rsq = v1**2+v2**2
          if (rsq > 0.0 .and. rsq < 1.0) exit
       end do
       rsq         = sqrt(-2.0_sp*log(rsq)/rsq)
       harvest     = v1*rsq
       g           = v2*rsq
       gaus_stored = .true.
    end if
  end subroutine random_gauss1
  subroutine random_gauss2(harvest)
    implicit none
    real(wp), dimension(:), intent(out) :: harvest
    real(wp), dimension(size(harvest))  :: rsq,v1,v2
    real(wp), allocatable, dimension(:), save :: g
    integer :: n,ng,nn,m
    integer, save :: last_allocated=0
    logical, save :: gaus_stored=.false.
    logical, dimension(size(harvest)) :: mask
    !$OMP THREADPRIVATE(g,last_allocated,gaus_stored)

    n = size(harvest)
    if (n /= last_allocated) then
       if (last_allocated /= 0) deallocate(g)
       allocate(g(n))
       last_allocated = n
       gaus_stored    = .false.
    end if
    if (gaus_stored) then
       harvest     = g
       gaus_stored = .false.
    else
       ng=1
       do
          if (ng > n) exit
          call random_number(v1(ng:n))
          call random_number(v2(ng:n))
          v1(ng:n)   = 2.0_sp*v1(ng:n)-1.0_sp
          v2(ng:n)   = 2.0_sp*v2(ng:n)-1.0_sp
          rsq(ng:n)  = v1(ng:n)**2+v2(ng:n)**2
          mask(ng:n) = (rsq(ng:n)>0.0 .and. rsq(ng:n)<1.0)
          call array_copy(pack(v1(ng:n),mask(ng:n)),v1(ng:),nn,m)
          v2(ng:ng+nn-1)  = pack(v2(ng:n),mask(ng:n))
          rsq(ng:ng+nn-1) = pack(rsq(ng:n),mask(ng:n))
          ng              = ng+nn
       end do
       rsq     = sqrt(-2.0_sp*log(rsq)/rsq)
       harvest = v1*rsq
       g       = v2*rsq
       gaus_stored = .true.
    end if
    contains
      subroutine array_copy(src,dest,n_copied,n_not_copied)
        real(WP), intent(in)  :: src(:)
        real(WP), intent(out) :: dest(:)
        integer,  intent(out) :: n_copied, n_not_copied

        n_copied         = min(size(src),size(dest))
        n_not_copied     = size(src)-n_copied
        dest(1:n_copied) = src(1:n_copied)
      end subroutine array_copy
  end subroutine random_gauss2
  subroutine random_gauss3(harvest)
    implicit none
    real(wp), intent(out) :: harvest(:,:)
    real(wp)              :: rand_g(size(harvest))
    integer :: n1,n2

    call random_gauss2(rand_g)
    n1 = size(harvest,1)
    n2 = size(harvest,2)
    harvest = reshape(rand_g,(/n1,n2/))
  end subroutine random_gauss3
  !------------------------------
  ! Trucated exponential distribution
  !------------------------------
  function rand_exp32(x0,xcut) result(expdev)
    implicit none
    real, intent(in) :: x0,xcut

    real :: expdev
    !local variables
    !real :: ra,rtemp
    !real, parameter :: eps = epsilon(0.0)

    !--- Truncated exponential
    ! the following algorithm gives inf number sometimes.
    !ra = rand_number()-0.5
    !rtemp = 2.0*ra*(1.0-exp(-xcut/x0))
    !if (abs(1.0-rtemp) > eps) then
    !if (rtemp /= 1.0) then
    !   expdev = -x0*sign(1.0,ra)*log(1.0-sign(1.0,ra)*rtemp)
    !else
    expdev = -x0 * log(rand_number())
    do while (expdev > xcut)
       expdev = -x0 * log(rand_number())
    enddo
    expdev = sign(1.0,rand_number()-0.5)*expdev
    !endif

    return
  end function rand_exp32
  function rand_exp64(x0,xcut) result(expdev)
    implicit none
    real(dp), intent(in) :: x0,xcut
    real(dp) :: expdev

    !local variables
    !real(dp) :: ra
    !--- Truncated exponential
    !ra = rand_number()-0.5_dp
    !expdev = -x0*sign(1.0_dp,ra)*log(1.0_dp-sign(1.0_dp,ra)*2.0_dp*ra*(1.0_dp-exp(-xcut/x0)))
    expdev = -x0 * log(rand_number())
    do while (expdev > xcut)
       expdev = -x0 * log(rand_number())
    enddo
    expdev = sign(1.0_dp,rand_number()-0.5_dp)*expdev

    return
  end function rand_exp64
  !--------------------
  ! random number generator for sech^2 distribution in z-direction.
  ! pdf = sech^2(z) when abs(z) < zmax
  !     = 0         otherwise
  function rand_sech2(zmax) result(r)
    implicit none
    real(dp), intent(in) :: zmax
    real(dp), save :: tanhzmax, zmax_save=-999.0_dp
    real(dp) :: r
    !$OMP THREADPRIVATE(tanhzmax,zmax_save)
    if (zmax_save /= zmax) then
       tanhzmax = tanh(zmax)
       zmax_save = zmax
    endif
    r = atanh(tanhzmax*(2.0_dp*rand_number()-1.0_dp))
    return
  end function rand_sech2
  !--------------------
  ! random number generator for exponential distribution in z-direction.
  ! pdf = exp(-abs(z)) when abs(z) < zmax
  !     = 0            otherwise
  function rand_zexp(zmax) result(r)
    implicit none
    real(dp), intent(in) :: zmax
    real(dp), save :: zmax_save=-999.0_dp, q_save
    real(dp) :: r
    !$OMP THREADPRIVATE(zmax_save,q_save)
    if (zmax_save /= zmax) then
       q_save = 1.0_dp - exp(-zmax)
       zmax_save = zmax
    endif
    r = rand_number()-0.5_dp
    r = -sign(1.0_dp,r)*log(1.0_dp-sign(1.0_dp,r)*2.0_dp*r*q_save)
    return
  end function rand_zexp
  !--------------------
  ! random number generator for exponential distribution in radial direction.
  ! p(r) dr = r x exp(-r) dr when r <= rmax
  !         = 0              otherwise
  ! see, Baes et al. (2003, MNRAS, 343, 1081)
  ! History:
  !      2015-12-07, bug-fixed
  function rand_rexp(rmax) result(r)
    implicit none
    real(dp), intent(in) :: rmax
    integer, parameter   :: n=2001
    real(dp), save       :: rarr(n), parr(n), dparr, qscale
    real(dp), save       :: rmax_save = -999.9_dp
    real(dp), parameter  :: exp1 = exp(1.0_dp)
    integer  :: i
    real(dp) :: r, p, z
    !$OMP THREADPRIVATE(rarr,parr,dparr,qscale,rmax_save)
    if (rmax_save /= rmax) then
       dparr   = 1.0_dp/(n-1.0_dp)
       parr(:) = (/ ((i-1.0_dp)*dparr, i=1,n) /)
       qscale  = 1.0_dp - (1.0_dp + rmax)*exp(-rmax)
       do i=2,n-1
          z = (parr(i)*qscale-1.0_dp)/exp1
          rarr(i) = -1.0_dp - halley_iteration(z)
       enddo
       ! do not use halley iteraction for calculation at the boundaries.
       ! at the boundaries, it can give NaN. (2015-12-07)
       rarr(1)   = 0.0_dp
       rarr(n)   = rmax
       rmax_save = rmax
    endif

    p = rand_number()
    i = floor(p/dparr)+1
    if (i < 1) then
       r = rarr(1)
    else if (i >= n) then
       r = rarr(n)
    else if (i==1) then
       r = p/dparr * rarr(i+1)
    else
       !r = (p-parr(i))/dparr * (rarr(i+1)-rarr(i)) + rarr(i)
       r = log(rarr(i+1)/rarr(i))/log(parr(i+1)/parr(i)) * log(p/parr(i)) + log(rarr(i))
       r = exp(r)
    endif
  end function rand_rexp
  !--------------------------
  ! Solver of the Lambert W function with Halley Interation Method.
  ! Input parameter range: -1/exp(1.0) < z < 0.0
  ! The Lambert W function is the inverse function of z = w x exp(w).
  !
  function halley_iteration(z) result(w)
    implicit none
    real(dp), intent(in) :: z
    real(dp), parameter :: huge_w=-huge(0.0_dp)
    real(dp) :: expw, wexpw, del, wplus1, w

    integer, parameter :: nc = 6
    real(dp), parameter :: exp1=exp(1.0_dp), exp2=-1.0_dp/exp1, eps_sp=epsilon(0.0_sp)
    real(dp), save :: c(nc) = (/ -221d0/8505d0, 769d0/17280d0, -43d0/540d0, &
                                      11d0/72d0, -1d0/3d0, 1d0 /)
    real(dp) :: p,L1,L2
    integer :: i
    !$OMP THREADPRIVATE(c)

    if (z == 0.0_dp) then
       w = huge_w
    else if (z == exp2) then
       w = -1.0_dp
    else if (z < -0.333_dp) then
       p = -sqrt(2.0_dp*(exp1*z+1.0_dp))
       w = 0.0_dp
       do i=1,nc
         w = (w + c(i))*p
       enddo
       w = w - 1.0_dp
    else if (z >= -0.333_dp .and. z <= -0.033_dp) then
       w = (-8.0960_dp + 391.0025*z)/(1.0_dp-82.9423_dp*z)
    else if (z >= -0.033_dp .and. z < 0.0_dp) then
       L1 = log(-z)
       L2 = log(-L1)
       w  = L1 - L2 + L2/L1 *(1.0_dp + (-2.0_dp+L2)/(2.0_dp*L1))
    endif

    wplus1 = w + 1.0_dp
    if (wplus1 /= 0.0_dp) then
       expw  = exp(w)
       wexpw = w * expw
       del   = wexpw - z
       do while(abs(del) >= eps_sp*abs(z))
          w = w - del/(expw * wplus1 - (wplus1 + 1.0_dp)*del/(2.0_dp*wplus1))
          wplus1 = w + 1.0_dp
          expw  = exp(w)
          wexpw = w * expw
          del   = wexpw - z
       enddo
    endif
  end function halley_iteration
  !-------
  ! Polar Generation of Random Variates with the t-distribution
  !    Bailey (1994, Mathematics of Computation, 62, 779)
  !    2016/02/06, Kwangil Seon, This algorithm is fater than the Kinderman & Monahan's ratio method.
  !-------
  function rand_t(dof) result(t)
  implicit none
  integer, intent(in) :: dof
  real(wp) :: t
  real(wp) :: u, v, w, r, c, u2

  do while(.true.)
     u  = 2.d0*rand_number() - 1.d0
     v  = 2.d0*rand_number() - 1.d0
     u2 = u*u
     w  = u2 + v*v
     !w  = u*u + v*v
     if (w <= 1d0) exit
  end do

  !c = u/sqrt(w)
  !r = sqrt(dof * (w**(-2d0/dof) - 1d0))
  !t = r*c
  c = u2/w
  r = dof * (w**(-2d0/dof) - 1d0)
  t = sqrt(r*c) * sign(1d0, u)
  end function

  !--------------------------------------------------------------
  ! Random Number for Scaled Inverse Chi-square Distribution
  ! 2016/02/16 Kwangil Seon
  function rand_scaled_inv_chi2(dof,var) result(r)
  implicit none
  integer,  intent(in) :: dof
  real(wp), intent(in) :: var
  real(wp) :: r
  r = (dof*var/2d0)/rand_gamma(dof/2d0)
  end function rand_scaled_inv_chi2
  !--------------------------------------------------------------
  !- Gamma Random Deviates for real alpha
  !- combination of two algorithms:
  !      (1) Numerical Recipes for integer alpha and
  !      (2) Martino & Luengo (2013, arXiv:1304.3800v3) for float alpha
  !- 2016/02/15 Kwangil Seon
  function gamdev_r(a) result(x)
  implicit none
  real(wp), intent(in) :: a
  real(wp) :: x
  real(wp) :: kp, bp, ratio
  integer  :: ap, i

  if (a < 1d0) then
     write(*,*) 'alpha should be larger than or equal to 1'
     stop
  endif

  ap = int(a)
  if (ap == a) then
     x = gamdev(ap)
  else
     if (a >= 1d0 .and. a < 2d0) then
        bp = 1d0/a
        kp = exp(a-1d0)*bp**(a-1d0)
     else
        bp = (ap-1d0)/(a-1d0)
        kp = exp(a-ap)*(1d0/(a-1d0))**(a-ap)
     endif

     do while(.true.)
        x     = gamdev(ap)
        ratio = kp * x**(a-ap)*exp(-(1d0-bp)*x)
        if (ratio >= rand_number()) exit
     enddo
  endif
  return
  end function gamdev_r

  !---------------------------------------------------------
  !- Gamma Random Deviates for integer alpha
  !- from Numerical Recipes
  !- 2016/02/15, Kwangil Seon
  !-      slightly modified to make it faster (ifort -fast)
  function gamdev(ia)
  implicit none
  integer,            intent(in) :: ia
  real(wp) :: gamdev
  real(wp) :: am,e,h,s,x,y,v(2),arr(5)
  integer :: i
  if (ia < 1) then
     write(*,*) 'ia should be larger than or equal to 1.'
     stop
  endif

  !if (ia < 6) then
    !call random_number(arr(1:ia))
    !x = -log(product(arr(1:ia)))
  if (ia <= 10) then
    x = 1d0
    do i=1,ia
       x = x*rand_number()
    enddo
    x  = -log(x)
  else
    do
      call random_number(v)
      v(2) = 2.0d0*v(2)-1.0d0
      if (dot_product(v,v) > 1.0d0) cycle
      y  = v(2)/v(1)
      am = ia-1
      s  = sqrt(2.0d0*am+1.0d0)
      x  = s*y+am
      if (x <= 0.0d0) cycle
      e = (1.0d0+y**2)*exp(am*log(x/am)-s*y)
      call random_number(h)
      if (h <= e) exit
    end do
  end if
  gamdev = x
  end function gamdev
  !==============================
  subroutine random_mvnorm(chole, x, dof_t)
  ! generates an n multivariate random normal vector using a Cholesky decomposition.
  !    Dagpunar, J. (1988) Principles of random variate generation
  ! input:
  !    n     = number of variates in vector
  !    avg   = vector of means
  !    chole = Cholescky decomposition calculated with make_mvcovmat2
  ! output:
  !     x    = multivariate random normal vector
  ! 2016-02-20, Kwangil Seon
  implicit none
  real(wp), intent(in)  :: chole(:)
  real(wp), intent(out) :: x(:)
  integer,optional,intent(in) :: dof_t
  real(wp) :: y

  ! local variables
  integer  :: i, j, n, n2, mvt_dof

  n = size(x)
  if (n*(n+1)/2 /= size(chole)) then
     write(*,*) 'something wrong in random_mvnorm'
     stop
  endif

  mvt_dof = -1
  if (present(dof_t)) then
     if (dof_t > 0) mvt_dof = dof_t
  endif

  n2   = 2*n
  x(:) = 0.0_wp
  do j = 1,n
    if (mvt_dof > 0) then
       y = rand_t(mvt_dof)
    else
       y = rand_gauss()
    endif
    do i = j, n
      x(i) = x(i) + chole((j-1)*(n2-j)/2 + i) * y
    enddo
  enddo
  end subroutine random_mvnorm

  subroutine make_mvcovmat2(x, avg, cov, chole)
  ! setup Cholesky decomposition matrix for a given set of input data x(1:n,1:np)
  ! input:
  !    n           = number of data sample
  !    np          = number of parameters
  !    x(1:n,1:np) = input data
  ! output:
  !    avg(1:np)            = vector of means
  !    chole(1:np*(np+1)/2) = lower triangular decomposition of variance matrix
  ! 2016-02-20, Kwangil Seon
  implicit none
  real(wp), intent(in)  :: x(:,:)
  real(wp), intent(out) :: avg(:), cov(:), chole(:)

  integer  :: n, np, ncov
  integer  :: i, j, k, m, n2
  real(wp) :: y, v

  n    = size(x,1)
  np   = size(x,2)
  ncov = np*(np+1)/2
  if (n == 1 .or. np /= size(avg) .or. ncov /= size(cov) .or. ncov /= size(chole)) then
     write(*,*) 'something wrong in make_mvcovmat2'
     stop
  endif

  ! Calculate the mean values
  do j=1, np
    avg(j) = sum(x(:,j))/real(n, kind=wp)
  enddo

  ! Calculate the upper triangle of the covariance matrix
  ! covar(i+j*(j-1)/2) = (i,j)th element (j >= i)
  do j=1,np
     do i=1,j
       k = i + j*(j-1)/2
       cov(k) = sum((x(:,i)-avg(i)) * (x(:,j)-avg(j))) / real(n-1, kind=wp)
     enddo
  enddo

  ! Calculate the lower triangle of the Cholesky decomposition of the covariance matrix
  ! chole((j-1)*(2*np-j)/2+i) = (i, j)th element of lower triangular decomposition of variance matrix (i >= j)
  chole(:) = 0.0d0
  n2       = 2*np
  chole(1) = sqrt(cov(1))
  y        = 1d0/chole(1)
  do j = 2,np
     chole(j) = cov(1+j*(j-1)/2) * y
  enddo

  do i = 2, np
     v = cov(i*(i-1)/2+i)
     do m = 1,i-1
        v = v - chole((m-1)*(n2-m)/2+i)**2
     enddo
     if (v < 0d0) return
     v = sqrt(v)
     y = 1d0/v
     chole((i-1)*(n2-i)/2+i) = v
     do j = i+1, np
        v = cov(j*(j-1)/2+i)
        do m = 1, i-1
           v = v - chole((m-1)*(n2-m)/2+i)*chole((m-1)*(n2-m)/2 + j)
        enddo
        chole((i-1)*(n2-i)/2 + j) = v*y
     enddo
  enddo
  end subroutine make_mvcovmat2

  subroutine make_mvcovmat(x,avg,cov,chole,update)
  ! setup Cholesky decomposition matrix for a given set of input data x(1:n,1:np)
  ! input:
  !    x(1:np) = input data
  ! output:
  !    avg(1:np)            = vector of means
  !    cov(1:np*(np+1)/2)   = upper triangular of covariance matrix
  !    chole(1:np*(np+1)/2) = lower triangular decomposition of variance matrix
  ! 2016-02-20, Kwangil Seon
  implicit none
  real(wp), intent(in)    :: x(:)
  real(wp), intent(inout) :: avg(:), cov(:), chole(:)
  logical,optional,intent(in) :: update

  integer  :: np, ncov
  integer  :: i, j, k, m, n2
  real(wp) :: y, v
  integer,  save :: ncall_covmat = 0, nupdate = 0
  real(wp), save, allocatable :: avg_old(:), cov_old(:)

  np   = size(x)
  ncov = np*(np+1)/2
  if (np /= size(avg) .or. ncov /= size(cov) .or. ncov /= size(chole)) then
     write(*,*) 'something wrong in make_mvcovmat2'
     stop
  endif

  if (present(update)) then
      if (update) nupdate = nupdate + 1
  endif

  if (nupdate == 1) then
     if (allocated(avg_old)) deallocate(avg_old)
     if (allocated(cov_old)) deallocate(cov_old)
     allocate(avg_old(np))
     allocate(cov_old(ncov))
     avg_old = avg
     cov_old = cov
     ncall_covmat = 1
  endif

  ncall_covmat = ncall_covmat + 1
  if (ncall_covmat == 1 .and. nupdate == 0) then
     if (allocated(avg_old)) deallocate(avg_old)
     if (allocated(cov_old)) deallocate(cov_old)
     allocate(avg_old(np))
     allocate(cov_old(ncov))
     avg = x
     cov = epsilon(0.0_wp)
  else
     ! Calculate the mean values
     avg(:) = (ncall_covmat-1.0_wp)/real(ncall_covmat,kind=wp) * avg_old(:) + x(:)/real(ncall_covmat,wp)

     ! Calculate the upper triangle of the covariance matrix
     ! covar(i+j*(j-1)/2) = (i,j)th element (j >= i)
     do j=1,np
        do i=1,j
          k = i + j*(j-1)/2
          cov(k) = (ncall_covmat-2.0_wp)/(ncall_covmat-1.0_wp)*cov_old(k) &
                   + x(i)*x(j)/(ncall_covmat-1.0_wp) &
                   - ncall_covmat/(ncall_covmat-1.0_wp) * avg(i)*avg(j) &
                   + avg_old(i)*avg_old(j)
          if (i == j) cov(k) = cov(k) + epsilon(0.0_wp)
        enddo
     enddo
  endif
  avg_old = avg
  cov_old = cov

  ! Calculate the lower triangle of the Cholesky decomposition of the covariance matrix
  ! chole((j-1)*(2*np-j)/2+i) = (i, j)th element of lower triangular decomposition of variance matrix (i >= j)
  chole(:) = 0.0d0
  n2       = 2*np
  chole(1) = sqrt(cov(1))
  y        = 1d0/chole(1)
  do j = 2,np
     chole(j) = cov(1+j*(j-1)/2) * y
  enddo

  do i = 2, np
     v = cov(i*(i-1)/2+i)
     do m = 1,i-1
        v = v - chole((m-1)*(n2-m)/2+i)**2
     enddo
     if (v < 0d0) return
     v = sqrt(v)
     y = 1d0/v
     chole((i-1)*(n2-i)/2+i) = v
     do j = i+1, np
        v = cov(j*(j-1)/2+i)
        do m = 1, i-1
           v = v - chole((m-1)*(n2-m)/2+i)*chole((m-1)*(n2-m)/2 + j)
        enddo
        chole((i-1)*(n2-i)/2 + j) = v*y
     enddo
  enddo
  end subroutine make_mvcovmat

!==============================
  function rand_permutation(n) result(perm)
  ! (unbiased) random permutation
  ! The "inside-out" shuffle version of the Fisher-Yates shuffle (also known as the Knuth shffle)
     implicit none
     integer, intent(in)   :: n
     integer, dimension(n) :: perm
     integer  :: i, j
     real(wp) :: u

     perm(1) = 1
     do i=2,n
        call random_number(u)
        j       = floor(i*u) + 1
        perm(i) = perm(j)
        perm(j) = i
     enddo
  end function rand_permutation

  function rand_cyclic_permutation(n) result(perm)
  ! Author: Kwang-Il Seon (2016-01-01)
  !         The condition "perm(i) /= i" (i = 1,...,n) is imposed.
  !         This is in fact a cyclic permutation. (Sattolo's algorithm)
     implicit none
     integer, intent(in)   :: n
     integer, dimension(n) :: perm
     integer  :: i, j
     real(wp) :: u

     perm(1) = 2
     if (n >= 2) perm(2) = 1
     do i=3,n
        call random_number(u)
        j       = floor((i-1)*u) + 1
        perm(i) = perm(j)
        perm(j) = i
     enddo
  end function rand_cyclic_permutation
  !------------------------------
  function rand_binomial(pp,ntrial) result(bnldev)
    implicit none
    real(wp), intent(in) :: pp
    integer,  intent(in) :: ntrial
    integer              :: bnldev
    ! local variables
    integer        :: j
    integer, save  :: nold=-1
    real(wp)       :: am,em,g,h,p,sq,t,y,arr(24)
    real(wp), save :: pc,plog,pclog,en,oldg,pold=-1.0
    !$OMP THREADPRIVATE(nold,pc,plog,pclog,en,oldg,pold)

    p  = merge(pp, 1.0_wp-pp, pp <= 0.5_wp )
    am = ntrial*p
    if (ntrial < 25) then
       call random_number(arr(1:ntrial))
       bnldev = count(arr(1:ntrial) < p)
    else if (am < 1.0) then
       g = exp(-am)
       t = 1.0
       do j=0,ntrial
          call random_number(h)
          t = t*h
          if (t < g) exit
       end do
       bnldev = merge(j,ntrial, j <= ntrial)
    else
       if (ntrial /= nold) then
          en   = ntrial
          oldg = gammln(en+1.0_wp)
          nold = ntrial
       end if
       if (p /= pold) then
          pc    = 1.0_wp-p
          plog  = log(p)
          pclog = log(pc)
          pold  = p
       end if
       sq = sqrt(2.0_wp*am*pc)
       do
          call random_number(h)
          y  = tan(PI*h)
          em = sq*y+am
          if (em < 0.0 .or. em >= en+1.0_wp) cycle
          em = int(em)
          t  = 1.2_wp*sq*(1.0_wp+y**2)*exp(oldg-gammln(em+1.0_wp)-&
               gammln(en-em+1.0_wp)+em*plog+(en-em)*pclog)
          call random_number(h)
          if (h <= t) exit
       end do
       bnldev = em
    end if
    if (p /= pp) bnldev = ntrial-bnldev
    contains
     !------------------------------
     function gammln(xx)
       implicit none
       real(wp), intent(in) :: xx
       real(dp) :: gammln
       real(dp) :: tmp,x,ser,y
       real(dp) :: stp = 2.5066282746310005_dp
       real(dp), dimension(6) :: coef = (/76.18009172947146_dp,&
                                         -86.50532032941677_dp,&
                                          24.01409824083091_dp,&
                                          -1.231739572450155_dp,&
                                           0.1208650973866179e-2_dp,&
                                          -0.5395239384953e-5_dp/)
       integer :: i
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
    end function gammln
  end function rand_binomial
  !------------------------------
  function rand_multinomial(p,ntrial) result(ix)
  ! Generates a multinomial random deviate.
  !    P(1:NP)  : the probability that an event will be classified into category I=1,...,NP.
  !    NTRIAL   : number of trials. (0 <= NTRIAL)
  !    IX(1:NP) : a random observation from the multinomial distribution. IX(:) >= 0 and sum(IX) = NTRIAL.
  !
    implicit none
    integer,  intent(in) :: ntrial
    real(wp), intent(in) :: p(:)
    integer :: ix(size(p))

    ! local variables
    integer  :: i, ntot, np
    real(wp) :: prob, ptot
    !integer :: rand_binomial

    !  Initialize variables.
    np    = size(p)
    ntot  = ntrial
    ptot  = 1.0_wp
    ix(:) = 0

    ! Generate the observation.
    do i = 1, np - 1
      prob  = p(i) / ptot
      ix(i) = rand_binomial(prob,ntot)
      ntot  = ntot - ix(i)
      if (ntot <= 0) then
        return
      end if
      ptot = ptot - p(i)
    end do
    ix(np) = ntot

    return
  end function rand_multinomial
  !------------------------------
  function rand_index1(p) result(ix)
    implicit none
    real(wp), intent(in) :: p(:)
    integer :: ix
    integer :: i, np
    integer :: ind(size(p))
  
    np = size(p)
    if (np == 1) then
       ix = 1
    else
       ind = rand_multinomial(p,1)
       do i=1,np
          if (ind(i) == 1) then
             ix = i
             exit
          endif
       enddo
    endif
    return
  end function rand_index1
  function rand_indexn(p,ntrial) result(ix)
    implicit none
    integer,  intent(in) :: ntrial
    real(wp), intent(in) :: p(:)
    integer :: ix(ntrial)
    integer :: i,j,np
    integer :: ind(size(p))
  
    np = size(p)
    if (np == 1) then
       ix(:) = 1
    else
       ind = rand_multinomial(p,ntrial)
       j = 1
       do i=1,np
          if (ind(i) >= 1) then
             ix(j:j+ind(i)-1) = i
             j = j + ind(i)
          endif
       enddo
       ix = ix(rand_permutation(ntrial))
    endif
    return
  end function rand_indexn
  !------------------------------
end module random
