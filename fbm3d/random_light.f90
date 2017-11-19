module random
!-----
!--- Random Number Generator Module
!--- Use Mersenne Twister (MT) random generator
!
! Author: Kwang-Il Seon
!
! Note:
!    - gfortran random_number uses KISS algorithm while ifort uses Lecuyer algorithm.
!    - intrinsic random_number routine is not thread-safe with openmp.
!    - the intial seed in each thread should be different from thread to thread.
!
! Usage:
!    call init_random_seed()     to initialize random seed
!    call random_number(harvest) to call as a subroutine
!    harvest = rand_number()     to call as a function
!         where harvest is a number (32- or 64-bit) or an 1-, 2-, or 3-dimentional array.
!
! History:
!   v1.00 (2017/11/18) first light version.
!-----
  use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64
  implicit none
  private
  integer, parameter  :: sp = real32
  integer, parameter  :: dp = real64
  integer, parameter  :: wp = dp
  real(wp), parameter :: PI = 3.141592653589793238462643383279502884197_wp
  real(wp), parameter :: TWOPI  = PI*2.0_wp
  real(wp), parameter :: HALFPI = PI/2.0_wp
  real(wp), parameter :: FOURPI = PI*4.0_wp

  public init_random_seed
  public random_seed, random_number, random_gauss
  public rand_number, rand_gauss

  !======= parameter for Mersenne Twister Random Number Generator (MT19937)
  ! Period parameters
  INTEGER, PARAMETER :: mt_n = 624, mt_n1 = mt_n+1, mt_m = 397, mt_mata = -1727483681 ! constant vector a
  INTEGER, PARAMETER :: mt_seed0 = 4357
  INTEGER, SAVE      :: mt(0:mt_n-1)     ! the array for the state vector
  INTEGER, SAVE      :: mti = mt_n1      ! mti == N+1 means mt(:) is not initialized
  !$OMP THREADPRIVATE(mt, mti)

  interface init_random_seed
     module procedure init_random_mt_seed
  end interface init_random_seed

  interface random_seed
     module procedure random_mt_seed
  end interface random_seed

  interface random_number
     module procedure random_mt32,random_mt64
  end interface random_number
  !=================================================================================

  interface random_gauss
     module procedure random_gauss1
  end interface random_gauss

  interface rand_gauss
     module procedure rand_gauss1
  end interface rand_gauss
contains
!==================================================================
  subroutine random_mt32(harvest)
     implicit none
     real(sp), intent(out) :: harvest
     harvest = rand_number()
  end subroutine random_mt32
  subroutine random_mt64(harvest)
     implicit none
     real(dp), intent(out) :: harvest
     harvest = rand_number()
  end subroutine random_mt64
  !-----------------------------------------------------------------------
  !   genrand() generates one pseudorandom real number (double) which is
  ! uniformly distributed on [0,1]-interval, for each call.
  ! sgenrand(seed) set initial values to the working area of 624 words.
  ! Before genrand(), sgenrand(seed) must be called once.  (seed is any 32-bit
  ! integer except for 0).
  !-----------------------------------------------------------------------
  ! Mersenne Twister Random Number Generator (MT19937)
  !   Matsumoto & Nishimura (1998, ACM Transactions of Modeling and Computer Simulation, 8, 3)
  !   Period = 2^19937 - 1 ~ 10^6001
  ! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
  !   genrand()      -> double precision function grnd()
  !   sgenrand(seed) -> subroutine sgrnd(seed)
  !                     integer seed
  !
  ! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
  ! When you use this, send an email to: matumoto@math.keio.ac.jp
  ! with an appropriate reference to your work.
  !-----------------------------------------------------------------------
  SUBROUTINE sgrnd(seed)
  ! This is the original version of the seeding routine.
  ! It was replaced in the Japanese version in C on 26 January 2002
  ! It is recommended that routine init_genrand is used instead.

  INTEGER, INTENT(IN)   :: seed
  !    setting initial seeds to mt[N] using the generator Line 25 of Table 1 in
  !    [KNUTH 1981, The Art of Computer Programming Vol. 2 (2nd Ed.), pp102]

  mt(0) = IAND(seed, -1)
  DO  mti=1,mt_n-1
    mt(mti) = IAND(69069 * mt(mti-1), -1)
  END DO

  RETURN
  END SUBROUTINE sgrnd

  !-----------------------------------------------------------------------
  SUBROUTINE init_mt(seed)
  ! This initialization is based upon the multiplier given on p.106 of the
  ! 3rd edition of Knuth, The Art of Computer Programming Vol. 2.
  ! This version assumes that integer overflow does NOT cause a crash.

  INTEGER, INTENT(IN)  :: seed
  INTEGER  :: latest

  mt(0)  = seed
  latest = seed
  DO mti = 1, mt_n-1
    latest  = IEOR( latest, ISHFT( latest, -30 ) )
    latest  = latest * 1812433253 + mti
    mt(mti) = latest
  END DO

  RETURN
  END SUBROUTINE init_mt
  !-----------------------------------------------------------------------
  FUNCTION rand_number() RESULT(fn_val)
  REAL (dp) :: fn_val

  INTEGER, SAVE :: mag01(0:1) = [0, mt_mata]
  !$OMP THREADPRIVATE(mag01)
  INTEGER       :: kk, y
  ! Period parameters
  INTEGER, PARAMETER :: lmask =  2147483647                            ! least significant r bits
  INTEGER, PARAMETER :: umask = -2147483647 - 1                        ! most significant w-r bits
  ! Tempering parameters
  INTEGER, PARAMETER  :: tmaskb= -1658038656, tmaskc= -272236544
  real(dp), parameter :: a = 1.0_real64/4294967297.0_real64
  real(dp), parameter :: b = 2147483649.0_real64/4294967297.0_real64
  !  ishft(i,n): If n > 0, shifts bits in i by n positions to left.
  !              If n < 0, shifts bits in i by n positions to right.
  !  iand (i,j): Performs logical AND on corresponding bits of i and j.
  !  ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
  !  ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
  !  statement functions
  integer :: tshftu, tshfts, tshftt, tshftl
  tshftu(y) = ISHFT(y,-11)
  tshfts(y) = ISHFT(y,7)
  tshftt(y) = ISHFT(y,15)
  tshftl(y) = ISHFT(y,-18)

  IF(mti >= mt_n) THEN       !  generate N words at one time
    IF(mti == mt_n+1) THEN   !  if sgrnd() has not been called,
      !CALL sgrnd(4357)       !  a default initial seed is used
      CALL init_mt(mt_seed0)  !  a default initial seed is used
    END IF

    DO  kk = 0, mt_n-mt_m-1
      y      = IOR(IAND(mt(kk),umask), IAND(mt(kk+1),lmask))
      mt(kk) = IEOR(IEOR(mt(kk+mt_m), ISHFT(y,-1)),mag01(IAND(y,1)))
    END DO
    DO  kk = mt_n-mt_m, mt_n-2
      y      = IOR(IAND(mt(kk),umask), IAND(mt(kk+1),lmask))
      mt(kk) = IEOR(IEOR(mt(kk+(mt_m-mt_n)), ISHFT(y,-1)),mag01(IAND(y,1)))
    END DO
    y          = IOR(IAND(mt(mt_n-1),umask), IAND(mt(0),lmask))
    mt(mt_n-1) = IEOR(IEOR(mt(mt_m-1), ISHFT(y,-1)),mag01(IAND(y,1)))
    mti        = 0
  END IF

  y   = mt(mti)
  mti = mti + 1
  y = IEOR(y, tshftu(y))
  y = IEOR(y, IAND(tshfts(y),tmaskb))
  y = IEOR(y, IAND(tshftt(y),tmaskc))
  y = IEOR(y, tshftl(y))

  !IF (y < 0) THEN
  !  fn_val = (DBLE(y) + 2.0D0**32) / (2.0D0**32 - 1.0D0)
  !ELSE
  !  fn_val = DBLE(y) / (2.0D0**32 - 1.0D0)
  !END IF
  fn_val = b + a * dble(y)

  RETURN
  END FUNCTION rand_number
!==================================
  subroutine random_mt_seed_(seed_size,put,get)
    ! size is an intrinsic function so that we need to define an intermediate routine
    ! which use the argument name "seed_size" instead of "size"
    ! random_kiss_seed has the optional argument "size"
    ! bug-fixed, 2017-06-25
    implicit none
    integer, optional, intent(out) :: seed_size
    integer, optional, intent(in)  :: put(:)
    integer, optional, intent(out) :: get(:)
    ! local variable
    integer :: n

    if (present(seed_size)) then
       seed_size = mt_n
    endif
    if (present(put)) then
       n              = minval([mt_n,size(put)])
       !seed_kiss(1:n) = put(1:n)
       mt(0:n-1) = put(1:n)
    endif
    if (present(get)) then
       n        = minval([mt_n,size(get)])
       !get(1:n) = seed_kiss(1:n)
       get(1:n) = mt(0:n-1)
    endif
  end subroutine random_mt_seed_
  subroutine random_mt_seed(size,put,get)
    implicit none
    integer, optional, intent(out) :: size
    integer, optional, intent(in)  :: put(:)
    integer, optional, intent(out) :: get(:)
    call random_mt_seed_(seed_size=size,put=put,get=get)
  end subroutine random_mt_seed
!-------------------------------------
  subroutine init_random_mt_seed(iseed)
#ifdef _OPENMP
    use omp_lib
#endif
#ifdef MPI
    use mpi
#endif
    implicit none
    integer, optional :: iseed
    ! local variable
    integer        :: i, n, un, istat, dt(8), pid
    !integer(int64) :: t
    integer        :: t
    integer        :: getpid
    logical        :: set_iseed
#ifdef MPI
    integer        :: myid, ierr
#endif

    set_iseed = .false.
    if (present(iseed)) then
       if (iseed /= 0) then
          set_iseed = .true.
          t         = iseed
        endif
    endif
    if (.not. set_iseed) then
       ! First try if the OS provides a random number generator
       ! return an integer uniform random deviate between [0 and huge(0)] from /dev/urandom if successful.
       open(newunit=un, file="/dev/urandom", access="stream", &
            form="unformatted", action="read", status="old", iostat=istat)
       if (istat == 0) then
          read(un) t
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
       end if
    else
    ! bug-fixed, 2017-06-25. "else" part was missing. I don't know why this part was missing.
       t = iseed
#ifdef _OPENMP
       ! OPENMP does not launch multiple instances.
       t = t + omp_get_thread_num()
#endif
#ifdef MPI
       call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
       t = t + myid
#endif
    endif
    call init_mt(t)
  end subroutine init_random_mt_seed
  !-------------------------------------
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
  !-----------------------------------------------------
end module random
