  program make_fbm3d
! updated slightly, 2017-11-18
  implicit none
  integer :: unit
  integer :: iseed,nx,ny,nz,out_mode
  real :: slope_ln,slope_gauss,sigma_lnrho
  character(len=100) :: model_infile
  character(len=100) :: outfile

  real :: bvalue,mach
  namelist /input/ iseed,bvalue,mach,nx,outfile,out_mode

  real :: p(5,4) = (/2.841e-01, -9.168e-01, -9.334e-01,  1.221e+00, -2.546e-01, &
                     8.173e-01,  5.994e-01,  1.326e+00, -1.125e+00,  2.119e-01, &
                     5.019e-02, -1.468e-01, -3.838e-01,  2.970e-01, -5.417e-02, &
                    -4.428e-03,  1.186e-02,  3.111e-02, -2.375e-02,  4.329e-03/)
  real :: bcoeff(4)
  integer :: i,j

! Comments added on 2013-09-17, Kwang-Il Seon
!
! bvalue = 1/3, 0.4, 1.0 for solenoidal, natural mixing, and compressive modes
! sigma_lnrho = sqrt(ln(1+(sigma_rho)^2)) = sqrt(1+bvalue^2 Mach^2)
! To generate 3D lognormal density field, we need two parameters
!   (1) sigma_lnrho : standard deviation of ln(rho), or equivalently standard deviation of Gaussian random field.
!   (2) slope_gauss : power-law indenx of power spectrum (or spectral slope) of Gaussian random field.
!
! See Seon (2012, ApJL, 761, L17) and Seon & Draine (2016, ApJ, 833, 201)

! default input parameters.
  iseed    = 0
  nx       = 128
  bvalue   = 0.4d0
  mach     = 1.0d0
  outfile  = 'M010b04.fits.gz'
  out_mode = 1

! read input parameters
  !read(5,input)
  if (command_argument_count() >= 1) then
     call get_command_argument(1, model_infile)
  else
     write(*,*) 'Usage: make_fbm3d input_file'
     write(*,*) '       see M005b040_001.in for a sample input'
     stop
  endif
  ! newunit specifier is instrodueced in Fortran 2008.
  open(newunit=unit,file=trim(model_infile),status='old')
  read(unit,input)
  close(unit)

! bvalue = 1/3, 0.4, or 1.0
  if (bvalue < 0.4) then
     bvalue = 1.0/3.0
  else if (bvalue > 0.4) then
     bvalue = 1.0
  endif
  if (bvalue /= 1.0/3.0 .and. bvalue /= 0.4 .and. bvalue /= 1.0) then
     write(*,*) 'bvalue should be 1/3, 0.4, or 1.0.'
     write(*,*) 'The input bvalue was ',bvalue
     stop
  endif

  sigma_lnrho = sqrt(alog(real(1.0d0+(1d0*bvalue*mach)**2)))
  if (bvalue <= 0.4) then
     slope_ln = 3.81d0*mach**(-0.16d0)
  else
     slope_ln = 3.81d0*mach**(-0.16d0) + 0.6d0
  endif

  slope_gauss = 0.0d0
  bcoeff(:)   = 0.0d0
  do j=1,4
     do i=1,5
        bcoeff(j) = bcoeff(j)  + p(i,j)    * sigma_lnrho**(i-1)
     enddo
     slope_gauss = slope_gauss + bcoeff(j) * (slope_ln)**(j-1)
  enddo

  ny = nx
  nz = nx
  print*,'B-value                : ',bvalue
  print*,'Mach                   : ',mach
  print*,'Power-law Index(Ln)    : ',slope_ln
  print*,'Power-law Index(Gauss) : ',slope_gauss
  print*,'Sigma_Ln(rho)          : ',sigma_lnrho

  call fbm3d(iseed,bvalue,mach,slope_ln,slope_gauss,sigma_lnrho,nx,ny,nz,outfile,out_mode)
  end program make_fbm3d
