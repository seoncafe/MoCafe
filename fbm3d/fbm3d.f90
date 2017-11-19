   subroutine fbm3d(iseed,bvalue,mach,slope_ln,slope_gauss,sigma_lnrho,nx,ny,nz,outfile,out_mode)
!-------------------------------------------------
! Generate Fractal Density in 3D using fractional Brownian motion algorithm.
! 2009/04, Kwangil Seon
!
! 2013-09-17, slightly modified.
!-------------------------------------------------

   use define, only : twopi, wp, sp, i8b
   use random
   implicit none

   integer, intent(in) :: nx,ny,nz,out_mode
   integer, intent(inout) :: iseed
   real, intent(in)  :: bvalue,mach,slope_ln,slope_gauss,sigma_lnrho
   character(len=100), intent(in) :: outfile

   integer, parameter :: FFTW_ESTIMATE = 64
   integer(kind=i8b) :: plan

   complex, allocatable, dimension(:,:,:) :: Ak
   real, allocatable, dimension(:,:,:) :: arr,phi,ang,gauss_k,Pk1,Pk2

   real :: kx(nx/2+1),ky(ny),kz(nz),kr,kxc,kyc,kzc,Anorm,kscale,stdev,mean
   integer :: i,j,k

   call init_random_seed(iseed)

   kxc = nx/2.0
   kyc = ny/2.0
   kzc = nz/2.0
   kscale = 1.0/(nx*ny*nz)

   do i=1,nx/2+1
     kx(i) = i-1
   enddo
   do j=1,ny/2+1
     ky(j) = j-1
   enddo
   do k=1,nz/2+1
     kz(k) = k-1
   enddo
   do j=1,ny/2-1
     ky(ny/2+1+j) = -kyc+j
   enddo
   do k=1,nz/2-1
     kz(nz/2+1+k) = -kyc+k
   enddo

   allocate(phi(nx,ny,nz),ang(nx/2+1,ny,nz))
   do i=1,nx
   do j=1,ny
   do k=1,nz
     phi(i,j,k) = rand_number()
   enddo
   enddo
   enddo

   ang(2:nx/2+1,2:ny,2:nz) = phi(2:nx/2+1,2:ny,2:nz)-phi(nx:nx/2+1:-1,ny:2:-1,nz:2:-1)
   ang(1,2:ny,2:nz)        = phi(1,2:ny,2:nz)-phi(1,ny:2:-1,nz:2:-1)
   ang(2:nx/2+1,1,2:nz)    = phi(2:nx/2+1,1,2:nz)-phi(nx:nx/2+1:-1,1,nz:2:-1)
   ang(2:nx/2+1,2:ny,1)    = phi(2:nx/2+1,2:ny,1)-phi(nx:nx/2+1:-1,ny:2:-1,1)
   ang(1,1,2:nz)           = phi(1,1,2:nz)-phi(1,1,nz:2:-1)
   ang(1,2:ny,1)           = phi(1,2:ny,1)-phi(1,ny:2:-1,1)
   ang(2:nx/2+1,1,1)       = phi(2:nx/2+1,1,1)-phi(nx:nx/2+1:-1,1,1)
   ang(1,1,1)              = 0.0
   deallocate(phi)
   ang = twopi*ang

   allocate(Ak(nx/2+1,ny,nz))
   do i=1,nx/2+1
   do j=1,ny
   do k=1,nz
     kr = sqrt(kx(i)**2+ky(j)**2+kz(k)**2)
     if (kr.gt.0.0) then
       Ak(i,j,k) = cmplx(cos(ang(i,j,k)),sin(ang(i,j,k)))*kr**(-slope_gauss/2.0)
     endif
   enddo
   enddo
   enddo
   deallocate(ang)

! Ak(1,1,1) -> average value of the resulting field
   Ak(1,1,1) = cmplx(0.0,0.0)

!----------------------------------------------------------------
   Anorm = 2.0*sum(dble(abs(Ak(2:nx/2,:,:))**2)) &
              +sum(dble(abs(Ak(1,:,:))**2))+sum(dble(abs(Ak(nx/2+1,:,:))**2))
   write(*,*) 'K_norm',Anorm
   Ak = Ak/sqrt(Anorm)

!----------------------------------------------------------------
   allocate(gauss_k(nx/2+1,ny,nz))
   do i=1,nx/2+1
   do j=1,ny
   do k=1,nz
      gauss_k(i,j,k) = rand_gauss()
   enddo
   enddo
   enddo
   Ak(:,:,:) = Ak(:,:,:) * gauss_k(:,:,:)
   deallocate(gauss_k)

!----------------------------------------------------------------
   allocate(Pk1(nx/2+1,ny,nz),Pk2(nx/2+1,ny,nz))
   Pk1 = abs(Ak(:,:,:))**2

!----------------------------------------------------------------
   allocate(arr(nx,ny,nz))
   call sfftw_plan_dft_c2r_3d(plan,nx,ny,nz,Ak,arr,FFTW_ESTIMATE)
   call sfftw_execute_dft_c2r(plan,Ak,arr)
   call sfftw_destroy_plan(plan)

   mean  = sum(arr)/size(arr)
   stdev = sqrt(sum((arr-mean)**2)/size(arr))
!   write(*,*) 'MEAN, STDEV',mean,stdev
   arr   = (arr - mean)/stdev
   mean  = sum(arr)/size(arr)
   stdev = sqrt(sum((arr-mean)**2)/size(arr))
   write(*,*) 'MEAN, STDEV',mean,stdev

! Transform to lognormal distribution
   arr = exp(arr*sigma_lnrho)

!----------------------------------------------------------------
   if (out_mode == 3) then
      call sfftw_plan_dft_r2c_3d(plan,nx,ny,nz,arr,Ak,FFTW_ESTIMATE)
      call sfftw_execute_dft_r2c(plan,arr,Ak)
      call sfftw_destroy_plan(plan)
      Ak = Ak*kscale
      Pk2 = abs(Ak(:,:,:))**2
   endif

!----------------------------------------------------------------
   call write_output(outfile,out_mode,iseed,bvalue,mach,slope_ln,slope_gauss,sigma_lnrho,nx,ny,nz,arr,Pk1,Pk2)
   deallocate(arr)
   deallocate(Pk1,Pk2)
   deallocate(Ak)

   end subroutine fbm3d
