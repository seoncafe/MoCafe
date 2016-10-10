module grid_mod
  use define, only : wp
  implicit none
  ! note: 1 - exp(-1) = 0.61312
  real(wp), parameter :: fraction_center = 0.6312_wp
  real(wp), parameter :: accuracy        = 1.0e-2_wp
  private fraction_center
contains
!---------------------------------------------------------------------------------
  subroutine grid_create(grid)
  use define
  use density_mod
  use omp_lib
#ifdef MPI
  use mpi
#endif
  implicit none
  type(grid_type), intent(inout) :: grid

! local variables
  integer :: i,j,k,id
  integer :: nxcen,nycen,nzcen
  real(kind=wp) :: r,z,opac,taueq,taupole
  real(kind=wp), allocatable :: eq_face(:)

#ifdef MPI
  integer :: ierr,myid
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
#endif
  !--- Set up grid faces.
    grid%nx     = par%nx
    grid%ny     = par%ny
    grid%nz     = par%nz
    grid%z_alpha= par%z_alpha
    grid%xmax   = par%xmax
    grid%ymax   = par%ymax
    grid%zmax   = par%zmax
    grid%xmin   = -par%xmax
    grid%ymin   = -par%ymax
    grid%zmin   = -par%zmax
    grid%xrange = grid%xmax - grid%xmin
    grid%yrange = grid%ymax - grid%ymin
    grid%zrange = grid%zmax - grid%zmin
    nullify(grid%xface)
    nullify(grid%yface)
    nullify(grid%zface)
    nullify(grid%opacity)
    if (.not. associated(grid%xface))   allocate(grid%xface(grid%nx+1))
    if (.not. associated(grid%yface))   allocate(grid%yface(grid%ny+1))
    if (.not. associated(grid%zface))   allocate(grid%zface(grid%nz+1))
    if (.not. associated(grid%opacity)) allocate(grid%opacity(grid%nx,grid%ny,grid%nz))

    nullify(grid%dz)
    if (.not. associated(grid%dz))   allocate(grid%dz(grid%nz))

    grid%dx     = grid%xrange/grid%nx
    grid%dy     = grid%yrange/grid%ny
    !grid%dz     = grid%zrange/grid%nz
    grid%xface(:) = [ ((i-1)*grid%dx + grid%xmin, i=1,grid%nx+1) ]
    grid%yface(:) = [ ((i-1)*grid%dy + grid%ymin, i=1,grid%ny+1) ]
    !grid%zface(:) = [ ((i-1)*grid%dz + grid%zmin, i=1,grid%nz+1) ]

    if (grid%z_alpha == -999.0_wp) then
       call find_alpha(zexp_tint,grid%z_alpha)
       !print*,'grid%z_alpha = ',grid%z_alpha
    endif

    if (grid%z_alpha == 1.0_wp .or. grid%z_alpha == 0.0_wp) then
       grid%dz(:)    = grid%zrange/grid%nz
       grid%zface(:) = [ (dble(i-1)*grid%dz(1) + grid%zmin, i=1,grid%nz+1) ]
    else
       if (allocated(eq_face)) deallocate(eq_face)
       allocate(eq_face(grid%nz+1))
       eq_face(:)    = [ (dble(i-1)*2.0_wp/grid%nz - 1.0_wp, i=1,grid%nz+1) ]
       grid%zface(:) = grid%zmax * (abs(eq_face(:)))**grid%z_alpha * sign(1.0_wp,eq_face(:))
       grid%zmin_eq  = (abs(grid%zmin))**(1.0d0/grid%z_alpha) * sign(1.0_wp, grid%zmin)
       grid%zmax_eq  = (abs(grid%zmax))**(1.0d0/grid%z_alpha) * sign(1.0_wp, grid%zmax)
       grid%dz_eq    = (grid%zmax_eq - grid%zmin_eq)/grid%nz
       grid%dz(:)    = grid%zface(2:grid%nz+1) - grid%zface(1:grid%nz)
       deallocate(eq_face)
    endif

!--- setup dust density
    nullify(grid%opacity)
    if (.not. associated(grid%opacity)) allocate(grid%opacity(grid%nx,grid%ny,grid%nz))
    grid%opacity(:,:,:) = 0.0_wp
!    !$OMP PARALLELDO &
!    !$OMP private(i,j,k,z,r,opac) shared(grid,dust) &
!    !$OMP schedule(static,1)
    do k=1,par%nz
      z = (grid%zface(k)+grid%zface(k+1))/2.0_wp
      do j=1,par%ny
        do i=1,par%nx
          r = sqrt((grid%xface(i)+grid%xface(i+1))**2 + (grid%yface(j)+grid%yface(j+1))**2)/2.0_wp
          do id=1,par%ndust
            call dust_density(par%dust(id),r,z,opac)
            grid%opacity(i,j,k) = grid%opacity(i,j,k) + opac
          enddo
        enddo
      enddo
    enddo
    !open(1,file='rho.dat',form='unformatted')
    !write(1) 3
    !write(1) par%nr,par%np,par%nz
    !write(1) grid%opacity
    !close(1)
!    !$OMP END PARALLELDO
!--- Calculate equatorial and polar optical depths
    nxcen = (grid%nx+1)/2
    nycen = (grid%ny+1)/2
    nzcen   = (grid%nz+1)/2
    taueq   = sum(grid%opacity(:,nycen,nzcen))*2.0_wp*grid%xmax/grid%ny
    taupole = sum(grid%opacity(nxcen,nycen,:))*2.0_wp*grid%zmax/grid%nz
#ifdef MPI
    if (myid==0) then
       write(6,'(a,f8.3)') 'tau_equat  : ',taueq
       write(6,'(a,f8.3)') 'tau_pole   : ',taupole
    endif
#else
    write(6,'(a,f8.3)') 'tau_equat  : ',taueq
    write(6,'(a,f8.3)') 'tau_pole   : ',taupole
#endif
  end subroutine grid_create

  !-----------------
  subroutine grid_destroy(grid)
  use define
  implicit none
  type(grid_type), intent(inout) :: grid
     if (associated(grid%xface))   deallocate(grid%xface)
     if (associated(grid%yface))   deallocate(grid%yface)
     if (associated(grid%zface))   deallocate(grid%zface)
     if (associated(grid%opacity)) deallocate(grid%opacity)
     if (associated(grid%dz))      deallocate(grid%dz)
  end subroutine grid_destroy

  ! calculate power-law index alpha for mesh.
  subroutine find_alpha(func,alpha)
  use define, only : wp
  implicit none
  real(wp), intent(out) :: alpha
  interface
    function func(x)
    use define, only : wp
    implicit none
    real(wp), intent(in) :: x
    real(wp) :: func
    end function func
  end interface
  real(wp) :: a1,a2,da, f1,f2,amid

  a1 = 1.0_wp
  a2 = a1 + 2.0_wp
  f1 = func(a1)
  f2 = func(a2)
  do while(f1*f2 >= 0.0_wp)
     a1 = a2
     a2 = a2 + 2.0_wp
     f1 = f2
     f2 = func(a2)
  enddo

  if (f1 < 0.0) then
    alpha = a1
    da    = a2-a1
  else
    alpha = a2
    da    = a1-a2
  end if

  do while (.true.)
    da   = da*0.5_wp
    amid = alpha + da
    f2   = func(amid)
    if (f2 <= 0.0_wp) alpha = amid
    if (abs(da) < abs(alpha)*accuracy .or. f2 == 0.0_wp) return
  end do
  end subroutine find_alpha

  ! ancilliary function to calculate the power-law index for z-axis grid
  function zexp_tint(alpha) result(tint)
  use define, only : wp, par
  implicit none
  real(wp), intent(in) :: alpha
  real(wp) :: tint
  real(wp) :: ncen, eq_cen, zcen, taumax
  integer  :: i

  ! note:
  ! (ncen, eq_cen) = ((nz/2+1), 0) if fraction_center = 0.0
  ! (ncen, eq_cen) = (nz+1,     1) if fraction_center = 1.0
  ncen   = par%nz/2.0_wp * fraction_center + (par%nz/2.0_wp+1.0_wp)
  !eq_cen = (ncen-(par%nz/2.0_wp+1.0_wp))/(par%nz/2.0_wp)
  eq_cen = fraction_center
  zcen   = par%zmax * (abs(eq_cen))**alpha

  taumax = 0.0_wp
  tint   = 0.0_wp
  do i=1,par%ndust
     taumax=taumax+par%dust(i)%tau_faceon*(1.0_wp-exp(-par%dust(i)%zmax/par%dust(i)%zscale))
     tint  =tint  +par%dust(i)%tau_faceon*(1.0_wp-exp(-zcen/par%dust(i)%zscale))
  enddo
  tint = tint - fraction_center * taumax
  !write(11,'(*(es14.6))') alpha,ncen,eq_cen,zcen,tint
  end function zexp_tint

!  ! ancilliary function to calculate the power-law index for r-axis grid
!  function rexp_tint(alpha) result(tint)
!  use define, only: wp, par
!  implicit none
!  real(wp), intent(in) :: alpha
!  real(wp) :: tint
!  real(wp) :: ncen, eq_cen, rcen, taumax
!  integer  :: i
!
!  ! note:
!  ! (ncen, eq_cen) = (1,    0) if fraction_center = 0.0
!  ! (ncen, eq_cen) = (nr+1, 1) if fraction_center = 1.0
!  ncen   = fraction_center * par%nr + 1.0_wp
!  !eq_cen = (ncen-1.0_wp)/par%nr
!  eq_cen = fraction_center
!  rcen   = par%rmax * (abs(eq_cen))**alpha
!
!  taumax = 0.0_wp
!  tint   = 0.0_wp
!  do i=1,par%ndust
!     taumax=taumax+par%dust(i)%tau_faceon*(par%dust(i)%rscale/par%dust(i)%zscale)*(1.0_wp-exp(-par%dust(i)%rmax/par%dust(i)%rscale))
!     tint  =tint  +par%dust(i)%tau_faceon*(par%dust(i)%rscale/par%dust(i)%zscale)*(1.0_wp-exp(-rcen/par%dust(i)%rscale))
!  enddo
!  tint = tint - fraction_center * taumax
!  !write(12,'(*(es14.6))') alpha,ncen,eq_cen,rcen,tint
!  end function rexp_tint

end module grid_mod
