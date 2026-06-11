module output_sum
  use define
  use mpi
  use memory_mod
contains
!--------------------------------------------------------------
  subroutine output_reduce(grid)
  implicit none
  type(grid_type), intent(inout) :: grid

  !--- local variables
  integer :: ierr, k

  call MPI_ALLREDUCE(MPI_IN_PLACE, par%nscatt_tot,  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

  do k=1,par%nobs
     if (par%use_tau_list) then
        call reduce_mem(observer(k)%scatt_agt)
        call reduce_mem(observer(k)%direc_t)
     else if (par%use_ag_list) then
        call reduce_mem(observer(k)%scatt_ag)
        call reduce_mem(observer(k)%direc)
     else
        call reduce_mem(observer(k)%scatt)
        call reduce_mem(observer(k)%direc)
     endif
     if (par%save_direc0) then
        call reduce_mem(observer(k)%direc0)
     endif
     if (par%use_stokes) then
        call reduce_mem(observer(k)%I)
        call reduce_mem(observer(k)%Q)
        call reduce_mem(observer(k)%U)
        call reduce_mem(observer(k)%V)
     endif
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end subroutine output_reduce
!--------------------------------------------------------------
  subroutine output_normalize(grid)
  implicit none
  type(grid_type), intent(inout) :: grid
  !--- local variable
  real(kind=wp) :: scale_factor, lum_fac
  integer :: k

  if (par%output_normalization(1:9) == 'intensity') then
      select case(trim(par%source_geometry))
      case ('external_sph')
         lum_fac        = fourpi*par%rmax**2 * pi
         par%luminosity = lum_fac * par%distance2cm**2
      case ('external_cyl')
         lum_fac        = (2.d0*pi*par%rmax**2 + 2.d0*pi*par%rmax*2.d0*par%zmax) * pi
         par%luminosity = lum_fac * par%distance2cm**2
      case ('external_rec')
         lum_fac        = 8.d0* (par%xmax*par%ymax + par%ymax*par%zmax + par%zmax*par%xmax) * pi
         par%luminosity = lum_fac * par%distance2cm**2
      case default
         lum_fac        = 1.d0
         par%output_normalization = 'luminosity'
      end select
  endif

  par%nscatt_tot = par%nscatt_tot /par%nphotons
  do k=1, par%nobs
     scale_factor   = par%no_photons*observer(k)%steradian_pix * par%distance2cm**2 / par%luminosity
     !if (par%source_geometry(1:8) == 'external' .and. par%output_normalization(1:9) == 'intensity') then
     !   scale_factor = scale_factor / lum_fac
     !endif

     if (par%use_tau_list) then
        observer(k)%scatt_agt = observer(k)%scatt_agt / scale_factor
        observer(k)%direc_t   = observer(k)%direc_t   / scale_factor
     else if (par%use_ag_list) then
        observer(k)%scatt_ag = observer(k)%scatt_ag / scale_factor
        observer(k)%direc    = observer(k)%direc / scale_factor
     else
        observer(k)%scatt = observer(k)%scatt / scale_factor
        observer(k)%direc = observer(k)%direc / scale_factor
     endif
     if (par%save_direc0) then
        observer(k)%direc0 = observer(k)%direc0 / scale_factor
     endif

     if (par%use_stokes) then
        observer(k)%I = observer(k)%I / scale_factor
        observer(k)%Q = observer(k)%Q / scale_factor
        observer(k)%U = observer(k)%U / scale_factor
        observer(k)%V = observer(k)%V / scale_factor
     endif
  enddo
  !if (observer%nxim == observer%nyim) call make_radial_profile(observer%I,observer%Iprof)
  end subroutine output_normalize
!--------------------------------------------------------------
!  subroutine make_radial_profile(image, prof)
!  implicit none
!  real(kind=wp), intent(in)  :: image(:,:)
!  real(kind=wp), allocatable, intent(out) :: prof(:)
!
!  integer       :: nx, ny
!  real(kind=wp) :: xcen, ycen
!  integer       :: i,j,ir,nr
!  real(kind=wp) :: xx,yy,rr,dr,rmin,rmax
!  integer, allocatable :: ncount(:)
!
!  nx   = size(image(:,1))
!  ny   = size(image(1,:))
!  xcen = (nx+1)/2.0_wp
!  ycen = (ny+1)/2.0_wp
!  nr   = int(maxval([xcen,ycen]))
!  rmax = nr
!  if ((nr/2)*2 == nr) then
!     dr   = rmax/nr
!     rmin = 0.0_wp
!  else
!     dr   = rmax/(nr-0.5_wp)
!     rmin = -dr/2.0_wp
!  endif
!
!  if (.not. allocated(prof))   allocate(prof(nr))
!  prof(:)   = 0.0_wp
!  if (.not. allocated(ncount)) allocate(ncount(nr))
!  ncount(:) = 0
!
!  do j=1,ny
!    yy = j - ycen
!    do i=1,nx
!      xx = i - xcen
!      rr = sqrt(xx**2 + yy**2)
!      ir = floor((rr-rmin)/dr) + 1
!      if (ir >= 1 .and. ir <= nr) then
!         ncount(ir) = ncount(ir) + 1
!         prof(ir)   = prof(ir)   + image(i,j)
!      endif
!    enddo
!  enddo
!
!  where(ncount > 0)
!     prof(:) = prof(:)/ncount(:)
!  endwhere
!  if (allocated(ncount)) deallocate(ncount)
!  end subroutine make_radial_profile
!--------------------------------------------------------------
end module output_sum
