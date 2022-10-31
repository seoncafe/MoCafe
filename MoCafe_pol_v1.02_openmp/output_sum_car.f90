module output_sum
  use define
  use memory_mod
contains
!--------------------------------------------------------------
  subroutine output_normalize(grid)
  implicit none
  type(grid_type), intent(inout) :: grid
  !--- local variable
  real(kind=wp) :: scale_factor

  par%nscatt_tot = par%nscatt_tot /par%no_photons
  scale_factor   = par%no_photons*observer%steradian_pix * par%distance2cm**2

  observer%scatt = observer%scatt / scale_factor
  observer%direc = observer%direc / scale_factor

  if (par%use_stokes) then
     observer%I = observer%I / scale_factor
     observer%Q = observer%Q / scale_factor
     observer%U = observer%U / scale_factor
     observer%V = observer%V / scale_factor
  endif
  end subroutine output_normalize
!--------------------------------------------------------------
  subroutine make_radial_profile(image, prof)
  implicit none
  real(kind=wp), intent(in)  :: image(:,:)
  real(kind=wp), allocatable, intent(out) :: prof(:)

  integer       :: nx, ny
  real(kind=wp) :: xcen, ycen
  integer       :: i,j,ir,nr
  real(kind=wp) :: xx,yy,rr,dr,rmin,rmax
  integer, allocatable :: ncount(:)

  nx   = size(image(:,1))
  ny   = size(image(1,:))
  xcen = (nx+1)/2.0_wp
  ycen = (ny+1)/2.0_wp
  nr   = int(maxval([xcen,ycen]))
  rmax = nr
  if ((nr/2)*2 == nr) then
     dr   = rmax/nr
     rmin = 0.0_wp
  else
     dr   = rmax/(nr-0.5_wp)
     rmin = -dr/2.0_wp
  endif

  if (.not. allocated(prof))   allocate(prof(nr))
  prof(:)   = 0.0_wp
  if (.not. allocated(ncount)) allocate(ncount(nr))
  ncount(:) = 0

  do j=1,ny
    yy = j - ycen
    do i=1,nx
      xx = i - xcen
      rr = sqrt(xx**2 + yy**2)
      ir = floor((rr-rmin)/dr) + 1
      if (ir >= 1 .and. ir <= nr) then
         ncount(ir) = ncount(ir) + 1
         prof(ir)   = prof(ir)   + image(i,j)
      endif
    enddo
  enddo

  where(ncount > 0)
     prof(:) = prof(:)/ncount(:)
  endwhere
  if (allocated(ncount)) deallocate(ncount)
  end subroutine make_radial_profile
!--------------------------------------------------------------
end module output_sum
