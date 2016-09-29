module detect
contains
  subroutine add_direct(photon,grid,observer,output)

  use define
  use raytrace, only : raytrace_to_edge_cyl
  implicit none

  type(photon_type),    intent(in) :: photon
  type(grid_type),      intent(in) :: grid
  type(observer_type),  intent(in) :: observer
  type(output_type), intent(inout) :: output

! local variables
  type(photon_type) :: pobs
  real(kind=wp) :: r2,r,vdet(3),wgt,tau
  integer :: ix,iy
  integer :: i

  pobs    = photon
  pobs%vx = (observer%x-photon%x)
  pobs%vy = (observer%y-photon%y)
  pobs%vz = (observer%z-photon%z)
  r2      = pobs%vx*pobs%vx + pobs%vy*pobs%vy + pobs%vz*pobs%vz
  r       = sqrt(r2)
  pobs%vx = pobs%vx/r
  pobs%vy = pobs%vy/r
  pobs%vz = pobs%vz/r

  call raytrace_to_edge_cyl(pobs,grid,tau)
  wgt = exp(-tau)/(fourpi*r2) * photon%lscale

  do i=1,3
     vdet(i) = observer%rmatrix(i,1) * pobs%vx + &
               observer%rmatrix(i,2) * pobs%vy + &
               observer%rmatrix(i,3) * pobs%vz
  enddo

  ix = floor(atan2(-vdet(1),vdet(3))*rad2deg/output%dx+output%nx/2.0_wp + par%xshift) + 1
  iy = floor(atan2(-vdet(2),vdet(3))*rad2deg/output%dy+output%ny/2.0_wp + par%yshift) + 1

  if (ix >= 1 .and. ix <= output%nx .and. iy >= 1 .and. iy <= output%ny) then
     if (par%output_mode == 1) output%tot(ix,iy)  =output%tot(ix,iy)  +wgt
     if (par%output_mode == 2) output%direc(ix,iy)=output%direc(ix,iy)+wgt
  endif

  !----------
  ! the model is axi-symmetric. (2016-03-19)
  if (par%left_right_fold .and. par%phase_angle == 0.0_wp) then
     pobs%vx = (observer%x+photon%x)
     pobs%vy = (observer%y-photon%y)
     pobs%vz = (observer%z-photon%z)
     !r2      = pobs%vx*pobs%vx + pobs%vy*pobs%vy + pobs%vz*pobs%vz
     !r       = sqrt(r2)
     !pobs%vx = pobs%vx/r
     !pobs%vy = pobs%vy/r
     !pobs%vz = pobs%vz/r

     do i=1,3
        vdet(i) = observer%rmatrix(i,1) * pobs%vx + &
                  observer%rmatrix(i,2) * pobs%vy + &
                  observer%rmatrix(i,3) * pobs%vz
     enddo
     ix = floor(atan2(-vdet(1),vdet(3))*rad2deg/output%dx+output%nx/2.0_wp + par%xshift) + 1
     iy = floor(atan2(-vdet(2),vdet(3))*rad2deg/output%dy+output%ny/2.0_wp + par%yshift) + 1

     if (ix >= 1 .and. ix <= output%nx .and. iy >= 1 .and. iy <= output%ny) then
        if (par%output_mode == 1) output%tot(ix,iy)  =output%tot(ix,iy)  +wgt
        if (par%output_mode == 2) output%direc(ix,iy)=output%direc(ix,iy)+wgt
     endif
  endif

  end subroutine add_direct
!--------------------------------------------------
  subroutine peelingoff(photon,grid,observer,output)

  use define
  use raytrace, only : raytrace_to_edge_cyl
  implicit none

  type(photon_type),   intent(in)    :: photon
  type(grid_type),     intent(in)    :: grid
  type(observer_type), intent(in)    :: observer
  type(output_type),   intent(inout) :: output

! local variables
  type (photon_type) :: pobs
  real(kind=wp) :: r2,r,vdet(3)
  real(kind=wp) :: cosa,phot,hgfac,tau
  integer :: ix,iy,i

  pobs    = photon
  pobs%vx = (observer%x-photon%x)
  pobs%vy = (observer%y-photon%y)
  pobs%vz = (observer%z-photon%z)
  r2      = pobs%vx*pobs%vx + pobs%vy*pobs%vy + pobs%vz*pobs%vz
  r       = sqrt(r2)
  pobs%vx = pobs%vx/r
  pobs%vy = pobs%vy/r
  pobs%vz = pobs%vz/r

  call raytrace_to_edge_cyl(pobs,grid,tau)
  cosa  = photon%vx*pobs%vx+photon%vy*pobs%vy+photon%vz*pobs%vz
  hgfac = par%gg1/(par%gg2-par%gg3*cosa)**1.5_wp/fourpi
  phot  = hgfac*photon%lscale/r2 * photon%wgt*exp(-tau)

  do i=1,3
     vdet(i) = observer%rmatrix(i,1) * pobs%vx + &
               observer%rmatrix(i,2) * pobs%vy + &
               observer%rmatrix(i,3) * pobs%vz
  enddo

  !--- Bin the photon into TAN image
  ix = floor(atan2(-vdet(1),vdet(3))*rad2deg/output%dx+output%nx/2.0_wp + par%xshift) + 1
  iy = floor(atan2(-vdet(2),vdet(3))*rad2deg/output%dy+output%ny/2.0_wp + par%yshift) + 1

  if (ix >= 1 .and. ix <= output%nx .and. iy >= 1 .and. iy <= output%ny) then
     !--- place weighted photon into image location
     if (par%output_mode == 1) output%tot(ix,iy)   = output%tot(ix,iy)  +phot
     if (par%output_mode == 2) output%scatt(ix,iy) = output%scatt(ix,iy)+phot
  endif

  !----------
  ! the model is axi-symmetric. (2016-03-19)
  if (par%left_right_fold .and. par%phase_angle == 0.0_wp) then
     pobs%vx = (observer%x+photon%x)
     pobs%vy = (observer%y-photon%y)
     pobs%vz = (observer%z-photon%z)
     !r2      = pobs%vx*pobs%vx + pobs%vy*pobs%vy + pobs%vz*pobs%vz
     !r       = sqrt(r2)
     !pobs%vx = pobs%vx/r
     !pobs%vy = pobs%vy/r
     !pobs%vz = pobs%vz/r

     !cosa  = -photon%vx*pobs%vx+photon%vy*pobs%vy+photon%vz*pobs%vz
     !hgfac = par%gg1/(par%gg2-par%gg3*cosa)**1.5_wp/fourpi
     !phot  = hgfac*photon%lscale/r2 * photon%wgt*exp(-tau)

     do i=1,3
        vdet(i) = observer%rmatrix(i,1) * pobs%vx + &
                  observer%rmatrix(i,2) * pobs%vy + &
                  observer%rmatrix(i,3) * pobs%vz
     enddo

     ix = floor(atan2(-vdet(1),vdet(3))*rad2deg/output%dx+output%nx/2.0_wp + par%xshift) + 1
     iy = floor(atan2(-vdet(2),vdet(3))*rad2deg/output%dy+output%ny/2.0_wp + par%yshift) + 1

     if (ix >= 1 .and. ix <= output%nx .and. iy >= 1 .and. iy <= output%ny) then
        if (par%output_mode == 1) output%tot(ix,iy)   = output%tot(ix,iy)  +phot
        if (par%output_mode == 2) output%scatt(ix,iy) = output%scatt(ix,iy)+phot
     endif
  endif

  end subroutine peelingoff
end module detect
