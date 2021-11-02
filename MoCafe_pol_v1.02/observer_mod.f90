module observer_mod
  use define
  use utility
  use memory_mod
  implicit none
contains

  !-----------------
  subroutine observer_create()
  use mpi
  implicit none
  real(kind=wp) :: cosa,cosb,cosg
  real(kind=wp) :: sina,sinb,sing
  real(kind=wp) :: dist_scale
  integer       :: ierr
  !-- 8 vertices of a cube
  real(kind=wp), parameter :: vertex_x(8) = [1.0,  1.0,  1.0, -1.0, -1.0, -1.0,  1.0,  -1.0]
  real(kind=wp), parameter :: vertex_y(8) = [1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  -1.0]
  real(kind=wp), parameter :: vertex_z(8) = [1.0, -1.0,  1.0,  1.0,  1.0, -1.0, -1.0,  -1.0]
  real(kind=wp) :: px, py, pz, kx, ky, kz, ang_x, ang_y, max_ang_x, max_ang_y
  integer       :: iv

  !--- setup observer and image plane
  ! alpha = phase angle
  ! beta  = inclination angle
  ! gamma = posision angle
  !         Position angle is measured counter-clock wise from x-axis of the detector plane.
  ! Rmaxtrix is the matrix for a rotation starting from the grid frame (initial fixed frame) to observer frame.
  ! rotation sequence is (1) alpha about z-axis,
  !                      (2) beta  about new y-axis, and
  !                      (3) gamma about new z-axis.
  !---
  if (is_finite(par%phase_angle))       par%alpha = -par%phase_angle
  if (is_finite(par%inclination_angle)) par%beta  = -par%inclination_angle
  if (is_finite(par%position_angle))    par%gamma = -par%position_angle

  !--- Default setting: when nothing is given, only a single observer located at a coordinate (0,0,100)*box size.
  if (.not.  (is_finite(par%alpha) .and. is_finite(par%beta) .and. is_finite(par%gamma)) .and. &
      .not.  (is_finite(par%obsx)  .and. is_finite(par%obsy) .and. is_finite(par%obsz))   ) then
     if (.not. is_finite(par%distance)) then
        par%distance = maxval([par%xmax, par%ymax, par%zmax]) * 100.0_wp
     endif
     par%obsx  = 0.0
     par%obsy  = 0.0
     par%obsz  = par%distance

     par%alpha = 0.0
     par%beta  = 0.0
     par%gamma = 0.0
  endif

  if (is_finite(par%obsx) .and. is_finite(par%obsy) .and. is_finite(par%obsz)) then
     if (.not. is_finite(par%distance)) then
        par%distance = sqrt(par%obsx**2 + par%obsy**2 + par%obsz**2)
        if (par%distance < 10.0_wp*maxval([par%xmax, par%ymax, par%zmax])) then
           par%distance = maxval([par%xmax, par%ymax, par%zmax]) * 100.0_wp
        endif
     endif

     if (.not. is_finite(par%gamma)) par%gamma = 0.0_wp

     !--- observer's coordinates represented in the grid system.
     dist_scale = par%distance/sqrt(par%obsx**2 + par%obsy**2 + par%obsz**2)
     observer%x = par%obsx * dist_scale
     observer%y = par%obsy * dist_scale
     observer%z = par%obsz * dist_scale

     cosb = observer%z/par%distance
     sinb = sqrt(1.0_wp - cosb**2)
     if (sinb == 0.0_wp) then
        cosa = 1.0_wp
        sina = 0.0_wp
     else
        cosa = observer%x/(par%distance*sinb)
        sina = observer%y/(par%distance*sinb)
     endif
     cosg = cos(par%gamma*deg2rad)
     sing = sin(par%gamma*deg2rad)
  else
     if (.not. is_finite(par%distance)) then
        par%distance = maxval([par%xmax, par%ymax, par%zmax]) * 100.0_wp
     endif

     !--- Note that the previous convention was to rotate the coordinate system starting from galaxy to observer.
     !--- To following the previous convention, we need to change the sign of the three angles.
     !--- comment added on 2020.08.15
     cosa = cos(par%alpha*deg2rad)
     sina = sin(par%alpha*deg2rad)
     cosb = cos(par%beta *deg2rad)
     sinb = sin(par%beta *deg2rad)
     cosg = cos(par%gamma*deg2rad)
     sing = sin(par%gamma*deg2rad)

     !--- observer's coordinates represented in the grid system.
     observer%x = par%distance*cosa*sinb
     observer%y = par%distance*sina*sinb
     observer%z = par%distance*cosb
  endif

  !--- rotation matrix from Grid system to Observer system.
  observer%rmatrix(1,1) =  cosa*cosb*cosg - sina*sing
  observer%rmatrix(1,2) =  sina*cosb*cosg + cosa*sing
  observer%rmatrix(1,3) = -sinb*cosg
  observer%rmatrix(2,1) = -cosa*cosb*sing - sina*cosg
  observer%rmatrix(2,2) = -sina*cosb*sing + cosa*cosg
  observer%rmatrix(2,3) =  sinb*sing
  observer%rmatrix(3,1) =  cosa*sinb
  observer%rmatrix(3,2) =  sina*sinb
  observer%rmatrix(3,3) =  cosb

  !--- parameters for observer's image plane
  !--- 2020.10.20, determine the image size that convers whole cloud's grid system.
  if (.not. (is_finite(par%dxim) .and. is_finite(par%dyim))) then
     !--- For a spherical geometry, the maximum angular size is always the radius, regardless of the observer direction.
     if (par%rmax > 0.0_wp) then
        par%dxim = atan2(par%rmax,par%distance)/(par%nxim/2.0_wp) * rad2deg
        par%dyim = atan2(par%rmax,par%distance)/(par%nyim/2.0_wp) * rad2deg
     !--- In general, the maxium angular size is determined by 8 vertices.
     else
        max_ang_x = -999.0
        max_ang_y = -999.0
        do iv = 1,8
           px = observer%x - vertex_x(iv) * par%xmax
           py = observer%y - vertex_y(iv) * par%ymax
           pz = observer%z - vertex_z(iv) * par%zmax
           kx = observer%rmatrix(1,1)*px + observer%rmatrix(1,2)*py + observer%rmatrix(1,3)*pz
           ky = observer%rmatrix(2,1)*px + observer%rmatrix(2,2)*py + observer%rmatrix(2,3)*pz
           kz = observer%rmatrix(3,1)*px + observer%rmatrix(3,2)*py + observer%rmatrix(3,3)*pz
           ang_x = abs(atan2(-kx, kz))
           ang_y = abs(atan2(-ky, kz))
           if (max_ang_x < ang_x) max_ang_x = ang_x
           if (max_ang_y < ang_y) max_ang_y = ang_y
        enddo
        if (par%nxim == par%nyim) then
           par%dxim = maxval([max_ang_x, max_ang_y]) / (par%nxim/2.0_wp) * rad2deg
           par%dyim = maxval([max_ang_x, max_ang_y]) / (par%nyim/2.0_wp) * rad2deg
        else
           if (.not. is_finite(par%dxim)) par%dxim = max_ang_x / (par%nxim/2.0_wp) * rad2deg
           if (.not. is_finite(par%dyim)) par%dyim = max_ang_y / (par%nyim/2.0_wp) * rad2deg
        endif
     endif
     !if (.not. is_finite(par%dxim)) par%dxim = atan2(par%xmax,par%distance)/(par%nxim/2.0_wp) * rad2deg
     !if (.not. is_finite(par%dyim)) par%dyim = atan2(par%ymax,par%distance)/(par%nyim/2.0_wp) * rad2deg
  endif
  observer%dxim = par%dxim
  observer%dyim = par%dyim
  observer%nxim = par%nxim
  observer%nyim = par%nyim

  observer%steradian_pix = observer%dxim * observer%dyim * (deg2rad**2)

  ! memory allocations of output images
  call create_mem(observer%scatt, [observer%nxim,observer%nyim])
  call create_mem(observer%direc, [observer%nxim,observer%nyim])
  if (par%use_stokes) then
     call create_mem(observer%I,  [observer%nxim,observer%nyim])
     call create_mem(observer%Q,  [observer%nxim,observer%nyim])
     call create_mem(observer%U,  [observer%nxim,observer%nyim])
     call create_mem(observer%V,  [observer%nxim,observer%nyim])
  endif
  if (par%sightline_tau) then
     call create_mem(observer%tau,[observer%nxim,observer%nyim])
  endif
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  !--- let's do not consider xyz_symmetry condition. (2017-08-19)
  !--- later, the xyz_symmetry condition will be implemented.
  !--- for a goemetry with xyz_symmetry, the peeling-off technique can be used only for special cases.
  par%xyz_symmetry = .false.
  end subroutine observer_create
  !-----------------
  subroutine observer_destroy()
  use define
  use mpi
  implicit none
  integer :: ierr

  call destroy_shared_mem_all()
  if (associated(observer%scatt)) deallocate(observer%scatt)
  if (associated(observer%direc)) deallocate(observer%direc)

  if (par%use_stokes) then
     if (associated(observer%I)) deallocate(observer%I)
     if (associated(observer%Q)) deallocate(observer%Q)
     if (associated(observer%U)) deallocate(observer%U)
     if (associated(observer%V)) deallocate(observer%V)
  endif
  end subroutine observer_destroy
  !---------------------------------------------------------------------------------
end module observer_mod
