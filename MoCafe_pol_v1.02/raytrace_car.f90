module raytrace
contains

  subroutine raytrace_to_edge_car(photon0,grid,tau)
!--- Find the cumulative optical depth to the edge of grid.
!--- Note that photon0.(x0,y0,z0) and photon0.(icell0,jcell0,kcell0) do not change.
!---
!--- Basic algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012
!--- Fixed minor bugs at grid boundaries: 2013-05-24
!--- Bug fixed at path: 2015-02-25

  use define
  implicit none

  type(photon_type), intent(in) :: photon0
  type(grid_type),   intent(in) :: grid
  real(kind=wp),    intent(out) :: tau

! local variables
  !-- note that exp(-tau_huge) = 0.0
  !-- at this large optical depth peeling-off contribution is practically zero.
  !-- Hence, no need to integrate further. This is not required for dust simulation (2020.11.08).
  !real(kind=wp), parameter :: tau_huge = 745.2_wp
  integer       :: icell,jcell,kcell
  integer       :: istep,jstep,kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: d
  real(kind=wp) :: tx,ty,tz,delx,dely,delz
  real(kind=wp) :: rhokap
  integer       :: idx_min

!--- (xp,yp,zp) and (icell,jcell,kcell) = the photon coordinates and cell indices
  xp = photon0%x
  yp = photon0%y
  zp = photon0%z
  kx = photon0%kx
  ky = photon0%ky
  kz = photon0%kz
  icell = photon0%icell
  jcell = photon0%jcell
  kcell = photon0%kcell

!--- tau = cumulative optical depth
!--- d   = cumulative path length
  tau        = 0.0_wp
  d          = 0.0_wp

  if (kx > 0.0_wp) then
     ! icell = nx+1 & xface(nx+1) = xp & kx > 0.0
     ! Ignore the case of going out of the grid system.
     if (icell > grid%nx .and. grid%xface(icell) <= xp) return
     istep = 1
     tx    = (grid%xface(icell+1)-xp)/kx
     delx  = grid%dx/kx
  else if (kx < 0.0_wp) then
     if (grid%xface(icell) == xp) then
        if (icell > 1) then
           ! icell = 2, 3, or, ... & xface(icell) = xp & kx < 0.0
           ! Skip the case of tau = 0
           icell = icell - 1
        else
           ! icell = 1 & xface(1) = xp & kx < 0.0
           ! Ignore the case of going out of the grid system.
           return
        endif
     endif
     istep = -1
     tx    = (grid%xface(icell)-xp)/kx
     delx  = -grid%dx/kx
  else
     istep = 0
     tx    = hugest
     delx  = hugest
  endif

  if (ky > 0.0_wp) then
     if (jcell > grid%ny .and. grid%yface(jcell) <= yp) return
     jstep = 1
     ty    = (grid%yface(jcell+1)-yp)/ky
     dely  = grid%dy/ky
  else if (ky < 0.0_wp) then
     if (grid%yface(jcell) == yp) then
        if (jcell > 1) then
           jcell = jcell - 1
        else
           return
        endif
     endif
     jstep = -1
     ty    = (grid%yface(jcell)-yp)/ky
     dely  = -grid%dy/ky
  else
     jstep = 0
     ty    = hugest
     dely  = hugest
  endif

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) <= zp) return
     kstep = 1
     tz    = (grid%zface(kcell+1)-zp)/kz
     delz  = grid%dz/kz
  else if (kz < 0.0_wp) then
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           return
        endif
     endif
     kstep = -1
     tz    = (grid%zface(kcell)-zp)/kz
     delz  = -grid%dz/kz
  else
     kstep = 0
     tz    = hugest
     delz  = hugest
  endif

!--- integrate through grid
  do while(.true.)
     rhokap = grid%rhokap(icell,jcell,kcell)

     idx_min = minloc([tx,ty,tz], dim=1)

     if (idx_min == 1) then
        tau   = tau + (tx - d) * rhokap
        d     = tx
        icell = icell + istep
        if (icell < 1 .or. icell > grid%nx) exit
        tx    = tx + delx
     else if (idx_min == 2) then
        tau   = tau + (ty - d) * rhokap
        d     = ty
        jcell = jcell + jstep
        if (jcell < 1 .or. jcell > grid%ny) exit
        ty = ty + dely
     else
        tau   = tau + (tz - d) * rhokap
        d     = tz
        kcell = kcell + kstep
        if (kcell < 1 .or. kcell > grid%nz) exit
        tz    = tz + delz
     endif
     !--- for dust simulation, this limitation is not necessary (2020.11.08).
     !if (tau >= tau_huge) exit
  enddo
  return
  end subroutine raytrace_to_edge_car

  subroutine raytrace_to_tau_car(photon,grid,tau_in)
!--- Find the coordinates and cell indices corresponding to an input optical depth tau_in.
!---
!--- The algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012-09-24
!--- Fixed minor bugs at grid boundaries: 2013-05-24

  use define
  implicit none

  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  real(kind=wp),     intent(in)    :: tau_in

! Local variables
  integer :: icell,jcell,kcell
  integer :: iold, jold, kold
  integer :: istep,jstep,kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: tau,d_overshoot,d
  real(kind=wp) :: tx,ty,tz,delx,dely,delz
  real(kind=wp) :: rhokap
  real(kind=wp) :: dold
  integer       :: idx_min

!--- xp, yp, zp = the current photon coordinates.
  xp = photon%x
  yp = photon%y
  zp = photon%z
  kx = photon%kx
  ky = photon%ky
  kz = photon%kz
  icell = photon%icell
  jcell = photon%jcell
  kcell = photon%kcell

  if (kx > 0.0_wp) then
     if (icell > grid%nx .and. grid%xface(icell) == xp) then
        ! icell = nx+1 & xface(nx+1) = xp & kx > 0.0
        ! Ignore the case of going out of the grid system.
        photon%inside = .false.
        return
     endif
     istep = 1
     tx   = (grid%xface(icell+1)-xp)/kx
     delx = grid%dx/kx
  else if (kx < 0.0_wp) then
     if (grid%xface(icell) == xp) then
        if (icell > 1) then
           ! icell = 2, 3, or, ... & xface(icell) = xp & kx < 0.0
           ! Skip the case of tau = 0
           icell = icell - 1
        else
           ! icell = 1 & xface(1) = xp & kx < 0.0
           ! Ignore the case of going out of the grid system.
           photon%inside = .false.
           return
        endif
     endif
     istep = -1
     tx   = (grid%xface(icell)-xp)/kx
     delx = -grid%dx/kx
  else
     istep = 0
     tx   = hugest
     delx = hugest
  endif

  if (ky > 0.0_wp) then
     if (jcell > grid%ny .and. grid%yface(jcell) == yp) then
        photon%inside = .false.
        return
     endif
     jstep = 1
     ty   = (grid%yface(jcell+1)-yp)/ky
     dely = grid%dy/ky
  else if (ky < 0.0_wp) then
     if (grid%yface(jcell) == yp) then
        if (jcell > 1) then
           jcell = jcell - 1
        else
           photon%inside = .false.
           return
        endif
     endif
     jstep = -1
     ty   = (grid%yface(jcell)-yp)/ky
     dely = -grid%dy/ky
  else
     jstep = 0
     ty   = hugest
     dely = hugest
  endif

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) == zp) then
        photon%inside = .false.
        return
     endif
     kstep = 1
     tz   = (grid%zface(kcell+1)-zp)/kz
     delz = grid%dz/kz
  else if (kz < 0.0_wp) then
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           photon%inside = .false.
           return
        endif
     endif
     kstep = -1
     tz   = (grid%zface(kcell)-zp)/kz
     delz = -grid%dz/kz
  else
     kstep = 0
     tz   = hugest
     delz = hugest
  endif

  !--- d   = cumulative path length
  !--- tau = cumulative optical depth
  tau    = 0.0_wp
  d      = 0.0_wp
  dold   = 0.0_wp
  iold   = icell
  jold   = jcell
  kold   = kcell

!--- integrate through grid
  do while(photon%inside)

     rhokap = grid%rhokap(icell,jcell,kcell)

     idx_min = minloc([tx,ty,tz], dim=1)

     if (idx_min == 1) then
        tau    = tau + (tx - d) * rhokap
        d      = tx
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d  = d - d_overshoot
           endif
           xp = xp + d * kx
           yp = yp + d * ky
           zp = zp + d * kz
           exit
        endif
        icell = icell + istep
        if (icell < 1 .or. icell > grid%nx) then
           photon%inside = .false.
           exit
        endif
        tx = tx + delx
     else if (idx_min == 2) then
        tau    = tau + (ty - d) * rhokap
        d      = ty
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d  = d - d_overshoot
           endif
           xp = xp + d * kx
           yp = yp + d * ky
           zp = zp + d * kz
           exit
        endif
        jcell = jcell + jstep
        if (jcell < 1 .or. jcell > grid%ny) then
           photon%inside = .false.
           exit
        endif
        ty = ty + dely
     else
        tau    = tau + (tz - d) * rhokap
        d      = tz
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d  = d - d_overshoot
           endif
           xp = xp + d * kx
           yp = yp + d * ky
           zp = zp + d * kz
           exit
        endif
        kcell = kcell + kstep
        if (kcell < 1 .or. kcell > grid%nz) then
           photon%inside = .false.
           exit
        endif
        tz = tz + delz
     endif

     iold = icell
     jold = jcell
     kold = kcell
     dold = d
  enddo

  if (.not. photon%inside) then
     icell = iold
     jcell = jold
     kcell = kold
  endif

!--- xp, yp, zp = the current photon coordinates.
  photon%x = xp
  photon%y = yp
  photon%z = zp
  photon%icell = icell
  photon%jcell = jcell
  photon%kcell = kcell

  return
  end subroutine raytrace_to_tau_car

  subroutine raytrace_to_tau_car_xyzsym(photon,grid,tau_in)
!--- Find the coordinates and cell indices corresponding to an input optical depth tau_in.
!---
!--- The algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012-09-24
!--- Fixed minor bugs at grid boundaries: 2013-05-24
!--- 2017-05-15, Reflected at xy-, xz- and yz-planes.

  use define
  implicit none

  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  real(kind=wp),     intent(in)    :: tau_in

! Local variables
  integer :: icell,jcell,kcell
  integer :: iold, jold, kold
  integer :: istep,jstep,kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: tau,d_overshoot,d
  real(kind=wp) :: tx,ty,tz,delx,dely,delz
  real(kind=wp) :: rhokap
  real(kind=wp) :: dold
  integer       :: idx_min

!--- xp, yp, zp = the current photon coordinates.
  xp = photon%x
  yp = photon%y
  zp = photon%z
  kx = photon%kx
  ky = photon%ky
  kz = photon%kz
  icell = photon%icell
  jcell = photon%jcell
  kcell = photon%kcell

  if (kx > 0.0_wp) then
     if (icell > grid%nx .and. grid%xface(icell) == xp) then
        ! icell = nx+1 & xface(nx+1) = xp & kx > 0.0
        ! Ignore the case of going out of the grid system.
        photon%inside = .false.
        return
     endif
     istep = 1
     tx    = (grid%xface(icell+1)-xp)/kx
     delx  = grid%dx/kx
  else if (kx < 0.0_wp) then
     istep = -1
     delx  = -grid%dx/kx
     if (grid%xface(icell) == xp) then
        if (icell > 1) then
           ! icell = 2, 3, or, ... & xface(icell) = xp & kx < 0.0
           ! Skip the case of tau = 0
           icell = icell - 1
        else
           ! icell = 1 & xface(1) = xp & kx < 0.0
           icell = grid%i0
           istep = 1
           kx    = abs(kx)
        endif
     endif
     if (istep == 1) then
        tx = (grid%xface(icell+1)-xp)/kx
     else
        tx = (grid%xface(icell)-xp)/kx
     endif
  else
     istep = 0
     tx   = hugest
     delx = hugest
  endif

  if (ky > 0.0_wp) then
     if (jcell > grid%ny .and. grid%yface(jcell) == yp) then
        photon%inside = .false.
        return
     endif
     jstep = 1
     ty    = (grid%yface(jcell+1)-yp)/ky
     dely  = grid%dy/ky
  else if (ky < 0.0_wp) then
     jstep = -1
     dely  = -grid%dy/ky
     if (grid%yface(jcell) == yp) then
        if (jcell > 1) then
           jcell = jcell - 1
        else
           jcell = grid%j0
           jstep = 1
           ky    = abs(ky)
        endif
     endif
     if (jstep == 1) then
        ty = (grid%yface(jcell+1)-yp)/ky
     else
        ty = (grid%yface(jcell)-yp)/ky
     endif
  else
     jstep = 0
     ty   = hugest
     dely = hugest
  endif

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) == zp) then
        photon%inside = .false.
        return
     endif
     kstep = 1
     tz   = (grid%zface(kcell+1)-zp)/kz
     delz = grid%dz/kz
  else if (kz < 0.0_wp) then
     kstep = -1
     delz  = -grid%dz/kz
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           kcell = grid%k0
           kstep = 1
           kz    = abs(kz)
        endif
     endif
     if (kstep == 1) then
        tz = (grid%zface(kcell+1)-zp)/kz
     else
        tz = (grid%zface(kcell)-zp)/kz
     endif
  else
     kstep = 0
     tz   = hugest
     delz = hugest
  endif

  !--- d   = cumulative path length
  !--- tau = cumulative optical depth
  tau    = 0.0_wp
  d      = 0.0_wp
  dold   = 0.0_wp
  iold   = icell
  jold   = jcell
  kold   = kcell

!--- integrate through grid
  do while(photon%inside)

     rhokap = grid%rhokap(icell,jcell,kcell)

     idx_min = minloc([tx,ty,tz], dim=1)

     if (idx_min == 1) then
        tau    = tau + (tx - d) * rhokap
        d      = tx
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d  = d - d_overshoot
           endif
           exit
        endif
        icell = icell + istep
        if (icell < 1) then
           icell = grid%i0
           istep = 1
           kx    = -kx
        else if (icell > grid%nx) then
           photon%inside = .false.
           exit
        endif
        tx     = tx + delx
     else if (idx_min == 2) then
        tau    = tau + (ty - d) * rhokap
        d      = ty
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d  = d - d_overshoot
           endif
           exit
        endif
        jcell = jcell + jstep
        if (jcell < 1) then
           jcell = grid%j0
           jstep = 1
           ky    = -ky
        else if (jcell > grid%ny) then
           photon%inside = .false.
           exit
        endif
        ty     = ty + dely
     else
        tau    = tau + (tz - d) * rhokap
        d      = tz
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d  = d - d_overshoot
           endif
           exit
        endif
        kcell = kcell + kstep
        if (kcell < 1) then
           !kcell = 1
           kcell = grid%k0
           kstep = 1
           kz    = -kz
        else if (kcell > grid%nz) then
           photon%inside = .false.
           exit
        endif
        tz    = tz + delz
     endif

     iold   = icell
     jold   = jcell
     kold   = kcell
     dold   = d
  enddo

  if (.not. photon%inside) then
     icell = iold
     jcell = jold
     kcell = kold
  endif

!--- xp, yp, zp = the current photon coordinates.
  ! here, photon%kx, etc should be used instead of kx,ky,kz
  photon%x  = photon%x + d*photon%kx
  photon%y  = photon%y + d*photon%ky
  photon%z  = photon%z + d*photon%kz
  if (photon%x < grid%xmin) photon%x = -photon%x
  if (photon%y < grid%ymin) photon%y = -photon%y
  if (photon%z < grid%zmin) photon%z = -photon%z
  photon%kx = kx
  photon%ky = ky
  photon%kz = kz
  photon%icell = icell
  photon%jcell = jcell
  photon%kcell = kcell

  return
  end subroutine raytrace_to_tau_car_xyzsym

  subroutine raytrace_to_tau_car_xyper(photon,grid,tau_in)
!--- Find the coordinates and cell indices corresponding to an input optical depth tau_in.
!---
!--- The algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012-09-24
!--- Fixed minor bugs at grid boundaries: 2013-05-24
!--- Period Boundary Codition around x, y axis for plane-paralle geometry (2016-12-07)

  use define
  implicit none

  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  real(kind=wp),     intent(in)    :: tau_in

! Local variables
  integer :: icell,jcell,kcell
  integer :: iold, jold, kold
  integer :: istep,jstep,kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: tau,d_overshoot,d
  real(kind=wp) :: tx,ty,tz,delx,dely,delz
  real(kind=wp) :: rhokap
  real(kind=wp) :: dold
  integer       :: idx_min

!--- xp, yp, zp = the current photon coordinates.
  xp = photon%x
  yp = photon%y
  zp = photon%z
  kx = photon%kx
  ky = photon%ky
  kz = photon%kz
  icell = photon%icell
  jcell = photon%jcell
  kcell = photon%kcell

  if (kx > 0.0_wp) then
     if (icell > grid%nx .and. grid%xface(icell) == xp) then
        ! icell = nx+1 & xface(nx+1) = xp & kx > 0.0
        icell = 1
        xp    = grid%xface(icell)
     endif
     istep = 1
     tx   = (grid%xface(icell+1)-xp)/kx
     delx = grid%dx/kx
  else if (kx < 0.0_wp) then
     if (grid%xface(icell) == xp) then
        if (icell > 1) then
           ! icell = 2, 3, or, ... & xface(icell) = xp & kx < 0.0
           ! Skip the case of tau = 0
           icell = icell - 1
        else
           ! icell = 1 & xface(1) = xp & kx < 0.0
           icell = grid%nx
           xp    = grid%xface(icell+1)
        endif
     endif
     istep = -1
     tx   = (grid%xface(icell)-xp)/kx
     delx = -grid%dx/kx
  else
     istep = 0
     tx   = hugest
     delx = hugest
  endif

  if (ky > 0.0_wp) then
     if (jcell > grid%ny .and. grid%yface(jcell) == yp) then
        jcell = 1
        yp    = grid%yface(jcell)
     endif
     jstep = 1
     ty   = (grid%yface(jcell+1)-yp)/ky
     dely = grid%dy/ky
  else if (ky < 0.0_wp) then
     if (grid%yface(jcell) == yp) then
        if (jcell > 1) then
           jcell = jcell - 1
        else
           jcell = grid%ny
           yp    = grid%yface(jcell+1)
        endif
     endif
     jstep = -1
     ty   = (grid%yface(jcell)-yp)/ky
     dely = -grid%dy/ky
  else
     jstep = 0
     ty   = hugest
     dely = hugest
  endif

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) == zp) then
        photon%inside = .false.
        return
     endif
     kstep = 1
     tz   = (grid%zface(kcell+1)-zp)/kz
     delz = grid%dz/kz
  else if (kz < 0.0_wp) then
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           photon%inside = .false.
           return
        endif
     endif
     kstep = -1
     tz   = (grid%zface(kcell)-zp)/kz
     delz = -grid%dz/kz
  else
     kstep = 0
     tz   = hugest
     delz = hugest
  endif

  !--- d   = cumulative path length
  !--- tau = cumulative optical depth
  d     = 0.0_wp
  tau   = 0.0_wp
  dold  = 0.0_wp
  iold  = icell
  jold  = jcell
  kold  = kcell

  !--- integrate through grid
  do while(photon%inside)

     rhokap = grid%rhokap(icell,jcell,kcell)

     idx_min = minloc([tx, ty, tz], dim=1)
     if (idx_min == 1) then
        tau    = tau + (tx - d) * rhokap
        d      = tx
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d  = d - d_overshoot
           endif
           xp = xp + d * kx
           yp = yp + d * ky
           zp = zp + d * kz
           exit
        endif
        icell = icell + istep
        if (icell < 1)       icell = grid%nx
        if (icell > grid%nx) icell = 1
        tx = tx + delx
     else if (idx_min == 2) then
        tau    = tau + (ty - d) * rhokap
        d      = ty
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d  = d - d_overshoot
           endif
           xp = xp + d * kx
           yp = yp + d * ky
           zp = zp + d * kz
           exit
        endif
        jcell = jcell + jstep
        if (jcell < 1)       jcell = grid%ny
        if (jcell > grid%ny) jcell = 1
        ty = ty + dely
     else
        tau    = tau + (tz - d) * rhokap
        d      = tz
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d  = d - d_overshoot
           endif
           xp = xp + d * kx
           yp = yp + d * ky
           zp = zp + d * kz
           exit
        endif
        kcell = kcell + kstep
        if (kcell < 1 .or. kcell > grid%nz) then
           photon%inside = .false.
           exit
        endif
        tz = tz + delz
     endif

     iold = icell
     jold = jcell
     kold = kcell
     dold = d
  enddo

  if (.not. photon%inside) then
     icell = iold
     jcell = jold
     kcell = kold
  endif

!--- xp, yp, zp = the current photon coordinates.
  ! fold back into the system (2016-12-11)
  photon%x  = xp - floor((xp-grid%xmin)/grid%xrange)*grid%xrange
  photon%y  = yp - floor((yp-grid%ymin)/grid%yrange)*grid%yrange
  photon%z  = zp
  photon%icell = icell
  photon%jcell = jcell
  photon%kcell = kcell

  return
  end subroutine raytrace_to_tau_car_xyper

  subroutine raytrace_to_tau_car_zonly(photon,grid,tau_in)
!--- Find the coordinates and cell indices corresponding to an input optical depth tau_in.
!---
!--- The algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012-09-24
!--- Fixed minor bugs at grid boundaries: 2013-05-24

  use define
  implicit none

  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  real(kind=wp),     intent(in)    :: tau_in

! Local variables
  integer :: icell,jcell,kcell
  integer :: iold, jold, kold
  integer :: kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: tau,d_overshoot,d
  real(kind=wp) :: tz,delz
  real(kind=wp) :: rhokap
  real(kind=wp) :: dold

!--- xp, yp, zp = the current photon coordinates.
  xp = photon%x
  yp = photon%y
  zp = photon%z
  kx = photon%kx
  ky = photon%ky
  kz = photon%kz
  icell = photon%icell
  jcell = photon%jcell
  kcell = photon%kcell

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) == zp) then
        photon%inside = .false.
        return
     endif
     kstep = 1
     tz   = (grid%zface(kcell+1)-zp)/kz
     delz = grid%dz/kz
  else if (kz < 0.0_wp) then
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           photon%inside = .false.
           return
        endif
     endif
     kstep = -1
     tz   = (grid%zface(kcell)-zp)/kz
     delz = -grid%dz/kz
  else
     kstep = 0
     tz   = hugest
     delz = hugest
  endif

  !--- d   = cumulative path length
  !--- tau = cumulative optical depth
  tau    = 0.0_wp
  d      = 0.0_wp
  dold   = 0.0_wp
  iold   = icell
  jold   = jcell
  kold   = kcell

!--- integrate through grid
  do while(photon%inside)

     rhokap = grid%rhokap(icell,jcell,kcell)

     tau    = tau + (tz - d) * rhokap
     d      = tz
     if (tau >= tau_in) then
        if (rhokap > 0.0) then
           d_overshoot = (tau - tau_in)/rhokap
           d  = d - d_overshoot
        endif
        xp = xp + d * kx
        yp = yp + d * ky
        zp = zp + d * kz
        exit
     endif
     kcell = kcell + kstep
     if (kcell < 1 .or. kcell > grid%nz) then
        photon%inside = .false.
        exit
     endif
     tz = tz + delz

     kold = kcell
     dold = d
  enddo

  if (.not. photon%inside) then
     kcell = kold
  endif

!--- xp, yp, zp = the current photon coordinates.
  photon%x = xp
  photon%y = yp
  photon%z = zp
  photon%kcell = kcell

  return
  end subroutine raytrace_to_tau_car_zonly

end module raytrace
