module raytrace
contains

  subroutine raytrace_to_edge_car(photon0,grid,tau)
!--- Find the cumulative optical depth to the edge of grid.
!--- Note that photon0.(x0,y0,z0) and photon0.(icell0,jcell0,kcell0) do not change.
!---
!--- The algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012
!--- Fixed minor bugs at grid boundaries: 2013-05-24

  use define
  implicit none

  type(photon_type), intent(in)  :: photon0
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(out) :: tau

! Local variables
  integer       :: icell,jcell,kcell
  integer       :: istep,jstep,kstep
  real(kind=wp) :: xp,yp,zp,vx,vy,vz
  real(kind=wp) :: d,tx,ty,tz,delx,dely,delz

!--- (xp,yp,zp) and (icell,jcell,kcell) = the photon coordinates and cell indices
  xp = photon0%x
  yp = photon0%y
  zp = photon0%z
  vx = photon0%vx
  vy = photon0%vy
  vz = photon0%vz
  icell = photon0%icell
  jcell = photon0%jcell
  kcell = photon0%kcell

!--- tau = cumulative optical depth
!--- d   = cumulative path length
  tau = 0.0_wp
  d   = 0.0_wp

  if (vx > 0.0_wp) then
     ! icell = nx+1 & xface(nx+1) = xp & vx > 0.0
     ! Ignore the case of going out of the grid system.
     if (icell > grid%nx .and. grid%xface(icell) == xp) return
     istep = 1
     tx   = (grid%xface(icell+1)-xp)/vx
     delx = grid%dx/vx
  else if (vx < 0.0_wp) then
     if (grid%xface(icell) == xp) then
        if (icell > 1) then
           ! icell = 2, 3, or, ... & xface(icell) = xp & vx < 0.0
           ! Skip the case of tau = 0
           icell = icell - 1
        else
           ! icell = 1 & xface(1) = xp & vx < 0.0
           ! Ignore the case of going out of the grid system.
           return
        endif
     endif
     istep = -1
     tx   = (grid%xface(icell)-xp)/vx
     delx = -grid%dx/vx
  else
     istep = 0
     tx   = hugest
     delx = hugest
  endif

  if (vy > 0.0_wp) then
     if (jcell > grid%ny .and. grid%yface(jcell) == yp) return
     jstep = 1
     ty   = (grid%yface(jcell+1)-yp)/vy
     dely = grid%dy/vy
  else if (vy < 0.0_wp) then
     if (grid%yface(jcell) == yp) then
        if (jcell > 1) then
           jcell = jcell - 1
        else
           return
        endif
     endif
     jstep = -1
     ty   = (grid%yface(jcell)-yp)/vy
     dely = -grid%dy/vy
  else
     jstep = 0
     ty   = hugest
     dely = hugest
  endif

  if (vz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) == zp) return
     kstep = 1
     tz   = (grid%zface(kcell+1)-zp)/vz
     delz = grid%dz/vz
  else if (vz < 0.0_wp) then
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           return
        endif
     endif
     kstep = -1
     tz   = (grid%zface(kcell)-zp)/vz
     delz = -grid%dz/vz
  else
     kstep = 0
     tz   = hugest
     delz = hugest
  endif

!--- integrate through grid
  do while(.true.)
     if (tx < ty) then
        if (tx < tz) then
           tau   = tau + (tx - d) * grid%opacity(icell,jcell,kcell)
           d     = tx
           icell = icell + istep
           if (icell < 1 .or. icell > grid%nx) exit
           tx = tx + delx
        else
           tau   = tau + (tz - d) * grid%opacity(icell,jcell,kcell)
           d     = tz
           kcell = kcell + kstep
           if (kcell < 1 .or. kcell > grid%nz) exit
           tz = tz + delz
           if (tx == tz) then
              icell = icell + istep
              if (icell < 1 .or. icell > grid%nx) exit
              tx = tx + delx
           endif
        endif
     else
        if (ty < tz) then
           tau   = tau + (ty - d) * grid%opacity(icell,jcell,kcell)
           d     = ty
           jcell = jcell + jstep
           if (jcell < 1 .or. jcell > grid%ny) exit
           ty = ty + dely
        else
           tau   = tau + (tz - d) * grid%opacity(icell,jcell,kcell)
           d     = tz
           kcell = kcell + kstep
           if (kcell < 1 .or. kcell > grid%nz) exit
           tz = tz + delz
           if (ty == tz) then
              jcell = jcell + jstep
              if (jcell < 1 .or. jcell > grid%ny) exit
              ty = ty + dely
           endif
        endif
        if (tx == ty) then
           icell = icell + istep
           if (icell < 1 .or. icell > grid%nx) exit
           tx = tx + delx
        endif
     endif
  enddo

  return
  end subroutine raytrace_to_edge_car

  subroutine raytrace_to_tau_car(photon,grid,tau_in,inside_grid)
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
  type(grid_type),   intent(in)    :: grid
  real(kind=wp),     intent(in)    :: tau_in
  logical,           intent(inout) :: inside_grid


! Local variables
  integer :: icell,jcell,kcell
  integer :: istep,jstep,kstep
  real(kind=wp) :: xp,yp,zp,vx,vy,vz
  real(kind=wp) :: tau,d_overshoot,d
  real(kind=wp) :: tx,ty,tz,delx,dely,delz

!--- xp, yp, zp = the current photon coordinates.
  xp = photon%x
  yp = photon%y
  zp = photon%z
  vx = photon%vx
  vy = photon%vy
  vz = photon%vz
  icell = photon%icell
  jcell = photon%jcell
  kcell = photon%kcell

!--- d   = cumulative path length
!--- tau = cumulative optical depth
  d   = 0.0_wp
  tau = 0.0_wp

  if (vx > 0.0_wp) then
     if (icell > grid%nx .and. grid%xface(icell) == xp) then
        ! icell = nx+1 & xface(nx+1) = xp & vx > 0.0
        ! Ignore the case of going out of the grid system.
        inside_grid = .false.
        return
     endif
     istep = 1
     tx   = (grid%xface(icell+1)-xp)/vx
     delx = grid%dx/vx
  else if (vx < 0.0_wp) then
     if (grid%xface(icell) == xp) then
        if (icell > 1) then
           ! icell = 2, 3, or, ... & xface(icell) = xp & vx < 0.0
           ! Skip the case of tau = 0
           icell = icell - 1
        else
           ! icell = 1 & xface(1) = xp & vx < 0.0
           ! Ignore the case of going out of the grid system.
           inside_grid = .false.
           return
        endif
     endif
     istep = -1
     tx   = (grid%xface(icell)-xp)/vx
     delx = -grid%dx/vx
  else
     istep = 0
     tx   = hugest
     delx = hugest
  endif

  if (vy > 0.0_wp) then
     if (jcell > grid%ny .and. grid%yface(jcell) == yp) then
        inside_grid = .false.
        return
     endif
     jstep = 1
     ty   = (grid%yface(jcell+1)-yp)/vy
     dely = grid%dy/vy
  else if (vy < 0.0_wp) then
     if (grid%yface(jcell) == yp) then
        if (jcell > 1) then
           jcell = jcell - 1
        else
           inside_grid = .false.
           return
        endif
     endif
     jstep = -1
     ty   = (grid%yface(jcell)-yp)/vy
     dely = -grid%dy/vy
  else
     jstep = 0
     ty   = hugest
     dely = hugest
  endif

  if (vz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) == zp) then
        inside_grid = .false.
        return
     endif
     kstep = 1
     tz   = (grid%zface(kcell+1)-zp)/vz
     delz = grid%dz/vz
  else if (vz < 0.0_wp) then
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           inside_grid = .false.
           return
        endif
     endif
     kstep = -1
     tz   = (grid%zface(kcell)-zp)/vz
     delz = -grid%dz/vz
  else
     kstep = 0
     tz   = hugest
     delz = hugest
  endif

!--- integrate through grid
  do while(inside_grid)
     if (tx < ty) then
        if (tx < tz) then
           tau = tau + (tx - d) * grid%opacity(icell,jcell,kcell)
           d   = tx
           if (tau >= tau_in) then
              if (grid%opacity(icell,jcell,kcell) > 0.0) then
                 d_overshoot = (tau - tau_in)/grid%opacity(icell,jcell,kcell)
                 d  = d - d_overshoot
              endif
              xp = xp + d * vx
              yp = yp + d * vy
              zp = zp + d * vz
              exit
           endif
           icell = icell + istep
           if (icell < 1 .or. icell > grid%nx) then
              inside_grid = .false.
              exit
           endif
           tx = tx + delx
        else
           tau = tau + (tz - d) * grid%opacity(icell,jcell,kcell)
           d   = tz
           if (tau >= tau_in) then
              if (grid%opacity(icell,jcell,kcell) > 0.0) then
                  d_overshoot = (tau - tau_in)/grid%opacity(icell,jcell,kcell)
                  d  = d - d_overshoot
              endif
              xp = xp + d * vx
              yp = yp + d * vy
              zp = zp + d * vz
              exit
           endif
           kcell = kcell + kstep
           if (kcell < 1 .or. kcell > grid%nz) then
              inside_grid = .false.
              exit
           endif
           tz = tz + delz
           if (tx == tz) then
              icell = icell + istep
              if (icell < 1 .or. icell > grid%nx) then
                 inside_grid = .false.
                 exit
              endif
              tx = tx + delx
           endif
        endif
     else
        if (ty < tz) then
           tau = tau + (ty - d) * grid%opacity(icell,jcell,kcell)
           d   = ty
           if (tau >= tau_in) then
              if (grid%opacity(icell,jcell,kcell) > 0.0) then
                 d_overshoot = (tau - tau_in)/grid%opacity(icell,jcell,kcell)
                 d  = d - d_overshoot
              endif
              xp = xp + d * vx
              yp = yp + d * vy
              zp = zp + d * vz
              exit
           endif
           jcell = jcell + jstep
           if (jcell < 1 .or. jcell > grid%ny) then
              inside_grid = .false.
              exit
           endif
           ty = ty + dely
        else
           tau = tau + (tz - d) * grid%opacity(icell,jcell,kcell)
           d   = tz
           if (tau >= tau_in) then
              if (grid%opacity(icell,jcell,kcell) > 0.0) then
                 d_overshoot = (tau - tau_in)/grid%opacity(icell,jcell,kcell)
                 d  = d - d_overshoot
              endif
              xp = xp + d * vx
              yp = yp + d * vy
              zp = zp + d * vz
              exit
           endif
           kcell = kcell + kstep
           if (kcell < 1 .or. kcell > grid%nz) then
              inside_grid = .false.
              exit
           endif
           tz = tz + delz
           if (ty == tz) then
              jcell = jcell + jstep
              if (jcell < 1 .or. jcell > grid%ny) then
                 inside_grid = .false.
                 exit
              endif
              ty = ty + dely
           endif
        endif
        if (tx == ty) then
           icell = icell + istep
           if (icell < 1 .or. icell > grid%nx) then
              inside_grid = .false.
              exit
           endif
           tx = tx + delx
        endif
     endif
  enddo

!--- xp, yp, zp = the current photon coordinates.
  photon%x  = xp
  photon%y  = yp
  photon%z  = zp
  photon%icell = icell
  photon%jcell = jcell
  photon%kcell = kcell

  return
  end subroutine raytrace_to_tau_car

end module raytrace
