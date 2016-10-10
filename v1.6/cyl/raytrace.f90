module raytrace
contains
!------------------------------
! Find next radial-cell index in cylindrical grid system
! 2013-06 K.-I. Seon
  subroutine next_r_cyl(grid,aa,sqrt_aa,xp,yp,vx,vy,rcell,tdeltr,rstep)
  use define
  implicit none

  type(grid_type), intent(in) :: grid
  real(kind=wp), intent(in)   :: aa,sqrt_aa,xp,yp,vx,vy
  real(kind=wp), intent(out)  :: tdeltr
  integer, intent(inout)      :: rcell
  integer, intent(out)        :: rstep

  real(kind=wp) :: bb,bb2,cc,cc1,det,rr,r_cospv
!--- radial direction
  bb  = xp*vx + yp*vy
  bb2 = bb*bb
  rr  = xp*xp + yp*yp
  r_cospv = -bb/sqrt_aa

  ! Note: aa is greater than or equal to 0.
  if (aa > 0.0_wp) then
     cc  = rr - grid%rface(rcell)**2
     ! Note: cc should be greater than or equal to 0.
     ! However, sometimes cc becomes negative numerically.
     if (cc <= 0.0_wp) cc = 0.0_wp
     det = bb2 - aa*cc

     ! condition for outer circle
     if (r_cospv <= sqrt(cc) .or. rcell == 1) then
        cc  = rr - grid%rface(rcell+1)**2
        det = bb2 - aa*cc
        tdeltr = (-bb+sqrt(det))/aa
        rstep  = 1
     else
     ! condition for inner circle
        tdeltr = (-bb-sqrt(det))/aa
        rstep  = -1
        if (tdeltr <= 0.0_wp) then
           rcell = rcell - 1
           cc1   = rr - grid%rface(rcell)**2
           if (cc1 <= 0.0_wp) cc1 = 0.0_wp
           !if (r_cospv <= sqrt(cc1) .or. rcell == 1) then
           if (r_cospv <= sqrt(cc1)) then
              tdeltr = (-bb+sqrt(det))/aa
              rstep  = 1
           else
              det    = bb2 - aa*cc1
              tdeltr = (-bb-sqrt(det))/aa
           endif
        endif
     endif
  else
     rstep  = 0
     tdeltr = hugest
  endif
  end subroutine next_r_cyl
!------------------------------
! Find next phi-cell index in cylindrical grid system
! 2013-06 K.-I. Seon
  subroutine next_p_cyl(grid,xp,yp,vx,vy,vz,pcell,tdeltp,pstep)
  use define
  implicit none

  type(grid_type), intent(in) :: grid
  real(kind=wp), intent(in) :: xp,yp,vx,vy,vz
  real(kind=wp), intent(out) :: tdeltp
  integer, intent(inout) :: pcell
  integer, intent(out) :: pstep

  real(kind=wp) :: bv1,bv2,pv

!--- phi- direction
  !--- Phi-direction
  ! vector products
  ! v, b, and p represent velocity, boundary, and position vectors.
  ! bv1 : b(i) x v
  ! bv2 : b(i+1) x v
  ! pv  : p x v
  if (grid%np > 1 .and. abs(vz) < 1.0_wp) then
     bv1 = grid%cosp(pcell  )*vy - grid%sinp(pcell  )*vx
     bv2 = grid%cosp(pcell+1)*vy - grid%sinp(pcell+1)*vx
     pv  = xp*vy - yp*vx
     if (bv2 > 0.0_wp .and. pv > 0.0_wp) then
        tdeltp = -(grid%cosp(pcell+1)*yp - grid%sinp(pcell+1)*xp) / bv2
        pstep = 1
     else if (bv1 < 0.0_wp .and. pv < 0.0_wp) then
        tdeltp = -(grid%cosp(pcell)*yp - grid%sinp(pcell)*xp) / bv1
        pstep = -1
        if (tdeltp <= 0.0_wp) then
           pcell = pcell - 1
           if (pcell == 0) pcell = grid%np
           bv1 = grid%cosp(pcell)*vy - grid%sinp(pcell)*vx
           !bug-fixed 2013-06
           !if (bv1 /= 0.0_wp) then
           if (bv1 < 0.0_wp) then
              tdeltp = -(grid%cosp(pcell)*yp - grid%sinp(pcell)*xp) / bv1
              pstep = -1
           else
              tdeltp = hugest
              pstep = 0
           endif
        endif
     else
        tdeltp = hugest
        pstep = 0
     endif
  else
     tdeltp = hugest
     pstep = 0
  endif
  end subroutine next_p_cyl
!------------------------------
! Find next z-cell index in a cartiesian or cylindrical grid system
! 2013-06 K.-I. Seon
! 2015-12-16 tdeltz removed.
  subroutine next_z_cyl(grid,zp,vz,zcell,tmaxz,zstep)
  use define
  implicit none

  type(grid_type), intent(in) :: grid
  real(kind=wp), intent(in) :: zp,vz
  real(kind=wp), intent(out) :: tmaxz
  integer, intent(inout) :: zcell
  integer, intent(out) :: zstep

!--- z direction
  if (vz > 0.0_wp) then
     if (zcell > grid%nz .and. grid%zface(zcell) == zp) return
     zstep  = 1
     tmaxz  = (grid%zface(zcell+1)-zp)/vz
  else if (vz < 0.0_wp) then
     if (grid%zface(zcell) == zp) then
        if (zcell > 1) then
           zcell = zcell - 1
        else
           return
        endif
     endif
     zstep  = -1
     tmaxz  = (grid%zface(zcell)-zp)/vz
  else
     zstep  = 0
     tmaxz  = hugest
  endif
  end subroutine next_z_cyl
!------------------------------
  subroutine raytrace_to_edge_cyl(photon0,grid,tau)
!--- Find the cumulative distance and optical depth to the edge of grid.

  use define
  implicit none

  type(photon_type), intent(in)  :: photon0
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(out) :: tau

! local variables
  integer :: rcell,pcell,zcell
  real(kind=wp) :: xp,yp,zp,vx,vy,vz
  real(kind=wp) :: d
  real(kind=wp) :: rcur,pcur

! local varaibles
  real(kind=wp) :: tmaxr,tmaxp,tmaxz,tdeltr,tdeltp,tdeltz
  integer :: rstep,pstep,zstep
! local variables for radial direction
  real(kind=wp) :: aa,sqrt_aa

!DIR$ ATTRIBUTES INLINE :: raytrace_mp_next_r_cyl
!DIR$ ATTRIBUTES INLINE :: raytrace_mp_next_p_cyl
!DIR$ ATTRIBUTES INLINE :: raytrace_mp_next_z_cyl
!--- xp, yp, zp = the current photon coordinates.
  xp = photon0%x
  yp = photon0%y
  zp = photon0%z
  vx = photon0%vx
  vy = photon0%vy
  vz = photon0%vz
  rcell = photon0%rcell
  pcell = photon0%pcell
  zcell = photon0%zcell

!--- d   = cumulative distance
!--- tau = cumulative optical depth
  tau = 0.0_wp
  d   = 0.0_wp

!--- radial direction
  aa = vx*vx + vy*vy
  sqrt_aa = sqrt(aa)
  call next_r_cyl(grid,aa,sqrt_aa,xp,yp,vx,vy,rcell,tmaxr,rstep)

!--- z direction
  call next_z_cyl(grid,zp,vz,zcell,tmaxz,zstep)

!--- integrate through grid
  do while(.true.)
    if (tmaxr < tmaxz) then
       tau   = tau + (tmaxr - d)*grid%opacity(rcell,pcell,zcell)
       d     = tmaxr
       rcell = rcell + rstep
       if (rcell > grid%nr) exit
       xp = photon0%x + d*vx
       yp = photon0%y + d*vy
       call next_r_cyl(grid,aa,sqrt_aa,xp,yp,vx,vy,rcell,tdeltr,rstep)
       tmaxr = tmaxr + tdeltr
    else
       tau   = tau + (tmaxz - d)*grid%opacity(rcell,pcell,zcell)
       d     = tmaxz
       zcell = zcell + zstep
       if (zcell < 1 .or. zcell > grid%nz) exit
!!
       tdeltz = grid%dz(zcell)/abs(vz)
       tmaxz = tmaxz + tdeltz
       if (tmaxr == tmaxz) then
          rcell = rcell + rstep
          if (rcell > grid%nr) exit
          xp = photon0%x + d*vx
          yp = photon0%y + d*vy
          call next_r_cyl(grid,aa,sqrt_aa,xp,yp,vx,vy,rcell,tdeltr,rstep)
          tmaxr = tmaxr + tdeltr
       endif
    endif
  enddo

  return
  end subroutine raytrace_to_edge_cyl
  subroutine raytrace_to_tau_cyl(photon,grid,tau_in,inside_grid)
!--- Find the cumulative distance and optical depth to the edge of grid.

  use define
  implicit none

  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(in)    :: grid
  real(kind=wp),     intent(in)    :: tau_in
  logical,           intent(inout) :: inside_grid

! local variables
  integer :: rcell,pcell,zcell
  real(kind=wp) :: xp,yp,zp,vx,vy,vz
  real(kind=wp) :: tau
  real(kind=wp) :: d_overshoot,d
  real(kind=wp) :: rcur,pcur

! local varaibles
  real(kind=wp) :: tmaxr,tmaxp,tmaxz,tdeltr,tdeltp,tdeltz
  integer :: rstep,pstep,zstep
! local variables for radial direction
  real(kind=wp) :: aa,sqrt_aa

!DIR$ ATTRIBUTES INLINE :: raytrace_mp_next_r_cyl
!DIR$ ATTRIBUTES INLINE :: raytrace_mp_next_p_cyl
!DIR$ ATTRIBUTES INLINE :: raytrace_mp_next_z_cyl
!--- xp, yp, zp = the current photon coordinates.
  xp = photon%x
  yp = photon%y
  zp = photon%z
  vx = photon%vx
  vy = photon%vy
  vz = photon%vz
  rcell = photon%rcell
  pcell = photon%pcell
  zcell = photon%zcell

!--- inside_grid = .true. means photon is in grid envelope
  inside_grid = .true.

!--- d   = cumulative distance
!--- tau = cumulative optical depth
  tau = 0.0_wp
  d   = 0.0_wp

!--- radial direction
  aa      = vx*vx + vy*vy
  sqrt_aa = sqrt(aa)
  rcur    = sqrt(xp*xp + yp*yp)
  call next_r_cyl(grid,aa,sqrt_aa,xp,yp,vx,vy,rcell,tmaxr,rstep)
  if (aa > 0.0_wp .and. grid%rface(grid%nr+1) == rcur) then
    inside_grid = .false.
    return
  endif

!--- z direction
  call next_z_cyl(grid,zp,vz,zcell,tmaxz,zstep)
  if ((vz > 0.0_wp .and. grid%zface(grid%nz+1) == zp) .or. &
      (vz < 0.0_wp .and. grid%zface(1) == zp)) then
    inside_grid = .false.
    return
  endif

!--- integrate through grid
  do while(inside_grid)
    if (tmaxr < tmaxz) then
       tau   = tau + (tmaxr - d)*grid%opacity(rcell,pcell,zcell)
       d     = tmaxr
       if (tau >= tau_in) then
          if (grid%opacity(rcell,pcell,zcell) > 0.0) then
             d_overshoot = (tau - tau_in)/grid%opacity(rcell,pcell,zcell)
             d = d - d_overshoot
          endif
          xp = photon%x + d*vx
          yp = photon%y + d*vy
          zp = photon%z + d*vz
          exit
       endif
       rcell = rcell + rstep
       if (rcell > grid%nr) then
          inside_grid = .false.
          exit
       endif
       xp = photon%x + d*vx
       yp = photon%y + d*vy
       call next_r_cyl(grid,aa,sqrt_aa,xp,yp,vx,vy,rcell,tdeltr,rstep)
       tmaxr = tmaxr + tdeltr
    else
       tau   = tau + (tmaxz - d)*grid%opacity(rcell,pcell,zcell)
       d     = tmaxz
       if (tau >= tau_in) then
          if (grid%opacity(rcell,pcell,zcell) > 0.0) then
             d_overshoot = (tau - tau_in)/grid%opacity(rcell,pcell,zcell)
             d = d - d_overshoot
          endif
          xp = photon%x + d*vx
          yp = photon%y + d*vy
          zp = photon%z + d*vz
          exit
       endif
       zcell = zcell + zstep
       if (zcell < 1 .or. zcell > grid%nz) then
          inside_grid = .false.
          exit
       endif
!!
       tdeltz = grid%dz(zcell)/abs(vz)
       tmaxz = tmaxz + tdeltz
       if (tmaxr == tmaxz) then
          rcell = rcell + rstep
          if (rcell > grid%nr) then
             inside_grid = .false.
             exit
          endif
          xp = photon%x + d*vx
          yp = photon%y + d*vy
          call next_r_cyl(grid,aa,sqrt_aa,xp,yp,vx,vy,rcell,tdeltr,rstep)
          tmaxr = tmaxr + tdeltr
       endif
    endif
  enddo

!--- xp, yp, zp = the current photon coordinates.
  photon%x  = xp
  photon%y  = yp
  photon%z  = zp
  photon%rcell = rcell
  photon%pcell = pcell
  photon%zcell = zcell

  return
  end subroutine raytrace_to_tau_cyl
end module raytrace
