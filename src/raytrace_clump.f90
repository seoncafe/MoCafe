module raytrace_clump_mod
!---------------------------------------------------------------------------
! Raytrace routines for the clumpy dust medium (par%grid_type = 'clump').
!
! The medium consists of N spherical dust clumps (radii cl_radius(:)) placed
! inside a sphere of radius sphere_R.  Outside the clumps the medium is vacuum.
! Each clump has its own grey dust opacity cl_rhokap(icl) (extinction per unit
! code length).  Because the opacity is grey there is NO frequency, Voigt
! profile, Doppler width, or velocity -- the clump subsystem is purely
! geometric plus a scalar opacity for each clump (LaRT's Lyman-alpha transfer is
! stripped; see AMR_CLUMPS_PLAN.md Part B).
!
!   photon%icell_clump: 0 = in vacuum; > 0 = index of the current clump.
!   A photon alternates between vacuum flight (find the next clump it enters)
!   and clump traversal (accumulate kappa_c * chord).
!---------------------------------------------------------------------------
  use define, only: wp
  use clump_mod
  implicit none
  private

  public :: raytrace_to_tau_clump
  public :: raytrace_to_edge_clump

contains

  !===========================================================================
  ! Internal helper: update photon%icell/jcell/kcell from position (keeps the
  ! Cartesian cell indices consistent for any code that reads them).
  !===========================================================================
  pure subroutine update_cell_idx(photon, grid)
  use define
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(in)    :: grid
  photon%icell = max(1, min(grid%nx, floor((photon%x - grid%xmin)/grid%dx) + 1))
  photon%jcell = max(1, min(grid%ny, floor((photon%y - grid%ymin)/grid%dy) + 1))
  photon%kcell = max(1, min(grid%nz, floor((photon%z - grid%zmin)/grid%dz) + 1))
  end subroutine update_cell_idx
  !===========================================================================

  !===========================================================================
  subroutine raytrace_to_tau_clump(photon, grid, tau_in)
  !---------------------------------------------------------------------------
  ! Advance photon through the clump medium until optical depth tau_in is
  ! accumulated or the photon leaves the sphere.
  ! Updates: photon%x/y/z, photon%icell_clump, photon%inside.
  !---------------------------------------------------------------------------
  use define
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  real(kind=wp),     intent(in)    :: tau_in

  real(kind=wp)  :: kx, ky, kz, tau_rem, kap, ds, t_seg, t_sp
  real(kind=wp)  :: te, tx2
  integer(int64) :: icl, icl_found, last_icl
  logical        :: found

  kx = photon%kx;  ky = photon%ky;  kz = photon%kz
  tau_rem  = tau_in
  last_icl = 0_int64   ! index of most-recently exited clump (avoid re-entry)

  do while (photon%inside)

     if (photon%icell_clump > 0) then
        !--- photon is inside clump icl
        icl   = int(photon%icell_clump, int64)
        t_seg = clump_exit_dist(photon%x, photon%y, photon%z, kx, ky, kz, icl)
        kap   = cl_rhokap(icl)

        if (tau_rem <= kap * t_seg) then
           !--- scatter inside this clump
           ds = tau_rem / max(kap, tiny(1.0_wp))
           photon%x = photon%x + ds * kx
           photon%y = photon%y + ds * ky
           photon%z = photon%z + ds * kz
           call update_cell_idx(photon, grid)
           return
        end if

        !--- traverse through this clump and exit into vacuum
        tau_rem  = tau_rem - kap * t_seg
        photon%x = photon%x + t_seg * kx
        photon%y = photon%y + t_seg * ky
        photon%z = photon%z + t_seg * kz

        last_icl           = icl
        photon%icell_clump = 0

        !--- check sphere exit
        if (photon%x**2 + photon%y**2 + photon%z**2 >= sphere_R**2) then
           photon%inside = .false.
           call update_cell_idx(photon, grid)
           return
        end if

     else
        !--- photon in vacuum: find next clump or sphere exit
        t_sp = sphere_exit_dist(photon%x, photon%y, photon%z, kx, ky, kz)
        if (t_sp <= 0.0_wp) then
           photon%inside = .false.
           return
        end if

        call find_next_clump(photon%x, photon%y, photon%z, kx, ky, kz, &
             last_icl, t_sp, te, tx2, icl_found, found)

        if (found) then
           !--- advance to clump entry (vacuum, no tau accumulated)
           te = max(0.0_wp, te)
           photon%x = photon%x + te * kx
           photon%y = photon%y + te * ky
           photon%z = photon%z + te * kz

           last_icl           = 0_int64   ! entering a new clump; no skip needed
           photon%icell_clump = int(icl_found)
           call update_cell_idx(photon, grid)
        else
           !--- no clump on this ray: skip to sphere boundary and exit
           photon%x = photon%x + t_sp * kx
           photon%y = photon%y + t_sp * ky
           photon%z = photon%z + t_sp * kz
           photon%inside = .false.
           call update_cell_idx(photon, grid)
           return
        end if

     end if

  end do

  end subroutine raytrace_to_tau_clump
  !===========================================================================

  !===========================================================================
  subroutine raytrace_to_edge_clump(photon0, grid, tau)
  !---------------------------------------------------------------------------
  ! Accumulate total dust optical depth from photon0 to the sphere edge along
  ! the photon direction.  photon0 is read-only.  Bound to the raytrace_to_edge
  ! procedure pointer in clump mode (used by the shared peeling-off routines).
  !---------------------------------------------------------------------------
  use define
  implicit none
  type(photon_type), intent(in)  :: photon0
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(out) :: tau

  real(kind=wp)  :: xp, yp, zp, kx, ky, kz
  real(kind=wp)  :: t_sp, t_seg, te, tx2
  integer(int64) :: icl_cur, icl_found
  logical        :: found

  tau     = 0.0_wp
  xp      = photon0%x;  yp = photon0%y;  zp = photon0%z
  kx      = photon0%kx; ky = photon0%ky; kz = photon0%kz
  icl_cur = int(photon0%icell_clump, int64)

  !--- first, if photon starts inside a clump, traverse it
  if (icl_cur > 0) then
     t_seg = clump_exit_dist(xp, yp, zp, kx, ky, kz, icl_cur)
     tau   = tau + cl_rhokap(icl_cur) * t_seg
     xp    = xp + t_seg * kx
     yp    = yp + t_seg * ky
     zp    = zp + t_seg * kz
     ! keep icl_cur as skip index for the first find_next_clump call
     if (xp**2 + yp**2 + zp**2 >= sphere_R**2) return
  end if

  !--- DDA through the remaining clumps to the sphere boundary
  do
     t_sp = sphere_exit_dist(xp, yp, zp, kx, ky, kz)
     if (t_sp <= 0.0_wp) exit

     call find_next_clump(xp, yp, zp, kx, ky, kz, icl_cur, t_sp, &
                           te, tx2, icl_found, found)
     if (.not. found) exit

     !--- advance to clump entry, accumulate tau through clump
     te = max(0.0_wp, te)
     xp = xp + te * kx;  yp = yp + te * ky;  zp = zp + te * kz

     t_seg = clump_exit_dist(xp, yp, zp, kx, ky, kz, icl_found)
     tau   = tau + cl_rhokap(icl_found) * t_seg
     xp    = xp + t_seg * kx;  yp = yp + t_seg * ky;  zp = zp + t_seg * kz

     icl_cur = icl_found  ! skip this clump in next search (just exited)
     if (xp**2 + yp**2 + zp**2 >= sphere_R**2) exit
  end do

  ! grid is intent(in) for procedure-pointer interface conformance; the clump
  ! raytrace needs no grid lookups (opacity lives in the clump arrays).
  if (.false.) tau = tau + grid%dx*0.0_wp

  end subroutine raytrace_to_edge_clump
  !===========================================================================

end module raytrace_clump_mod
