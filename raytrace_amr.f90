module raytrace_amr_mod
!---------------------------------------------------------------------------
! AMR octree raytrace for MoCafe (dust-only), ported from
! LaRT_v2.00/raytrace_amr.f90 with all Lyman-alpha machinery stripped.
!
!   Positions and path lengths are in CODE UNITS; amr_grid%rhokap(il) is the
!   grey dust opacity per code unit, so tau = rhokap*ds is dimensionless.
!   Octant traversal uses the precomputed neighbor table (amr_next_leaf,
!   O(1) per face crossing).  Face index iface: 1=+x 2=-x 3=+y 4=-y 5=+z 6=-z.
!
! Bound to the raytrace_to_tau / raytrace_to_edge procedure pointers for
! grid_type='amr'; scattering and peeling-off are SHARED with the Cartesian
! dust path (dust scattering is leaf-independent).  See AMR_CLUMPS_PLAN.md
! Part A.  Boundary handling honors par%xy_periodic / xy_symmetry /
! xyz_symmetry (no frequency shift for grey dust); a space-filling octree has
! no refinement gaps, so amr_next_leaf alone suffices for the transport walk.
!---------------------------------------------------------------------------
  use octree_mod
  implicit none
  private

  public :: raytrace_to_tau_amr
  public :: raytrace_to_edge_amr

  real(wp), parameter :: tau_huge = 745.2_wp  ! exp(-tau_huge) ~ 0 in double

contains

  !=========================================================================
  subroutine raytrace_to_tau_amr(photon, grid, tau_in)
    use define
    implicit none
    type(photon_type), intent(inout) :: photon
    type(grid_type),   intent(inout) :: grid   ! unused for AMR; interface conformance
    real(wp),          intent(in)    :: tau_in

    integer  :: il, il_new, icell, iface
    real(wp) :: x, y, z, kx, ky, kz
    real(wp) :: tau, t_exit, d_step, rhokap

    x  = photon%x;    y  = photon%y;    z  = photon%z
    kx = photon%kx;   ky = photon%ky;   kz = photon%kz
    il = photon%icell_amr

    if (il <= 0) then
      il = amr_find_leaf(x, y, z)
      if (il <= 0) then
        photon%inside = .false.
        return
      end if
    end if

    tau = 0.0_wp
    do while (photon%inside)
      icell = amr_grid%icell_of_leaf(il)
      call amr_cell_exit(icell, x, y, z, kx, ky, kz, t_exit, iface)
      rhokap = amr_grid%rhokap(il)

      if (tau + t_exit * rhokap >= tau_in) then
        if (rhokap > 0.0_wp) then
          d_step = (tau_in - tau) / rhokap
        else
          d_step = t_exit
        end if
        x = x + d_step * kx
        y = y + d_step * ky
        z = z + d_step * kz
        tau = tau_in
        exit
      end if

      tau = tau + t_exit * rhokap
      x = x + t_exit * kx
      y = y + t_exit * ky
      z = z + t_exit * kz

      il_new = amr_next_leaf(icell, iface, x, y, z)

      if (il_new <= 0) then
        if (par%xy_periodic) then
          select case (iface)
          case (1); x = x - amr_grid%xrange
          case (2); x = x + amr_grid%xrange
          case (3); y = y - amr_grid%yrange
          case (4); y = y + amr_grid%yrange
          end select
          if (iface <= 4) il_new = amr_find_leaf(x, y, z)
        else if (par%xyz_symmetry) then
          select case (iface)
          case (2);  kx = -kx;  il_new = il
          case (4);  ky = -ky;  il_new = il
          case (6);  kz = -kz;  il_new = il
          end select
        end if
      end if

      if (il_new <= 0) then
        photon%inside = .false.
        exit
      end if
      il = il_new
    end do

    photon%x  = x;   photon%y  = y;   photon%z  = z
    photon%kx = kx;  photon%ky = ky;  photon%kz = kz
    photon%icell_amr = il
  end subroutine raytrace_to_tau_amr

  !=========================================================================
  ! Integrate the dust optical depth from photon0 to the box edge.  photon0
  ! is read-only.  Bound to raytrace_to_edge for the peeling-off estimator.
  !=========================================================================
  subroutine raytrace_to_edge_amr(photon0, grid, tau)
    use define
    implicit none
    type(photon_type), intent(in)  :: photon0
    type(grid_type),   intent(in)  :: grid
    real(wp),          intent(out) :: tau

    integer  :: il, il_new, icell, iface
    real(wp) :: x, y, z, kx, ky, kz, t_exit, rhokap

    x  = photon0%x;    y  = photon0%y;    z  = photon0%z
    kx = photon0%kx;   ky = photon0%ky;   kz = photon0%kz
    il = photon0%icell_amr

    tau = 0.0_wp
    if (il <= 0) then
      il = amr_find_leaf(x, y, z)
      if (il <= 0) return
    end if

    do
      icell  = amr_grid%icell_of_leaf(il)
      call amr_cell_exit(icell, x, y, z, kx, ky, kz, t_exit, iface)
      rhokap = amr_grid%rhokap(il)
      tau = tau + t_exit * rhokap
      if (tau >= tau_huge) return

      x = x + t_exit * kx
      y = y + t_exit * ky
      z = z + t_exit * kz
      il_new = amr_next_leaf(icell, iface, x, y, z)

      if (il_new <= 0) then
        if (par%xy_periodic) then
          select case (iface)
          case (1); x = x - amr_grid%xrange
          case (2); x = x + amr_grid%xrange
          case (3); y = y - amr_grid%yrange
          case (4); y = y + amr_grid%yrange
          end select
          if (iface <= 4) il_new = amr_find_leaf(x, y, z)
        else if (par%xyz_symmetry) then
          select case (iface)
          case (2);  kx = -kx;  il_new = il
          case (4);  ky = -ky;  il_new = il
          case (6);  kz = -kz;  il_new = il
          end select
        end if
      end if

      if (il_new <= 0) return
      il = il_new
    end do
  end subroutine raytrace_to_edge_amr

end module raytrace_amr_mod
