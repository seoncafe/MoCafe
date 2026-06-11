module grid_mod_clump
!---------------------------------------------------------------------------
! Grid setup for the clumpy dust medium (par%grid_type = 'clump').
!
! A sphere of radius par%rmax holds N non-overlapping spherical dust clumps of
! radius par%clump_radius; the inter-clump medium is vacuum.  Approach: build
! the clump population (init_clumps), build the Cartesian bounding-box grid
! (grid_create) only so the observer/output machinery runs unchanged, then
! zero grid%rhokap (all opacity lives inside the clumps) and set the realized
! system-level dust scalars from the clump arrays.
!
! Dust-only slim of LaRT_v2.00/grid_mod_clump.f90: MoCafe's grid_type has no
! Dfreq/voigt_a/vf fields, so the frequency/velocity grid fills are dropped.
!---------------------------------------------------------------------------
  use grid_mod
  use clump_mod
  use iofile_mod, only: io_file_extension
  implicit none

  public :: grid_create_clump, grid_destroy_clump

contains

  !===========================================================================
  subroutine grid_create_clump(grid)
  !---------------------------------------------------------------------------
  ! Set up the Cartesian shell grid and initialize the clump population.
  !---------------------------------------------------------------------------
  use define
  use mpi
  implicit none
  type(grid_type), intent(inout) :: grid

  integer            :: ierr
  character(len=128) :: saved_distance_unit
  real(kind=wp)      :: saved_distance2cm

  if (par%rmax <= 0.0_wp) then
     if (mpar%p_rank == 0) write(*,*) 'ERROR: par%rmax must be > 0 for clump medium'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if

  !--- Force the Cartesian box to enclose the sphere [-rmax, +rmax]^3.
  par%xmax = par%rmax
  par%ymax = par%rmax
  par%zmax = par%rmax
  if (par%nx < 1) par%nx = 11
  if (par%ny < 1) par%ny = 11
  if (par%nz < 1) par%nz = 11
  par%xyz_symmetry = .false.
  par%xy_periodic  = .false.
  par%z_symmetry   = .false.

  !--- Initialize clumps FIRST: this consumes the system-level target
  !    (par%taumax / par%tauhomo) and par%distance2cm to back-solve the dust
  !    opacity, places clumps via RSA, and builds the CSR acceleration grid.
  call init_clumps(par%rmax)

  !--- grid_create sets up the bounding box (faces, dx, xmin, ...) and a
  !    placeholder rhokap, but it also (a) resets par%distance_unit/distance2cm
  !    to ''/1.0 when no density file is present and (b) recomputes
  !    par%tauhomo/taumax from that placeholder grid.  Save/restore the
  !    distance unit and recompute the clump scalars afterwards so the
  !    clump-derived values stand.
  saved_distance_unit = par%distance_unit
  saved_distance2cm   = par%distance2cm
  call grid_create(grid)
  par%distance_unit   = saved_distance_unit
  par%distance2cm     = saved_distance2cm

  !--- Opacity lives only inside clumps; the Cartesian grid is vacuum.
  if (mpar%h_rank == 0) grid%rhokap(:,:,:) = 0.0_wp
  call MPI_BARRIER(mpar%hostcomm, ierr)

  !--- Realized system-level dust scalars from the clump population.
  call compute_clump_scalars(par%tauhomo, par%taumax)

  if (mpar%p_rank == 0) then
     write(*,'(a,es14.5)') ' Clump derived: tauhomo = ', par%tauhomo
     write(*,'(a,es14.5)') ' Clump derived: taumax  = ', par%taumax
     write(*,'(a,f12.5)')  ' Clump grid: rmax     = ', par%rmax
     write(*,'(a,3i5)')    ' Clump grid: nx,ny,nz = ', grid%nx, grid%ny, grid%nz
  end if

  !--- Optionally save the clump population for reuse (clump_input_file).
  if (par%save_clump_info .and. mpar%p_rank == 0) then
     call write_clumps_info(trim(par%base_name)//'_clumps'// &
                            trim(io_file_extension(par%file_format)))
  end if

  end subroutine grid_create_clump
  !===========================================================================

  !===========================================================================
  subroutine grid_destroy_clump(grid)
  !---------------------------------------------------------------------------
  ! grid_destroy already frees ALL shared-memory windows (including the clump
  ! and CSR arrays) via destroy_shared_mem_all(), so we only nullify the clump
  ! pointers here -- calling destroy_clumps() would invoke a second collective
  ! free of already-freed windows.
  !---------------------------------------------------------------------------
  use define
  implicit none
  type(grid_type), intent(inout) :: grid
  call grid_destroy(grid)
  nullify(cl_x, cl_y, cl_z, cl_radius, cl_radius2, cl_rhokap, cg_start, cg_list)
  N_clumps = 0_int64
  end subroutine grid_destroy_clump
  !===========================================================================

end module grid_mod_clump
