  program main

  use define
  use setup_mod
  use observer_mod
  use grid_mod
  use grid_mod_clump
  use grid_mod_amr
  use sightline_tau_mod
  use output_sum
  use peelingoff_mod
  use jtally_mod,   only : jtally_setup, jtally_reduce, jtally_write, jt_on
  use dustemis_mod, only : setup_dustemis, compute_dustemis, write_dustemis
  use lucy_mod,     only : run_lucy_iteration
  use bw_mod,       only : setup_bw, run_bw, bw_finalize, bw_Tmap, bw_write
  use allsky_mod,   only : setup_allsky, allsky_reduce, allsky_write, allsky_on
  use utility
  use mpi

  !--- Parameter
  implicit none
  type(grid_type)   :: grid
  real(kind=wp)     :: dtime
  logical           :: use_bw

  !--- for MPI
  integer :: ierr

  call MPI_INIT(ierr)
  call time_stamp(dtime)

  !--- Initial set up
  call read_input()
  call setup_procedure()
  select case (trim(par%grid_type))
  case ('clump')
     call grid_create_clump(grid)
  case ('amr')
     call grid_create_amr(grid)
  case default
     call grid_create(grid)
  end select
  call observer_create()
  if (par%sightline_tau) call make_sightline_tau(grid)
  use_bw = par%use_dustemis .and. trim(par%dust_emission_method) == 'bw01'
  if (par%save_jlam)     call jtally_setup(grid)
  if (par%use_dustemis .and. .not. use_bw) call setup_dustemis(grid)
  if (use_bw)            call setup_bw(grid)
  if (par%allsky)        call setup_allsky()

  !--- Dust emission (Stage 3, Mode 1 Lucy): energy iterations first, which
  !--- tally J_lambda and converge the per-cell emission (dust self-absorption
  !--- if par%dust_niter > 1).  Leaves jt_sum = converged total J and turns the
  !--- tally off so the subsequent imaging pass does not overwrite it.
  if (par%use_dustemis .and. .not. use_bw) then
     call time_stamp(dtime)
     if (mpar%p_rank == 0) write(6,'(a,f8.3,a)') '---> Lucy energy iterations...  @ ', dtime/60.0_wp, ' mins'
     call run_lucy_iteration(grid)
  endif

  !--- Run Main Calculation.  Bjorkman & Wood (Mode 2) replaces the standard
  !--- imaging pass with its own scatter/absorb/reemit transport; otherwise the
  !--- normal stellar imaging pass runs (peel on).
  call time_stamp(dtime)
  if (mpar%p_rank == 0) write(6,'(a,f8.3,a)') '---> Now starting simulation...  @ ', dtime/60.0_wp, ' mins'
  if (use_bw) then
     call run_bw(grid)
  else
     call run_simulation(grid)
  endif

  !--- Final Output
  if (mpar%p_rank == 0) write(6,'(a)') 'Now Gathering Results...'
  call output_reduce(grid)
  if (use_bw) then
     call bw_finalize()
     call bw_Tmap(grid)
     call bw_write(grid)
  else if (par%use_dustemis) then
     !--- emission already computed/converged inside run_lucy_iteration.
     call write_dustemis(grid)
     if (par%save_jlam) call jtally_write(grid)   ! converged total J
  else if (par%save_jlam) then
     call jtally_reduce()
     call jtally_write(grid)
  endif
  if (par%allsky) then
     call allsky_reduce()
     call allsky_write()
  endif

  call time_stamp(dtime)
  if (mpar%p_rank == 0) then
     call output_normalize(grid)
     par%exetime = dtime/60.0_wp
     call write_output(trim(par%out_file),grid)

     write(6,'(a,es12.4)')  'Average Number of scattering : ', par%nscatt_tot
     write(6,'(a,f8.3,a)')  'Total Excution Time          : ', par%exetime,' mins'
     write(6,'(2a)')        ' >>> STOP  @ ', get_date_time()
  endif

  call grid_destroy(grid)
  call observer_destroy()
  call MPI_FINALIZE(ierr)
  stop
  end
