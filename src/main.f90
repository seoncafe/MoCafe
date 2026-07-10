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
  use jtally_mod,   only : jtally_setup, jtally_reduce, jtally_write
  use dustemis_mod, only : setup_dustemis, compute_dustemis, write_dustemis
  use utility
  use mpi

  !--- Parameter
  implicit none
  type(grid_type)   :: grid
  real(kind=wp)     :: dtime

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
  if (par%save_jlam)     call jtally_setup(grid)
  if (par%use_dustemis)  call setup_dustemis(grid)

  !--- Run Main Calculation
  call time_stamp(dtime)
  if (mpar%p_rank == 0) write(6,'(a,f8.3,a)') '---> Now starting simulation...  @ ', dtime/60.0_wp, ' mins'
  call run_simulation(grid)

  !--- Final Output
  if (mpar%p_rank == 0) write(6,'(a)') 'Now Gathering Results...'
  call output_reduce(grid)
  if (par%save_jlam) then
     call jtally_reduce()
     !--- dust emission (Stage 3) uses the RAW jt_sum, so compute it before
     !--- jtally_write converts jt_sum in place to J_lambda.
     if (par%use_dustemis) then
        call compute_dustemis(grid)
        call write_dustemis(grid)
     endif
     call jtally_write(grid)
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
