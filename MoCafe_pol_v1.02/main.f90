  program main

  use define
  use setup_mod
  use observer_mod
  use grid_mod
  use output_sum
  use peelingoff_mod
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
  call grid_create(grid)
  call observer_create()
  if (par%sightline_tau) call sightline_tau(grid)

  !--- Run Main Calculation
  call time_stamp(dtime)
  if (mpar%p_rank == 0) write(6,'(a,f8.3,a)') '---> Now starting simulation...  @ ', dtime/60.0_wp, ' mins'
  call run_simulation(grid)

  !--- Final Output
  if (mpar%p_rank == 0) write(6,'(a)') 'Now Gathering Results...'
  call output_reduce(grid)

  call time_stamp(dtime)
  if (mpar%p_rank == 0) then
     call output_normalize(grid)
     par%exetime = dtime/60.0_wp
     call write_output(trim(par%out_file),grid)

     write(6,'(a,es12.4)')  'Average Number of scattering : ', par%nscatt_tot
     write(6,'(a,f8.3,a)')  'Total Excution Time          : ', par%exetime,' mins'
  endif

  call grid_destroy(grid)
  call observer_destroy()
  call MPI_FINALIZE(ierr)
  stop
  end
