  program main

  use define
  use utility
  use setup_mod
  use observer_mod
  use grid_mod
  use output_sum
  use peelingoff_mod
#ifdef MPI
  use mpi
#endif

  !--- Parameter
  implicit none
  type(grid_type) :: grid
  real(kind=wp)   :: dtime

#ifdef MPI
  !--- for MPI
  integer :: ierr
#endif

  !--- Initialize MPI + Openmp
  !--- call MPI_INIT_THREAD(MPI_THREAD_FUNNELED, MPI_provided,ierr)
  !--- Please use MPI_THREAD_FUNNELED if MPI_provided < MPI_THREAD_MULTIPLE (comment added on 2020-11-08).
  !--- Then, par%use_master_slave will be set to false in setup_MPI_openMP_parameters.
#ifdef MPI
  call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, mpar%MPI_provided,ierr)
#endif
  call time_stamp(dtime)

  !--- Initial set up
  call read_input()
  call setup_procedure()

  !===== Start OPENMP
  !--- It is safe to do global settings outside of the parallel construct. (comment added, 2020.10.19).
  !--- Warning: Don't even try to parallelize global settings. You will be messed up with the settings.
  call grid_create(grid)
  call observer_create()
  !$OMP PARALLEL default(shared)
  if (par%sightline_tau) call sightline_tau(grid)
  !$OMP END PARALLEL

  !--- Run Main Calculation
  !$OMP PARALLEL default(shared)
  call run_simulation(grid)
  !$OMP END PARALLEL

  !--- Final Output
#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
  if (mpar%rank == 0) write(6,'(a)') 'Now Gathering Results...'
  call output_reduce(grid)

  call time_stamp(dtime)
  if (mpar%rank == 0) then
     call output_normalize(grid)
     par%exetime = dtime/60.0_wp
     call write_output(trim(par%out_file),grid)

     write(6,'(a,es12.4)')  'Average Number of scattering : ', par%nscatt_tot
     write(6,'(a,f8.3,a)')  'Total Excution Time          : ', par%exetime,' mins'
  endif

  call grid_destroy(grid)
  call observer_destroy()
#ifdef MPI
  call MPI_FINALIZE(ierr)
#endif MPI
  stop
  end
