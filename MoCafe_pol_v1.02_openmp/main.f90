  program main

  use define
  use utility
  use setup_mod
  use observer_mod
  use grid_mod
  use output_sum
  use peelingoff_mod

  !--- Parameter
  implicit none
  type(grid_type) :: grid
  real(kind=wp)   :: dtime

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

  call time_stamp(dtime)
  call output_normalize(grid)
  par%exetime = dtime/60.0_wp
  call write_output(trim(par%out_file),grid)

  write(6,'(a,es12.4)')  'Average Number of scattering : ', par%nscatt_tot
  write(6,'(a,f8.3,a)')  'Total Excution Time          : ', par%exetime,' mins'

  call grid_destroy(grid)
  call observer_destroy()
  stop
  end
