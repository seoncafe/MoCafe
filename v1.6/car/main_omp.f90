  program main

  use define
  use random
  use scatter
  use photon_mod
  use fits_module
  use setup_mod
  use grid_mod
  use output_mod
  use mathlib
  use raytrace, only: raytrace_to_edge_car, raytrace_to_tau_car
  use detect
  use omp_lib

!--- Parameters
  implicit none
  type(photon_type)   :: photon
  type(grid_type)     :: grid
  type(observer_type) :: observer
  type(output_type)   :: output, output_sum

  integer(kind=i8b) :: ip,my_nphotons
  integer           :: nscatt
  real(kind=wp)     :: totscatt
  real(kind=wp)     :: tau,tau0
  integer           :: i, status
  real(kind=dp)     :: time1,time2
  logical           :: inside_grid

!--- Initial set up
  time1 = omp_get_wtime()
  call setup(photon,grid,observer)
  status = -999
  totscatt = 0

  !$OMP PARALLEL &
  !$OMP default(shared) &
  !$OMP firstprivate(photon) &
  !$OMP private(my_nphotons,output)
  call output_create(output)
  par%nthreads = omp_get_num_threads()
  call init_random_seed()
  my_nphotons = 0

!--- Photon loop
  !$OMP DO &
  !$OMP private(tau0,tau,nscatt,inside_grid) &
  !$OMP reduction(+:totscatt) &
  !$OMP schedule(dynamic,1)
  do ip=1,par%nphotons
     my_nphotons = my_nphotons + 1
     !--- Release photon
     call generate_photon(photon,grid)

     !--- Add direct photon
     call add_direct(photon,grid,observer,output)

     nscatt = 0
     inside_grid = .true.

     do while (inside_grid)
        if (nscatt == 0) then
           !--- Find optical depth (tau0) to edge of grid
           !--- Force photon to scatter at optical depth tau before edge of grid
           call raytrace_to_edge_car(photon,grid,tau0)
           photon%wgt = 1.0_wp - exp(-tau0)
           tau = -log(1.0_wp - rand_number()*photon%wgt)
        else
           tau = -log(rand_number())
        endif

        !--- Find scattering location of tau
        call raytrace_to_tau_car(photon,grid,tau,inside_grid)

        if (inside_grid) then
           photon%wgt = photon%wgt * par%albedo
           call peelingoff(photon,grid,observer,output)
           call scatt(photon)
           nscatt = nscatt + 1
        endif
     enddo

     totscatt = totscatt + nscatt
     if (mod(ip,par%nprint) == 0 .or. ip == par%nphotons) then
        time2 = omp_get_wtime()
        write(6,'(i12,a,es13.5,a)') ip,' packets completed: ',(time2-time1),' secs'
     endif
  enddo
  !$OMP END DO
  !$OMP SINGLE
  call grid_destroy(grid)
  call output_create(output_sum)
  !$OMP END SINGLE
  call output_gather_omp(output,output_sum,my_nphotons)
  !$OMP END PARALLEL

!--- Put results into output files in units of "input luminosity unit"/cm^2/sr
  call write_output(par%out_file,output_sum)
  call output_destroy(output_sum)

  time2 = omp_get_wtime()
  write(6,'(a,f7.2)')   'Average Number of scattering: ',dble(totscatt)/dble(par%nphotons)
  write(6,'(a,f8.3,a)') 'Total Excution Time: ', (time2-time1),' secs'

  stop
  end
