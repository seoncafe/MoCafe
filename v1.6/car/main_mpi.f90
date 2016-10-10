  program main

!--- Modules
  use define
  use random
  use scatter
  use photon_mod
  use fits_module
  use grid_mod
  use setup_mod
  use output_mod
  use mathlib
  use raytrace, only: raytrace_to_edge_car, raytrace_to_tau_car
  use detect
  use mpi
!TEST
!  use write2D

!--- Parameters
  implicit none

  type(photon_type)   :: photon
  type(grid_type)     :: grid
  type(observer_type) :: observer
  type(output_type)   :: output, output_sum

  integer(kind=i8b) :: ip
  integer           :: nscatt
  real(kind=wp)     :: totscatt
  real(kind=wp)     :: tau,tau0
  integer           :: i, status
  real(kind=dp)     :: time1,time2
  logical           :: inside_grid
!TEST
!  character(len=100) :: fname

!--- Variables for MPI
  integer :: ierr,myid,nsize
  !integer :: nsize
  integer(kind=i8b) :: remainder,my_nphotons
  real(kind=wp)     :: my_totscatt

!--- Initial set up
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid,     ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, par%nthreads, ierr)

!--- Initial set up
  time1 = MPI_WTIME()

  call setup(photon,grid,observer)
  status   = -999

  call output_create(output)
  call init_random_seed(par%iseed)

!--- Photon loop
  my_nphotons = int(par%nphotons/par%nthreads)
  remainder   = mod(par%nphotons,par%nthreads)
  if (myid < remainder) my_nphotons = my_nphotons + 1

  do ip=1,my_nphotons
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

     my_totscatt = my_totscatt + nscatt
     if (myid == 0 .and. mod(ip,int(par%nprint/par%nthreads)) == 0) then
        time2 = MPI_WTIME()
        write(6,'(i12,a,es13.5,a)') par%nprint*int(ip/int(par%nprint/par%nthreads)),' photons calculated: ',(time2-time1),' sec'
     endif
  enddo

  call grid_destroy(grid)
!----TEST
!  write(fname,'(a,i2.2,a)') 'out',myid,'.fits.gz'
!  call write_2D(fname,output%tot)
!----
  if (myid == 0) call output_create(output_sum)
  call output_gather_mpi(output,output_sum,my_nphotons)
  totscatt = 0.0_wp
  call MPI_REDUCE(my_totscatt, totscatt, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if (myid == 0) then
     !--- Put results into output files in units of "input luminosity unit"/cm^2/sr
     call write_output(par%out_file,output_sum)
     call output_destroy(output_sum)
     time2 = MPI_WTIME()
     write(6,'(a,f7.2)')   'Average Number of scattering: ',totscatt/dble(par%nphotons)
     write(6,'(a,f8.3,a)') 'Total Excution Time: ', (time2-time1),' secs'
  endif

  call MPI_FINALIZE(ierr)

  stop
  end
