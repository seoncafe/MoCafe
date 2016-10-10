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
  use raytrace, only: raytrace_to_edge_cyl, raytrace_to_tau_cyl
  use detect

!--- Parameters
  implicit none

  integer, parameter  :: nchunk = 8
  type(photon_type)   :: photon
  type(grid_type)     :: grid
  type(observer_type) :: observer
  type(output_type)   :: output(nchunk),output_sum

  integer(kind=i8b)   :: ichunk,my_nphotons(nchunk)
  integer(kind=i8b) :: ip
  real(kind=wp)     :: totscatt
  integer           :: nscatt
  real(kind=wp)     :: tau,tau0
  integer           :: i, status
  real(kind=wp)     :: time1,time2
  logical           :: inside_grid

!--- Initial set up
  call cpu_time(time1)
  call setup(photon,grid,observer)
  call init_random_seed()
  status = -999
  totscatt = 0
  par%nthreads = nchunk

  do i=1,nchunk
     call output_create(output(i))
     my_nphotons(i) = 0
  enddo

!--- Photon loop
  do ip=1,par%nphotons
     ichunk = (ip-1)/(par%nphotons/nchunk)+1
     if (ichunk > nchunk) ichunk = ichunk - nchunk
     my_nphotons(ichunk) = my_nphotons(ichunk) + 1
     ! this is to compare with parallel version.

     !--- Release photon
     call generate_photon(photon,grid)

     !--- Add direct photon
     call add_direct(photon,grid,observer,output(ichunk))

     nscatt = 0
     inside_grid = .true.

     do while (inside_grid)
        if (nscatt == 0) then
           !--- Find optical depth (tau0) to edge of grid
           !--- Force photon to scatter at optical depth tau before edge of grid
           call raytrace_to_edge_cyl(photon,grid,tau0)
           photon%wgt = 1.0_wp - exp(-tau0)
           tau = -log(1.0_wp - rand_number()*photon%wgt)
        else
           tau = -log(rand_number())
        endif

        !--- Find scattering location of tau
        call raytrace_to_tau_cyl(photon,grid,tau,inside_grid)

        if (inside_grid) then
           photon%wgt = photon%wgt * par%albedo
           call peelingoff(photon,grid,observer,output(ichunk))
           call scatt(photon)
           nscatt = nscatt + 1
        endif
     enddo

     totscatt = totscatt + nscatt
     if (mod(ip,par%nprint) == 0 .or. ip == par%nphotons) then
        call cpu_time(time2)
        write(6,'(i12,a,f8.3,a)') ip,' photons completed: ',(time2-time1),' secs'
     endif
  enddo

  call grid_destroy(grid)
  call output_create(output_sum)
  do i=1,nchunk
     if (i/=nchunk) call output_gather_serial(output(i),output_sum,my_nphotons(i),do_final=.false.)
     if (i==nchunk) call output_gather_serial(output(i),output_sum,my_nphotons(i),do_final=.true.)
  enddo

!--- Put results into output files in units of "input luminosity unit"/cm^2/sr
  call write_output(par%out_file,output_sum)
  call output_destroy(output_sum)

  call cpu_time(time2)
  write(6,'(a,f7.2)')   'Average Number of scattering: ',dble(totscatt)/dble(par%nphotons)
  write(6,'(a,f8.3,a)') 'Total Excution Time: ', (time2-time1),' secs'

  stop
  end
