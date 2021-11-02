module setup_mod
contains
!+++++++++++++++++++++++++++++++++++++++++++
  subroutine read_input
  use define
  use utility
  use mpi
  implicit none

! local variables
  character(len=128) :: model_infile, arg
  integer :: unit, ierr, status

!--- Read in parameters from params.par using namelist command
  namelist /parameters/ par

! read parameters
  if (command_argument_count() >= 1) then
     call get_command_argument(1, model_infile)
  else
     call get_command_argument(0, arg)
     write(*,*) 'Usage: ',trim(arg),' input_file.'
     stop
  endif

  ! newunit specifier is instrodueced in Fortran 2008.
  open(newunit=unit,file=trim(model_infile),status='old')
  read(unit,parameters)
  close(unit)

  par%nprint = par%no_print
  if (par%nprint >= par%no_photons) par%nprint = par%no_photons/10
  if (par%no_photons < 10) par%nprint = 1

  par%nphotons   = par%no_photons
  par%nscatt_tot = 0.0_wp
  if (par%nx == 1 .or. par%ny == 1 .or. par%nz == 1) par%xyz_symmetry = .false.

  ! if rmax > 0, then the system is regarded as a sphere.
  if (par%rmax > 0.0_wp) then
     par%xmax = par%rmax
     par%ymax = par%rmax
     par%zmax = par%rmax
     if (par%nr > 1) then
        par%nx = par%nr
        par%ny = par%nr
        par%nz = par%nr
     endif
     par%nx   = maxval([par%nx,par%ny,par%nz])
     par%ny   = par%nx
     par%nz   = par%nx
  endif

  !--- MPI-related parameters
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpar%nproc, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpar%p_rank, ierr)
  !The third argument, key, determines the ordering (rank) within each new communicator.
  !call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, mpar%hostcomm, ierr)
  call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, mpar%p_rank, MPI_INFO_NULL, mpar%hostcomm, ierr)
  call MPI_COMM_RANK(mpar%hostcomm, mpar%h_rank, ierr)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, mpar%h_rank, mpar%p_rank, mpar%SAME_HRANK_COMM, ierr)
  call MPI_COMM_SIZE(mpar%SAME_HRANK_COMM, mpar%SAME_HRANK_NPROC, ierr)
  if (par%use_shared_memory) then
     par%use_shared_memory = .false.
     write(*,'(a)') 'par%use_shared_memory = .true. is not allowed.'
  endif

  if (len_trim(par%scatt_mat_file) == 0) par%use_stokes = .false.
  if (par%use_stokes) then
     call setup_scattering_matrix(par%scatt_mat_file)
  endif

  select case(trim(par%distance_unit))
     case ('kpc')
        par%distance2cm = kpc2cm
     case ('pc')
        par%distance2cm = pc2cm
     case ('au')
        par%distance2cm = au2cm
     case ('')
        par%distance2cm = 1.0_wp
     case default
        par%distance2cm = kpc2cm
  end select

  if (len_trim(par%out_file) == 0) then
     par%base_name = trim(get_base_input_name(model_infile))
     par%out_file  = trim(par%base_name)//'_obs.fits.gz'
  endif

  if (mpar%p_rank == 0) then
     write(*,'(a)')        ''
     write(*,'(3a)')       '+++++ ',trim(model_infile),' +++++'
     write(*,'(a,L1)')     'Use_Shared_Memory         : ', par%use_shared_memory
     write(*,'(a,L1)')     'Use_Master_Slave          : ', par%use_master_slave
     write(*,'(a,2f7.4)')  'Dust Parameters(a, g)     : ', par%albedo, par%hgg
     write(*,'(a,es12.3)') 'Dust Extinction per H     : ', par%cext_dust
     write(6,'(a,i14)')    'Total photons             : ', par%nphotons
  endif

  return
  end subroutine read_input
  !---------------------------------------
  subroutine setup_scattering_matrix(scatt_mat_file)
  use define
  use mathlib
  use memory_mod
  use random
  use mpi
  implicit none
  character(len=*), intent(in) :: scatt_mat_file

  character(len=128) :: string_tmp
  real(kind=wp) :: lambda,cext,albedo,hgg
  real(kind=wp) :: S11_norm
  integer       :: i, unit, ierr

  !--- setup Mueller Matrix
  !--- In this code, we will assume spherical dust grains.
  !--- For non-spherical grains, we also need to consider componets (i.e., S22, S21, S13, etc)
  if (len_trim(scatt_mat_file) > 0) then
     if (mpar%h_rank == 0) then
        open(newunit=unit,file=trim(scatt_mat_file),status='old')
        read(unit,*) string_tmp
        read(unit,*) lambda,cext,albedo,hgg,scatt_mat%nPDF
        !read(unit,*) lambda,albedo,hgg,scatt_mat%nPDF
        par%albedo  = albedo
        par%hgg     = hgg
        ! par%lambda0 and lambda is in units of micron.
        par%lambda0 = lambda
        read(unit,*) string_tmp
     endif
     call MPI_BCAST(par%albedo,     1, MPI_DOUBLE_PRECISION, 0, mpar%hostcomm, ierr)
     call MPI_BCAST(par%hgg,        1, MPI_DOUBLE_PRECISION, 0, mpar%hostcomm, ierr)
     call MPI_BCAST(par%lambda0,    1, MPI_DOUBLE_PRECISION, 0, mpar%hostcomm, ierr)
     call MPI_BCAST(scatt_mat%nPDF, 1, MPI_INTEGER,          0, mpar%hostcomm, ierr)

     call create_shared_mem(scatt_mat%coss,      [scatt_mat%nPDF])
     call create_shared_mem(scatt_mat%S11,       [scatt_mat%nPDF])
     call create_shared_mem(scatt_mat%S12,       [scatt_mat%nPDF])
     call create_shared_mem(scatt_mat%S33,       [scatt_mat%nPDF])
     call create_shared_mem(scatt_mat%S34,       [scatt_mat%nPDF])
     call create_shared_mem(scatt_mat%alias,     [scatt_mat%nPDF-1])
     call create_shared_mem(scatt_mat%phase_PDF, [scatt_mat%nPDF-1])

     if (mpar%h_rank == 0) then
        do i=1,scatt_mat%nPDF
           read(unit,*) scatt_mat%coss(i),scatt_mat%S11(i),scatt_mat%S12(i),scatt_mat%S33(i),scatt_mat%S34(i)
        enddo
        close(unit)

        ! Normalize the Mueller matrixes
        S11_norm      = calc_Integral(scatt_mat%coss,scatt_mat%S11)
        scatt_mat%S11 = scatt_mat%S11/S11_norm
        scatt_mat%S12 = scatt_mat%S12/S11_norm
        scatt_mat%S33 = scatt_mat%S33/S11_norm
        scatt_mat%S34 = scatt_mat%S34/S11_norm

        ! setup alias (2021.08.30)
        do i=1,scatt_mat%nPDF-1
           scatt_mat%phase_PDF(i) = (scatt_mat%S11(i) + scatt_mat%S12(i))/2.0_wp
        enddo
        scatt_mat%phase_PDF = scatt_mat%phase_PDF / sum(scatt_mat%phase_PDF)
        call random_alias_setup(scatt_mat%phase_PDF, scatt_mat%alias)
     endif
  else
     scatt_mat%nPDF = 0
  endif
  end subroutine setup_scattering_matrix
  !---------------------------------------
  subroutine setup_procedure
  use define
  use raytrace
  use peelingoff_mod
  use scatter_mod
  use run_simulation_mod
  use write_mod
  implicit none

  !--- Initialize Random Number Generator
  call init_random_seed(par%iseed)

  !--- procedure pointer for raytrace routine
  raytrace_to_edge => raytrace_to_edge_car

  if (par%xyz_symmetry) then
     raytrace_to_tau => raytrace_to_tau_car_xyzsym
  else if (par%xy_periodic) then
     if (par%nx == 1 .and. par%ny == 1) then
        raytrace_to_tau => raytrace_to_tau_car_zonly
     else
        raytrace_to_tau => raytrace_to_tau_car_xyper
     endif
  else
     raytrace_to_tau => raytrace_to_tau_car
  endif

  !--- procedure pointer for peelingoff_mod
  if (par%use_stokes) then
    peeling_direct_photon    => peeling_direct_photon_nostokes
    peeling_scattered_photon => peeling_scattered_photon_stokes
  else
    peeling_direct_photon    => peeling_direct_photon_nostokes
    peeling_scattered_photon => peeling_scattered_photon_nostokes
  endif

  !--- procedure pointer for scattering routine
  if (par%use_stokes) then
     scattering => scatter_dust_stokes
  else
     scattering => scatter_dust_nostokes
  endif

  !--- procedure pointer for simulation run
  if (par%use_master_slave) then
     run_simulation => run_master_slave
  else
     run_simulation => run_equal_number
  endif

  !--- procedure pointer to write output
  write_output => write_output_car

  end subroutine setup_procedure
  !=================================================
end module setup_mod
