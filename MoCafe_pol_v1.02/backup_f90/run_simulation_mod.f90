module run_simulation_mod
  use define
  use utility
  use mathlib
  use grid_mod
  use photon_mod
  use random
  use mpi
contains
  subroutine run_master_slave(grid)
  implicit none
  type(grid_type), intent(inout) :: grid

  !--- local variables
  type(photon_type)   :: photon
  integer(kind=int64) :: ip
  real(kind=wp) :: tau, tau0
  real(kind=wp) :: time1,time2

  !--- for MPI
  integer :: ierr
  integer :: master, irank
  integer :: status(MPI_STATUS_SIZE)
  integer :: ans, tag
  integer :: worker
  integer(kind=int64) :: numsent

  !--- Photon loop
  time1    = MPI_WTIME()
  master   = 0
  if (mpar%p_rank == master) then
     numsent = 0
     ! note 1,2,...,nproc-1
     do irank = 1, min(par%nphotons,mpar%nproc-1)
        tag     = 1
        numsent = numsent+1
        call MPI_SEND(numsent,1,MPI_INTEGER8,irank,tag,MPI_COMM_WORLD,ierr)
     enddo
     do ip = 1, par%nphotons
        call MPI_RECV(ans,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
        worker  = status(MPI_SOURCE)
        if (numsent < par%nphotons) then
           tag     = 1
           numsent = numsent+1
           call MPI_SEND(numsent,1,MPI_INTEGER8,worker,tag,MPI_COMM_WORLD,ierr)
        else
           tag = 0
           call MPI_SEND(MPI_BOTTOM,0,MPI_INTEGER8,worker,tag,MPI_COMM_WORLD,ierr)
        endif
        if (mod(ip,int(par%nprint,int64)) == 0) then
           time2 = MPI_WTIME()
           write(6,'(i14,a,f8.3,a)') ip,' photons calculated: ',(time2-time1)/60.0_wp,' mins'
        endif
     enddo
  else
     do while(.true.)
        call MPI_RECV(ip,1,MPI_INTEGER8,master,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
        if (status(MPI_TAG) == 0) exit
        !++++++++++++++++++++++++++++++++
        !--- Main Part of the simulation
        !--- Release photon
        call gen_photon(grid,photon)

        do while(photon%inside)
           !--- Find scattering location of tau
           if (photon%nscatt == 0) then
              !--- Force photon to scatter at optical depth tau before edge of grid
              call raytrace_to_edge(photon,grid,tau0)
              photon%wgt = (1.0_wp - exp(-tau0))*photon%wgt
              tau = -log(1.0_wp - rand_number()*photon%wgt)
           else
              tau = -log(rand_number())
           endif
           call raytrace_to_tau(photon,grid,tau)

           if (photon%inside) then
              call scattering(photon,grid)
           endif
        enddo
        !--- End of the main part of the simulation
        !++++++++++++++++++++++++++++++++
        ans = 1
        tag = 1
        call MPI_SEND(ans,1,MPI_INTEGER,master,tag,MPI_COMM_WORLD,ierr)
     enddo
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end subroutine run_master_slave

  !============================================
  subroutine run_equal_number(grid)
  implicit none
  type(grid_type), intent(inout) :: grid
  !--- local variables
  type(photon_type)   :: photon
  integer(kind=int64) :: ip
  real(kind=wp) :: tau, tau0
  real(kind=wp) :: time1,time2

  !--- for MPI
  integer :: ierr

  !--- Photon loop
  time1 = MPI_WTIME()
  do ip=mpar%p_rank+1,par%nphotons,mpar%nproc
     !++++++++++++++++++++++++++++++++
     !--- Main Part of the simulation
     !--- Release photon
     call gen_photon(grid,photon)

     do while(photon%inside)
        !--- Find scattering location of tau
        if (photon%nscatt == 0) then
           !--- Force photon to scatter at optical depth tau before edge of grid
           call raytrace_to_edge(photon,grid,tau0)
           photon%wgt = (1.0_wp - exp(-tau0))*photon%wgt
           tau = -log(1.0_wp - rand_number()*photon%wgt)
        else
           tau = -log(rand_number())
        endif
        call raytrace_to_tau(photon,grid,tau)

        if (photon%inside) then
           call scattering(photon,grid)
        endif
     enddo
     !--- End of the main part of the simulation
     !++++++++++++++++++++++++++++++++

     if (mod(ip,par%nprint) == 0 .or. ip == par%nphotons) then
        time2 = MPI_WTIME()
        write(6,'(es14.5,a,f8.3,a)') dble(ip),' photons calculated: ',(time2-time1)/60.0_wp,' mins'
     endif
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
end subroutine run_equal_number

end module run_simulation_mod
