module run_simulation_mod
  use define
  use grid_mod
  use photon_mod
  use random
  use scatter_mod
#ifdef MPI
  use mpi
#endif
  use omp_lib
contains
  !============================================
  !--- This routine works now. But, further tests are required. (comment added 2020.11.02).
  subroutine run_master_slave(grid)
  implicit none
  type(grid_type), intent(inout) :: grid

  !--- local variables
  type(photon_type)   :: photon
  integer(kind=int64) :: ip, ii
  real(kind=wp) :: tau0, tau, wgt1, nscatt_tot
  real(kind=wp) :: dtime

  !--- for MPI & OpenMP
  integer :: ierr
  integer :: my_threadid, master, irank, ithread, ithread0
  integer :: ans, tag, worker
#ifdef MPI
  integer :: status(MPI_STATUS_SIZE)
  integer(kind=int64) :: numsent, numreceived

  !--- Initialize
  my_threadid = omp_get_thread_num()
  master      = 0
  nscatt_tot  = 0.0

  !--- master part
  if (mpar%rank == master .and. my_threadid == 0) then
     numsent = 0
     !--- note 0,1,2,...,num_nodes-1
     do irank = 0, mpar%num_nodes-1
        ithread0 = 0
        if (irank == 0) ithread0 = 1
        do ithread = ithread0, mpar%num_threads(irank+1)-1
           tag     = ithread
           numsent = numsent + par%num_send_at_once
           if (numsent > par%nphotons) numsent = -1
           !write(*,'(a,i6,a,2i4)') '-- Hello, sending ip = ', numsent, '   to (rank, thread) = ', irank, tag
           call MPI_SEND(numsent,1,MPI_INTEGER8,irank,tag,MPI_COMM_WORLD,ierr)
        enddo
     enddo
     do ip = 1, par%nphotons, par%num_send_at_once
        call MPI_RECV(ans,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
        worker = status(MPI_SOURCE)
        tag    = status(MPI_TAG)
        if (numsent < par%nphotons) then
           numsent = numsent + par%num_send_at_once
           !write(*,'(a,i6,a,2i4)') '-- Hello, sending ip = ', numsent, '   to (rank, thread) = ', worker, tag
           call MPI_SEND(numsent,1,MPI_INTEGER8,worker,tag,MPI_COMM_WORLD,ierr)
           if (mod(numsent, par%nprint) == 0) then
              call time_stamp(dtime)
              write(6,'(es14.5,a,f8.3,a)') dble(numsent),' photons calculated: ',dtime/60.0_wp,' mins'
           endif
        else
           call MPI_SEND(-1_int64,1,MPI_INTEGER8,worker,tag,MPI_COMM_WORLD,ierr)
        endif
     enddo
  !--- worker (slave) part
  else
     do while(.true.)
        call MPI_RECV(numreceived,1,MPI_INTEGER8,master,my_threadid,MPI_COMM_WORLD,status,ierr)
        if (numreceived == -1) exit
        do ii = 1, par%num_send_at_once
           ip = numreceived - par%num_send_at_once + ii
           !++++++++++++++++++++++++++++++++
           !--- Main Part of the simulation
           !--- Release photon
           call gen_photon(grid,photon)

           do while(photon%inside)
              !--- Find scattering location of tau
              if (photon%nscatt == 0) then
                 !--- Force photon to scatter at optical depth tau before edge of grid
                 call raytrace_to_edge(photon,grid,tau0)
                 wgt1       = 1.0_wp - exp(-tau0)
                 photon%wgt = photon%wgt * wgt1
                 tau        = -log(1.0_wp - rand_number()*wgt1)
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
           nscatt_tot = nscatt_tot + photon%nscatt
           !++++++++++++++++++++++++++++++++
        enddo
        ans = 1
        tag = my_threadid
        call MPI_SEND(ans,1,MPI_INTEGER,master,tag,MPI_COMM_WORLD,ierr)
     enddo
  endif
  !$OMP BARRIER
  !$OMP ATOMIC
  par%nscatt_tot = par%nscatt_tot + nscatt_tot
  !$OMP BARRIER
#endif
  end subroutine run_master_slave
  !============================================
  subroutine run_equal_number(grid)
  implicit none
  type(grid_type), intent(inout) :: grid

  !--- local variables
  type(photon_type)   :: photon
  integer(kind=int64) :: ip
  real(kind=wp) :: tau0, tau, wgt1
#ifdef __GFORTRAN__
  !--- reduction variable should be private in outer context? what about intel compiler?
  real(kind=wp), save :: nscatt_tot
#else
  real(kind=wp) :: nscatt_tot
#endif
  real(kind=wp) :: dtime

  !--- Photon loop
  nscatt_tot = 0.0
  !$OMP DO &
  !$OMP private(wgt1,tau0,tau,photon,dtime,ip) &
  !$OMP reduction(+:nscatt_tot) &
  !$OMP schedule(dynamic,1)
  do ip = mpar%rank + 1, par%nphotons, mpar%num_nodes
     !++++++++++++++++++++++++++++++++
     !--- Main Part of the simulation
     !--- Release photon
     call gen_photon(grid,photon)

     do while(photon%inside)
        !--- Find scattering location of tau
        if (photon%nscatt == 0) then
           !--- Force photon to scatter at optical depth tau before edge of grid
           call raytrace_to_edge(photon,grid,tau0)
           wgt1       = 1.0_wp - exp(-tau0)
           photon%wgt = photon%wgt * wgt1
           tau        = -log(1.0_wp - rand_number()*wgt1)
        else
           tau        = -log(rand_number())
        endif
        call raytrace_to_tau(photon,grid,tau)

        if (photon%inside) then
           call scattering(photon,grid)
        endif
     enddo
     !--- End of the main part of the simulation
     !++++++++++++++++++++++++++++++++
     nscatt_tot = nscatt_tot + photon%nscatt

     if (mod(ip, par%nprint) == 0 .or. ip == par%nphotons) then
        call time_stamp(dtime)
        write(6,'(es14.5,a,f8.3,a)') dble(ip),' photons calculated: ',dtime/60.0_wp,' mins'
     endif
  enddo
  !$OMP END DO
  !--- Note that reduced value is stored in MASTER thread.
  !$OMP MASTER
  par%nscatt_tot = nscatt_tot
  !$OMP END MASTER
  end subroutine run_equal_number
!--------------------------------------------------

end module run_simulation_mod
