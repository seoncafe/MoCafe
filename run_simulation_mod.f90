module run_simulation_mod
  use define
  use grid_mod
  use photon_mod
  use random
  use scan_mod, only : scan_tau_step
  use mpi
contains
  subroutine run_master_slave(grid)
  implicit none
  type(grid_type), intent(inout) :: grid

  !--- local variables
  type(photon_type)   :: photon
  integer(kind=int64) :: ip, ii
  real(kind=wp) :: tau, tau0, wgt1
  real(kind=wp) :: dtime

  !--- for MPI
  integer :: ierr
  integer :: master, worker, irank, ans, tag
  integer :: status(MPI_STATUS_SIZE)
  integer(kind=int64) :: numsent, numreceived, numdone

  !--- Photon loop
  master   = 0
  if (mpar%p_rank == master) then
     numsent = 0_int64
     ! note 1,2,...,nproc-1
     do irank = 1, min(par%nphotons,int(mpar%nproc-1,int64))
        tag     = 1
        numsent = numsent + par%num_send_at_once
        call MPI_SEND(numsent,1,MPI_INTEGER8,irank,tag,MPI_COMM_WORLD,ierr)
     enddo
     !--- Send immediate termination to workers that received no initial batch
     !    (prevents deadlock when nphotons < nproc-1).
     do irank = int(min(par%nphotons,int(mpar%nproc-1,int64))) + 1, mpar%nproc - 1
        tag = 0
        call MPI_SEND(MPI_BOTTOM,0,MPI_INTEGER8,irank,tag,MPI_COMM_WORLD,ierr)
     enddo
     numdone = 0_int64
     do ip = 1, par%nphotons, par%num_send_at_once
        call MPI_RECV(ans,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
        worker  = status(MPI_SOURCE)
        !--- cap printed counter at nphotons (the final chunk may be partial)
        numdone = min(numdone + par%num_send_at_once, par%nphotons)
        if (mod(numdone,int(par%nprint,int64)) == 0) then
           call time_stamp(dtime)
           write(6,'(es14.5,a,f8.3,a)') dble(numdone),' photons calculated: ',dtime/60.0_wp,' mins'
        endif
        if (numsent < par%nphotons) then
           tag     = 1
           numsent = numsent + par%num_send_at_once
           call MPI_SEND(numsent,1,MPI_INTEGER8,worker,tag,MPI_COMM_WORLD,ierr)
        else
           tag = 0
           call MPI_SEND(MPI_BOTTOM,0,MPI_INTEGER8,worker,tag,MPI_COMM_WORLD,ierr)
        endif
     enddo
  else
     do while(.true.)
        call MPI_RECV(numreceived,1,MPI_INTEGER8,master,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
        if (status(MPI_TAG) == 0) exit
        do ii = 1, par%num_send_at_once
           ip = numreceived - par%num_send_at_once + ii
           !--- skip indices that exceed the requested total (final chunk may
           !    be partial when nphotons is not a multiple of num_send_at_once)
           if (ip > par%nphotons) exit
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
                 !--- tau0 = 0 & wgt1 = 0 occurs when photons are generated at
                 !--- locations where rho = 0 (e.g. external illumination through
                 !--- a void). Force the photon out of the grid in that case.
                 if (tau0 > 0.0_wp) then
                    tau = -log(1.0_wp - rand_number()*wgt1)
                 else
                    tau = hugest
                 endif
              else
                 tau = -log(rand_number())
              endif
              call raytrace_to_tau(photon,grid,tau)

              if (photon%inside) then
                 !--- hand the free-path optical depth of this segment to the
                 !--- tau scan (no-op when par%use_tau_list = .false.).
                 if (par%use_tau_list) scan_tau_step = tau
                 call scattering(photon,grid)
              endif
           enddo
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
  real(kind=wp) :: tau, tau0, wgt1
  real(kind=wp) :: dtime

  !--- for MPI
  integer :: ierr

  !--- Photon loop
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
           wgt1       = 1.0_wp - exp(-tau0)
           photon%wgt = photon%wgt * wgt1
           !--- tau0 = 0 protection (see run_master_slave for rationale).
           if (tau0 > 0.0_wp) then
              tau = -log(1.0_wp - rand_number()*wgt1)
           else
              tau = hugest
           endif
        else
           tau = -log(rand_number())
        endif
        call raytrace_to_tau(photon,grid,tau)

        if (photon%inside) then
           !--- hand the free-path optical depth of this segment to the
           !--- tau scan (no-op when par%use_tau_list = .false.).
           if (par%use_tau_list) scan_tau_step = tau
           call scattering(photon,grid)
        endif
     enddo
     !--- End of the main part of the simulation
     !++++++++++++++++++++++++++++++++

     if (mod(ip,int(par%nprint,int64)) == 0 .or. ip == par%nphotons) then
        call time_stamp(dtime)
        write(6,'(es14.5,a,f8.3,a)') dble(ip),' photons calculated: ',dtime/60.0_wp,' mins'
     endif
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
end subroutine run_equal_number

end module run_simulation_mod
