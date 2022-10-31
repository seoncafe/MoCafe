module run_simulation_mod
  use define
  use grid_mod
  use photon_mod
  use random
  use scatter_mod
  use omp_lib
contains
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
  do ip = 1, par%nphotons
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
