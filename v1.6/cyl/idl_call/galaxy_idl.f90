subroutine galaxy_idl(no_photons,nprint4,hgg,albedo,luminosity,&
                      dust1,tauface1,dust_rscale1,dust_zscale1,dust_rmax1,dust_zmax1,&
                      dust2,tauface2,dust_rscale2,dust_zscale2,dust_rmax2,dust_zmax2,&
                      rmax,pmax,zmax,nr,np,nz,&
                      disk_name,disk_rscale,disk_zscale,disk_rmax,disk_zmax,&
                      bulge_name,sersic_index,Reff,axial_ratio,BulgeToDisk,&
                      inclination_angle,position_angle,phase_angle,distance,&
                      nxim,nyim,dxim,dyim,im_scatt,im_direc,im_scatt_sig,im_direc_sig,&
                      psf_file,left_right_fold)
  use define
  use random
  use grid_mod
  use setup_idl_mod
  use output_mod
  use scatter
  use photon_mod
  use fits_module
  use mathlib
  use raytrace, only: raytrace_to_edge_cyl, raytrace_to_tau_cyl
  use detect
  use omp_lib

!--- Parameters
  implicit none
! input/output variables
  real(kind=wp),      intent(in) :: no_photons,hgg,albedo,luminosity
  integer,            intent(in) :: nprint4
  character(len=128), intent(in) :: dust1,dust2,disk_name,bulge_name
  real(kind=wp),      intent(in) :: tauface1,dust_rscale1,dust_zscale1,dust_rmax1,dust_zmax1,&
                                    tauface2,dust_rscale2,dust_zscale2,dust_rmax2,dust_zmax2,&
                                    disk_rscale,disk_zscale,disk_rmax,disk_zmax,&
                                    sersic_index,Reff,axial_ratio,BulgeToDisk
  real(kind=wp),      intent(in) :: rmax,pmax,zmax
  integer,            intent(in) :: nr,np,nz,nxim,nyim
  real(kind=wp),      intent(in) :: inclination_angle,position_angle,phase_angle,distance,dxim,dyim
  real(kind=sp),   intent(inout) :: im_scatt(nxim,nyim),im_direc(nxim,nyim),&
                                    im_scatt_sig(nxim,nyim),im_direc_sig(nxim,nyim)
  character(len=128), intent(in) :: psf_file
  integer,            intent(in) :: left_right_fold

! local variables
  type(photon_type)   :: photon
  type(grid_type)     :: grid
  type(observer_type) :: observer
  type(output_type)   :: output, output_sum

  logical           :: inside_grid
  integer(kind=i8b) :: ip, my_nphotons
  real(kind=wp)     :: totscatt
  integer           :: nscatt
  real(kind=wp)     :: tau,tau0
  integer           :: i
  real(kind=dp)     :: time1,time2

!--- Initial set up
  time1 = omp_get_wtime()
  call setup_par(no_photons,nprint4,hgg,albedo,luminosity,&
                 dust1,tauface1,dust_rscale1,dust_zscale1,dust_rmax1,dust_zmax1,&
                 dust2,tauface2,dust_rscale2,dust_zscale2,dust_rmax2,dust_zmax2,&
                 rmax,pmax,zmax,nr,np,nz,&
                 disk_name,disk_rscale,disk_zscale,disk_rmax,disk_zmax,&
                 bulge_name,sersic_index,Reff,axial_ratio,BulgeToDisk,&
                 inclination_angle,position_angle,phase_angle,distance,&
                 nxim,nyim,dxim,dyim,&
                 psf_file,left_right_fold)
  call setup_idl(photon,grid,observer)
  totscatt = 0

  !$OMP PARALLEL &
  !$OMP default(shared) &
  !$OMP firstprivate(photon) &
  !$OMP private(my_nphotons,output)
  call output_create(output)
  call init_random_seed()
  my_nphotons  = 0
  par%nthreads = omp_get_num_threads()

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
           call peelingoff(photon,grid,observer,output)
           call scatt(photon)
           nscatt = nscatt + 1
        endif
     enddo

     totscatt = totscatt + nscatt
     if (par%nprint > 0 .and. (mod(ip,par%nprint) == 0 .or. ip == par%nphotons)) then
        time2 = omp_get_wtime()
        write(6,'(i11,a,f8.3,a)') ip,' packets completed: ',(time2-time1),' secs'
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
  im_scatt     = output_sum%scatt
  im_direc     = output_sum%direc
  im_scatt_sig = output_sum%scatt_sig
  im_direc_sig = output_sum%direc_sig

  call output_destroy(output_sum)
  time2 = omp_get_wtime()
  write(6,'(a,f7.2)')   'Average Number of scattering: ',dble(totscatt)/dble(par%nphotons)
  write(6,'(a,f8.3,a)') 'Total Excution Time: ', (time2-time1),' secs'

  return
end subroutine galaxy_idl
