subroutine galaxy_car_idl(no_photons,nprint,hgg_in,albedo_in,luminosity,iseed4,&
                          dust1,tauface1,dust_rscale1,dust_zscale1,dust_rmax1,dust_zmax1,&
                          dust2,tauface2,dust_rscale2,dust_zscale2,dust_rmax2,dust_zmax2,&
                          xmax,ymax,zmax,nx,ny,nz,&
                          source_name,source_rscale,source_zscale,source_rmax,source_zmax,&
                          obs_angles,obs_distance,nxim,nyim,dxim,dyim,im_scatt,im_direc)

  use define
  use random
  use scatter
#ifdef _OPENMP
  use omp_lib
#endif

!--- Parameters
  implicit none
! input/output variables
  real(kind=wp), intent(in) :: no_photons,hgg_in,albedo_in,luminosity
  integer, intent(in) :: nprint
  character(len=100), intent(in) :: dust1,dust2,source_name
  real(kind=wp), intent(in) :: tauface1,dust_rscale1,dust_zscale1,dust_rmax1,dust_zmax1,&
                               tauface2,dust_rscale2,dust_zscale2,dust_rmax2,dust_zmax2,&
                               source_rscale,source_zscale,source_rmax,source_zmax
  real(kind=wp), intent(in) :: xmax,ymax,zmax
  integer, intent(in) :: nx,ny,nz,nxim,nyim
!  real(kind=wp), intent(inout) :: im_scatt(nxim,nyim),im_direc(nxim,nyim)
  real(kind=sp), intent(inout) :: im_scatt(nxim,nyim),im_direc(nxim,nyim)
  real(kind=wp), intent(in) :: obs_angles(3),obs_distance,dxim,dyim
  integer, intent(in) :: iseed4

! local variables
  logical :: inside_grid

  integer(kind=i8b) :: nphotons, ip
  real(kind=wp) :: totscatt
  integer :: nscatt

  type(photon_type)   :: photon
  type(grid_type)     :: grid
  type(dust_type)     :: dust
  type(source_type)   :: source
  type(observer_type) :: observer
  type(output_type)   :: output, output_sum

  real(kind=wp)     :: wgt,tau,tau0
  integer(kind=iwp) :: iseed,seed(4)
#ifndef _OPENMP
  integer, parameter :: nthreads = NTHREADS
  integer(kind=iwp) :: seed_array(4,nthreads)
  integer :: j
#endif

  integer :: i
#ifdef _OPENMP
  real(kind=dp) :: time1,time2
#else
  real(kind=wp) :: time1,time2
#endif

!--- Initial set up
#ifdef _OPENMP
  time1 = omp_get_wtime()
#else
  call cpu_time(time1)
#endif
  iseed = iseed4
  call setup_idl(no_photons,hgg_in,albedo_in,luminosity,&
                 dust1,tauface1,dust_rscale1,dust_zscale1,dust_rmax1,dust_zmax1,&
                 dust2,tauface2,dust_rscale2,dust_zscale2,dust_rmax2,dust_zmax2,&
                 xmax,ymax,zmax,nx,ny,nz,&
                 source_name,source_rscale,source_zscale,source_rmax,source_zmax,&
                 obs_angles,obs_distance,nxim,nyim,dxim,dyim,&
                 photon,grid,source,observer,output,nphotons)

  write(6,'(a,2i11)') 'Total photons: ',nphotons
  totscatt = 0

  !$OMP PARALLEL &
  !$OMP firstprivate(output,photon) &
  !$OMP private(seed) &
  !$OMP shared(nscatt_max,wgt_min,nphotons,iseed,grid,observer,dust,source,output_sum)
  call setup_pointers(output)
#ifdef _OPENMP
  !$OMP DO private(i) ordered
  do i=1,omp_get_num_threads()
     !$OMP ORDERED
     call set_kiss(iseed,seed)
     !write(6,'(a,i3,5i21)') 'I, ISEED, SEED: ',i,iseed,seed
     !$OMP END ORDERED
  enddo
  !$OMP END DO
#else
  do i=1,nthreads
     call set_kiss(iseed,seed)
     seed_array(:,i) = seed(:)
     write(6,'(a,i3,5i21)') 'I, ISEED, SEED: ',i,iseed,seed
  enddo
#endif

!--- Photon loop
  !$OMP DO &
  !$OMP private(wgt,tau0,tau,nscatt,inside_grid) &
  !$OMP reduction(+:totscatt) &
  !$OMP schedule(static,1)
  do ip=1,nphotons
#ifndef _OPENMP
     j = mod(ip,nthreads)
     if (j == 0) j = nthreads
     seed(:) = seed_array(:,j)
#endif
     !--- Release photon
     call generate_photon(source,grid,photon,seed)

     !--- Add direct photon
     call add_direct(photon,grid,observer,output)

     !--- Find optical depth (tau0) to edge of grid
     call raytrace_to_edge_car(photon,grid,tau0)

     if (tau0 > eps_sp) then
        wgt = (1.0_wp - exp(-tau0))

        !--- Force photon to scatter at optical depth tau before edge of grid
        tau = -log(1.0_wp - random_kiss(seed)*wgt)
        inside_grid = .true.

        !--- Find scattering location of tau
        call raytrace_to_tau_car(photon,grid,tau,inside_grid)
        nscatt = 1

        do while (inside_grid .and. (nscatt <= nscatt_max))
        !do while (inside_grid .and. (wgt >= wgt_min))
           totscatt = totscatt + nscatt*wgt
           wgt = wgt * albedo
           call peelingoff(photon,grid,observer,output,wgt)
           call scatt(photon%vx,photon%vy,photon%vz,seed)
           tau = -log(random_kiss(seed))
           call raytrace_to_tau_car(photon,grid,tau,inside_grid)
           if (inside_grid) nscatt = nscatt + 1
        enddo
        !totscatt = totscatt + nscatt*wgt
     endif
     if (nprint > 0 .and. (mod(ip,nprint) == 0 .or. ip == nphotons)) then
#ifdef _OPENMP
        time2 = omp_get_wtime()
#else
        call cpu_time(time2)
#endif
        write(6,'(i11,a,f8.3,a)') ip,' photons completed: ',(time2-time1),' secs'
     endif
#ifndef _OPENMP
     seed_array(:,j) = seed(:)
#endif
  enddo
  !$OMP END DO

  !$OMP SINGLE
  output_sum = output
  call setup_pointers(output_sum)
  !$OMP END SINGLE
  !$OMP CRITICAL
  output_sum%scatt(:,:) = output_sum%scatt(:,:) + output%scatt(:,:)
  output_sum%direc(:,:) = output_sum%direc(:,:) + output%direc(:,:)
  !$OMP END CRITICAL

  !$OMP END PARALLEL
  write(6,'(a,f7.2)') 'Average Number of scattering: ',dble(totscatt)/dble(nphotons)
  im_scatt(:,:) = output_sum%scatt(:,:)
  im_direc(:,:) = output_sum%direc(:,:)

#ifdef _OPENMP
  time2 = omp_get_wtime()
#else
  call cpu_time(time2)
#endif
  write(6,'(a,f8.3,a)') 'Total Excution Time: ', (time2-time1),' secs'

  return
end subroutine galaxy_car_idl
