module jtally_mod
!--- Per-cell mean-intensity tally J_lambda(cell) for MoCafe v2.00 (Stage 2).
!--- Lucy (1999) pathlength estimator: every actual photon path segment
!--- contributes wgt*dl to its (wavelength-bin, cell) slot; the mean
!--- intensity follows as J_lambda = E_packet * Sum(wgt*dl) / (4 pi V dlam).
!---
!--- The forced first scattering biases the free-path distribution of the
!--- nscatt = 0 flight, so that flight is tallied ANALYTICALLY instead:
!--- raytrace_to_edge_car already walks the full ray to the grid edge, and
!--- (with jt_first = .true.) accumulates the exact expectation
!---    wgt * Int exp(-tau_lambda(l)) dl
!---      = wgt * (exp(-s*tau_in) - exp(-s*tau_out)) / (rhokap*s)
!--- per cell (zero-variance direct component).  Flights with nscatt >= 1
!--- sample the standard exponential free path and are tallied along the
!--- actual walked segments in raytrace_to_tau_car (unbiased).
!---
!--- The tally array is per-MPI-rank private and MPI_REDUCEd once at the
!--- end of the run.  Stage 2 scope: Cartesian grid, SED mode only.
  use define
  implicit none
  public

  logical :: jt_on    = .false.  ! master switch (par%save_jlam, set in jtally_setup)
  logical :: jt_first = .false.  ! analytic first-flight tally inside raytrace_to_edge_car
  real(kind=wp), pointer :: jt_sum(:,:,:,:) => null()  ! (nlambda, nx, ny, nz): Sum(wgt*dl) [code length]
  real(kind=wp) :: jt_eabs = 0.0_wp  ! independent absorbed-energy counter: Sum wgt*(1-albedo) at scatterings

contains
  !---------------------------------------------------------------
  subroutine jtally_setup(grid)
  use sed_mod,    only : sed_nlam
  use memory_mod, only : create_mem
  implicit none
  type(grid_type), intent(in) :: grid
  real(kind=wp) :: mem_gb

  mem_gb = real(sed_nlam,wp)*grid%nx*grid%ny*grid%nz*8.0_wp/1024.0_wp**3
  if (mpar%p_rank == 0) then
     write(*,'(a,f8.3,a)') 'J_lambda tally: ', mem_gb, ' GB per MPI rank'
     if (mem_gb > 4.0_wp) write(*,'(a)') &
        'WARNING: J_lambda tally is large; consider fewer wavelength bins or cells.'
  endif
  call create_mem(jt_sum, [sed_nlam, grid%nx, grid%ny, grid%nz])
  jt_sum(:,:,:,:) = 0.0_wp
  jt_eabs = 0.0_wp
  jt_on   = .true.
  end subroutine jtally_setup

  !---------------------------------------------------------------
  subroutine jtally_reduce()
  use mpi
  implicit none
  integer :: ierr, k
  !--- ALLREDUCE (not reduce-to-0) so every rank holds the full tally: the
  !--- dust-emission stage (dustemis_mod) distributes cells across all ranks
  !--- and reads jt_sum locally.  Reduce one wavelength slab at a time to keep
  !--- each MPI count within a 32-bit int for large (nlam x ncell) grids.
  if (.not. jt_on) return
  do k = 1, size(jt_sum, 4)
     !--- jt_sum(:,:,:,k) is contiguous (leading dims full), so no copy.
     call MPI_ALLREDUCE(MPI_IN_PLACE, jt_sum(:,:,:,k), int(size(jt_sum(:,:,:,k))), &
                        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE, jt_eabs, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  end subroutine jtally_reduce

  !---------------------------------------------------------------
  !--- convert the tally to J_lambda and write '<base>_jlam.<ext>'.
  !--- J_lambda units: [luminosity unit] / dist_cm^2 / um / sr, where
  !--- dist_cm = par%distance2cm (1 when no physical unit is given).
  !--- Also writes the wavelength-integrated J_bol(x,y,z) and, for the
  !--- energy-conservation check, prints the absorbed luminosity from the
  !--- tally (A) against the independent event counter (B).
  subroutine jtally_write(grid)
  use sed_mod,    only : sed_nlam, sed_wave, sed_dwave, sed_sext, sed_albedo, sed_cext_ref
  use iofile_mod
  use utility,    only : get_base_name
  implicit none
  type(grid_type), intent(in) :: grid
  type(io_file_type) :: file
  character(len=192) :: filename
  real(kind=wp), allocatable :: jbol(:,:,:)
  real(kind=wp) :: E_p, vol, fac, eabs_A, eabs_B
  integer :: status, il, i, j, k

  if (.not. jt_on .or. mpar%p_rank /= 0) return

  E_p = par%luminosity/dble(par%nphotons)
  vol = grid%dx*grid%dy*grid%dz

  !--- energy-conservation check:
  !--- A = absorbed luminosity from the J tally
  !---   = E_p * Sum_{il,cell} jt_sum * rhokap(cell) * s_ext(il) * (1-albedo(il))
  !--- B = E_p * Sum wgt*(1-albedo) over scattering events (independent).
  eabs_A = 0.0_wp
  do k = 1, grid%nz
  do j = 1, grid%ny
  do i = 1, grid%nx
     if (grid%rhokap(i,j,k) > 0.0_wp) then
        eabs_A = eabs_A + grid%rhokap(i,j,k) * &
                 sum(jt_sum(:,i,j,k)*sed_sext(:)*(1.0_wp - sed_albedo(:)))
     endif
  enddo
  enddo
  enddo
  eabs_A = eabs_A * E_p
  eabs_B = jt_eabs * E_p
  write(*,'(a)')        '--- J_lambda tally: energy conservation check ---'
  write(*,'(a,es14.6)') 'absorbed L (pathlength tally, A): ', eabs_A
  write(*,'(a,es14.6)') 'absorbed L (event counter,    B): ', eabs_B
  if (eabs_B > 0.0_wp) write(*,'(a,f10.6)') 'ratio A/B                       : ', eabs_A/eabs_B

  !--- convert Sum(wgt*dl) -> J_lambda.
  fac = E_p / (fourpi*vol*par%distance2cm**2)
  do il = 1, sed_nlam
     jt_sum(il,:,:,:) = jt_sum(il,:,:,:) * (fac/sed_dwave(il))
  enddo
  allocate(jbol(grid%nx, grid%ny, grid%nz))
  do k = 1, grid%nz
  do j = 1, grid%ny
  do i = 1, grid%nx
     jbol(i,j,k) = sum(jt_sum(:,i,j,k)*sed_dwave(:))
  enddo
  enddo
  enddo

  status = 0
  filename = trim(get_base_name(par%out_file))//'_jlam'//trim(io_file_extension(par%file_format))
  call io_open_new(file, trim(filename), status)
  call io_append_image(file, jt_sum, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','J_lambda','J(lambda,x,y,z) mean intensity',status)
  call io_put_keyword(file,'J_UNIT','luminosity/dist_cm^2/um/sr','J_lambda unit',status)
  call io_put_keyword(file,'SED_NLAM', sed_nlam,       'number of wavelength bins',     status)
  call io_put_keyword(file,'SED_LREF', par%lambda_ref, 'reference wavelength [um]',     status)
  call io_put_keyword(file,'SED_CREF', sed_cext_ref,   'C_ext/H at lambda_ref [cm^2/H]',status)
  call io_put_keyword(file,'TOT_LUM',  par%luminosity, 'total luminosity',              status)
  call io_put_keyword(file,'DIST_CM',  par%distance2cm,'distance unit (cm)',            status)
  call io_put_keyword(file,'nphotons', par%no_photons, 'number of photons',             status)
  call io_put_keyword(file,'taumax',   par%taumax,     'tau_max at lambda_ref',         status)
  call io_put_keyword(file,'EABS_A',   eabs_A, 'absorbed L (pathlength tally)',         status)
  call io_put_keyword(file,'EABS_B',   eabs_B, 'absorbed L (event counter)',            status)
  call io_append_image(file, jbol, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','J_bol','wavelength-integrated J(x,y,z)',status)
  call io_put_keyword(file,'J_UNIT','luminosity/dist_cm^2/sr','J_bol unit',status)
  call io_append_image(file, sed_wave, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','Wavelength','bin centers [um]',status)
  call io_append_image(file, sed_dwave, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','Dwavelength','bin widths [um]',status)
  call io_close(file, status)
  write(*,'(2a)') 'J_lambda written to: ', trim(filename)

  deallocate(jbol)
  end subroutine jtally_write

end module jtally_mod
