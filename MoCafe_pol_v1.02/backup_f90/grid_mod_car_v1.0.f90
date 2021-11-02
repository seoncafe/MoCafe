module grid_mod
  use define
  use memory_mod
  !use read_mod
  implicit none
contains

  !-----------------
  subroutine observer_create()
    use mpi
    implicit none
    real(kind=wp) :: cosa,cosb,cosg
    real(kind=wp) :: sina,sinb,sing
    integer :: ierr

    !--- setup observer and image plane
    ! alpha = phase angle
    ! beta  = inclination angle
    ! gamma = posision angle
    !         Position angle is measured counter-clock wise from x-axis of the detector plane.
    ! Rmaxtrix is the inverse matrix of the rotation matrix from the observer frame to the grid frame.
    ! rotation sequence is (1) position_angle (gamma) about z-axis,
    !                      (2) inclination_angle (beta) about (new) x-axis, and
    !                      (3) alpha about (new) z-axis.
    cosa = cos(par%phase_angle*deg2rad)
    sina = sin(par%phase_angle*deg2rad)
    cosb = cos(par%inclination_angle*deg2rad)
    sinb = sin(par%inclination_angle*deg2rad)
    cosg = cos(par%position_angle*deg2rad)
    sing = sin(par%position_angle*deg2rad)

    observer%rmatrix(1,1) =  cosa*cosg - sina*cosb*sing
    observer%rmatrix(1,2) = -sina*cosg - cosa*cosb*sing
    observer%rmatrix(1,3) =  sinb*sing
    observer%rmatrix(2,1) =  cosa*sing + sina*cosb*cosg
    observer%rmatrix(2,2) = -sina*sing + cosa*cosb*cosg
    observer%rmatrix(2,3) = -sinb*cosg
    observer%rmatrix(3,1) =  sina*sinb
    observer%rmatrix(3,2) =  cosa*sinb
    observer%rmatrix(3,3) =  cosb

    if (par%distance <= 0.0_wp) then
       par%distance = maxval([par%xmax, par%ymax, par%zmax]) * 100.0_wp
    endif
    ! observer's coordinates represented in the grid system.
    observer%x = par%distance*sina*sinb
    observer%y = par%distance*cosa*sinb
    observer%z = par%distance*cosb

    ! observer's image plane
    observer%nxim = par%nxim
    observer%nyim = par%nyim
    if (par%dxim > 0.0_wp) then
       observer%dxim = par%dxim
    else
       observer%dxim = atan2(par%xmax,par%distance)/(observer%nxim/2.0_wp) * rad2deg
    endif
    if (par%dyim > 0.0_wp) then
       observer%dyim = par%dyim
    else
       observer%dyim = atan2(par%ymax,par%distance)/(observer%nyim/2.0_wp) * rad2deg
    endif
    observer%steradian_pix = observer%dxim * observer%dyim * (deg2rad**2)

    ! memory allocations of output images
    call create_mem(observer%scatt, [observer%nxim,observer%nyim], observer%w_scatt, mpar%h_rank, mpar%hostcomm, ierr)
    call create_mem(observer%direc, [observer%nxim,observer%nyim], observer%w_direc, mpar%h_rank, mpar%hostcomm, ierr)
    if (par%use_stokes) then
       call create_mem(observer%I,  [observer%nxim,observer%nyim], observer%w_I,     mpar%h_rank, mpar%hostcomm, ierr)
       call create_mem(observer%Q,  [observer%nxim,observer%nyim], observer%w_Q,     mpar%h_rank, mpar%hostcomm, ierr)
       call create_mem(observer%U,  [observer%nxim,observer%nyim], observer%w_U,     mpar%h_rank, mpar%hostcomm, ierr)
       call create_mem(observer%V,  [observer%nxim,observer%nyim], observer%w_V,     mpar%h_rank, mpar%hostcomm, ierr)
       observer%I = 0.0_wp
       observer%Q = 0.0_wp
       observer%U = 0.0_wp
       observer%V = 0.0_wp
    endif
    if (par%peel_tau) then
       call create_mem(observer%tau,[observer%nxim,observer%nyim], observer%w_tau,   mpar%h_rank, mpar%hostcomm, ierr)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    !--- let's do not consider xyz_symmetry condition. (2017-08-19)
    !--- later, the xyz_symmetry condition will be implemented.
    par%xyz_symmetry = .false.
  end subroutine observer_create
  !-----------------
  subroutine observer_destroy()
  use define
  use mpi
  implicit none
  integer :: ierr

  if (par%use_shared_memory) then
     call MPI_WIN_FENCE(0,observer%w_scatt, ierr)
     call MPI_WIN_FENCE(0,observer%w_direc, ierr)

     call MPI_WIN_FREE(observer%w_scatt,    ierr)
     call MPI_WIN_FREE(observer%w_direc,    ierr)
  else
     if (associated(observer%scatt)) deallocate(observer%scatt)
     if (associated(observer%direc)) deallocate(observer%direc)
  endif

  if (par%use_stokes) then
     if (par%use_shared_memory) then
        call MPI_WIN_FENCE(0,observer%w_I, ierr)
        call MPI_WIN_FENCE(0,observer%w_Q, ierr)
        call MPI_WIN_FENCE(0,observer%w_U, ierr)
        call MPI_WIN_FENCE(0,observer%w_V, ierr)

        call MPI_WIN_FREE(observer%w_I,    ierr)
        call MPI_WIN_FREE(observer%w_Q,    ierr)
        call MPI_WIN_FREE(observer%w_U,    ierr)
        call MPI_WIN_FREE(observer%w_V,    ierr)
     else
        if (associated(observer%I)) deallocate(observer%I)
        if (associated(observer%Q)) deallocate(observer%Q)
        if (associated(observer%U)) deallocate(observer%U)
        if (associated(observer%V)) deallocate(observer%V)
    endif
  endif
  end subroutine observer_destroy
  !---------------------------------------------------------------------------------
  subroutine grid_create(grid)
  use define
  use mpi
  use read_mod
  implicit none
  type(grid_type), intent(inout) :: grid

  ! local variables
  integer :: status, unit
  integer :: i,j,k
  integer :: nxcen,nycen,nzcen
  real(kind=wp) :: taueq,taupole
  real(kind=wp) :: dens
  real(kind=wp) :: xx,yy,zz,rr
  character(len=128) :: dens_out

  integer :: ierr

  ! external density, temperature, and velocity field files.
  if (len_trim(par%density_file) > 0) then
     !call get_dimension(trim(par%density_file),par%nx,par%ny,par%nz,par%xmax,par%ymax,par%zmax,status,&
     call get_dimension(trim(par%density_file),par%nx,par%ny,par%nz,status,&
                                               reduce_factor=par%reduce_factor)
     if (status /= 0 .and. mpar%p_rank == 0) then
        write(*,*) 'Density File ',trim(par%density_file),' does not have Dimensional Information!'
        write(*,*) 'You need to specify xmax, ymax, zmax, nx, ny, nz in input parameter file!'
     endif
  endif

  !--- grid's cell faces.
  grid%nx   = par%nx
  grid%ny   = par%ny
  grid%nz   = par%nz
  grid%xmax = par%xmax
  grid%ymax = par%ymax
  grid%zmax = par%zmax
  grid%rmax = par%rmax

  ! grid%i0, grid%j0, grid%k0 = cell indices when the photon met a boundary and refelected into the grid system.
  ! center of grid system is always located at (x,y,z) = (0,0,0).
  if (par%xyz_symmetry .or. par%z_symmetry) then
     if ((par%nz/2)*2 == par%nz) then
        grid%dz   = grid%zmax/grid%nz
        grid%zmin = 0.0_wp
        grid%k0   = 1
     else
        grid%dz   = grid%zmax/(grid%nz-0.5_wp)
        grid%zmin = -grid%dz/2.0_wp
        grid%k0   = 2
     endif
  else
     grid%dz   = 2.0_wp*par%zmax/par%nz
     grid%zmin = -grid%zmax
     grid%k0   = 0
  endif
  if (par%xyz_symmetry) then
     if ((par%nx/2)*2 == par%nx) then
        grid%dx   = grid%xmax/grid%nx
        grid%xmin = 0.0_wp
        grid%i0   = 1
     else
        grid%dx   = grid%xmax/(grid%nx-0.5_wp)
        grid%xmin = -grid%dx/2.0_wp
        grid%i0   = 2
     endif
     if ((par%ny/2)*2 == par%ny) then
        grid%dy   = grid%ymax/grid%ny
        grid%ymin = 0.0_wp
        grid%j0   = 1
     else
        grid%dy   = grid%ymax/(grid%ny-0.5_wp)
        grid%ymin = -grid%dy/2.0_wp
        grid%j0   = 2
     endif
  else
     grid%dx   = 2.0_wp*par%xmax/par%nx
     grid%dy   = 2.0_wp*par%ymax/par%ny
     grid%xmin = -grid%xmax
     grid%ymin = -grid%ymax
     grid%i0   = 0
     grid%j0   = 0
  endif
  grid%xrange = grid%xmax - grid%xmin
  grid%yrange = grid%ymax - grid%ymin
  grid%zrange = grid%zmax - grid%zmin

  !--- center indices.
  if (par%xyz_symmetry) then
     nxcen = 1
     nycen = 1
     nzcen = 1
  else
     nxcen = (grid%nx+1)/2
     nycen = (grid%ny+1)/2
     nzcen = (grid%nz+1)/2
  endif

  ! allocate shared memories for xface, yface, zface
  call create_shared_mem(grid%xface, [grid%nx+1], grid%w_xface, mpar%h_rank, mpar%hostcomm, ierr)
  call create_shared_mem(grid%yface, [grid%ny+1], grid%w_yface, mpar%h_rank, mpar%hostcomm, ierr)
  call create_shared_mem(grid%zface, [grid%nz+1], grid%w_zface, mpar%h_rank, mpar%hostcomm, ierr)
  if (mpar%h_rank == 0) then
     grid%xface(:) = [ ((i-1)*grid%dx + grid%xmin, i=1,grid%nx+1) ]
     grid%yface(:) = [ ((j-1)*grid%dy + grid%ymin, j=1,grid%ny+1) ]
     grid%zface(:) = [ ((k-1)*grid%dz + grid%zmin, k=1,grid%nz+1) ]
  endif

  ! allocate shared memories for rhokap
  call create_shared_mem(grid%rhokap, [grid%nx,grid%ny,grid%nz], grid%w_rhokap, mpar%h_rank, mpar%hostcomm, ierr)
  if (mpar%h_rank == 0) then
     grid%rhokap(:,:,:) = 0.0_wp
  endif
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  !---- Density Setup
  if (mpar%h_rank == 0) then
     if (len_trim(par%density_file) == 0) then
        !--- Note distance unit = 1 when no physical distance unit is provied. (2017-06-27)
        par%distance_unit  = ''
        par%distance2cm    = 1.0_wp
        grid%rhokap(:,:,:) = 1.0_wp
     else
        call read_3D(trim(par%density_file),grid%rhokap,reduce_factor=par%reduce_factor,centering=par%centering)
        do k=1,grid%nz
           grid%rhokap(:,:,k) = grid%rhokap(:,:,k) * par%cext_dust * par%distance2cm
        enddo
     endif
  endif

  ! Spherical Geometry?
  if (mpar%h_rank == 0) then
     if (par%rmax > 0.0_wp) then
       do k=1,grid%nz
         zz = (grid%zface(k)+grid%zface(k+1))/2.0_wp
         do j=1,grid%ny
           yy = (grid%yface(j)+grid%yface(j+1))/2.0_wp
           do i=1,grid%nx
             xx = (grid%xface(i)+grid%xface(i+1))/2.0_wp
             rr = sqrt(xx**2 + yy**2 + zz**2)
             if (rr > par%rmax) grid%rhokap(i,j,k) = 0.0_wp
           enddo
         enddo
       enddo
    endif
  endif

  ! Normalize density if density is not given in a physical unit (hydrogen number per cm^3).
  ! Note that normalization is done even when the density is not a constant.
  if (mpar%h_rank == 0) then
     if (par%taumax > 0.0_wp .and. len_trim(par%distance_unit) == 0) then
        dens = sum(grid%rhokap)/count(grid%rhokap > 0.0)
        do k=1,grid%nz
           grid%rhokap(:,:,k) = grid%rhokap(:,:,k)*(par%taumax/grid%zmax)/dens
        enddo
     endif
  endif
  !if (par%peel_tau .and. mpar%p_rank == 0) then
  !   dens_out = trim(par%base_name)//'_tau.fits.gz'
  !   call fits_open_new(unit,dens_out,status)
  !   call fits_append_image(unit,grid%rhokap,status)
  !   call fits_put_keyword(unit,'delz', grid%dz, 'delta z',status)
  !   call fits_close(unit,status)
  !endif
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  !--- Calculate equatorial and polar optical depths
  if (mpar%p_rank == 0) then
     if (par%xyz_symmetry) then
        taueq        = sum(grid%rhokap(:,nycen,nzcen))*grid%dx
        taupole      = sum(grid%rhokap(nxcen,nycen,:))*grid%dz
        if ((grid%nx/2)*2 /= grid%nx) taueq   = taueq   - grid%rhokap(1,nycen,nzcen)*(grid%dx/2.0_wp)
        if ((grid%nz/2)*2 /= grid%nz) taupole = taupole - grid%rhokap(nxcen,nycen,1)*(grid%dz/2.0_wp)
     else
        taueq        = sum(grid%rhokap(:,nycen,nzcen))*(grid%dx/2.0_wp)
        taupole      = sum(grid%rhokap(nxcen,nycen,:))*(grid%dz/2.0_wp)
     endif
  endif
  call MPI_BCAST(par%taumax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  !--- print out important parameters.
  if (mpar%p_rank == 0) then
     write(6,'(a,es16.3)') 'tau_equat (dust) : ',taueq
     write(6,'(a,es16.3)') 'tau_pole  (dust) : ',taupole
  endif
  end subroutine grid_create
  !------------------------------------
  subroutine grid_destroy(grid)
  use define
  use mpi
  implicit none
  type(grid_type), intent(inout) :: grid
  integer :: ierr

  call MPI_WIN_FENCE(0,grid%w_xface,  ierr)
  call MPI_WIN_FENCE(0,grid%w_yface,  ierr)
  call MPI_WIN_FENCE(0,grid%w_zface,  ierr)
  call MPI_WIN_FENCE(0,grid%w_rhokap, ierr)

  call MPI_WIN_FREE(grid%w_xface,   ierr)
  call MPI_WIN_FREE(grid%w_yface,   ierr)
  call MPI_WIN_FREE(grid%w_zface,   ierr)
  call MPI_WIN_FREE(grid%w_rhokap,  ierr)
  end subroutine grid_destroy
  !=================================================
  function get_base_name(filename) result(base_name)
  character(len=*), intent(in) :: filename
  character(len=128) :: base_name
  integer :: i1,i2
  i1        = index(filename,'/', back=.true.)+1
  i2        = index(filename,'.fits')-1
  base_name = trim(filename(i1:i2))
  return
  end function get_base_name
end module grid_mod
