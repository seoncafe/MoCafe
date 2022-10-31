module grid_mod
  use define
  use utility
  use memory_mod
  !use read_mod
  implicit none
contains

  !---------------------------------------------------------------------------------
  subroutine grid_create(grid)
  use define
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
     if (status /= 0) then
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

  ! allocate memories for xface, yface, zface
  call create_mem(grid%xface, [grid%nx+1])
  call create_mem(grid%yface, [grid%ny+1])
  call create_mem(grid%zface, [grid%nz+1])
  grid%xface(:) = [ ((i-1)*grid%dx + grid%xmin, i=1,grid%nx+1) ]
  grid%yface(:) = [ ((j-1)*grid%dy + grid%ymin, j=1,grid%ny+1) ]
  grid%zface(:) = [ ((k-1)*grid%dz + grid%zmin, k=1,grid%nz+1) ]

  ! allocate memories for rhokap
  call create_mem(grid%rhokap, [grid%nx,grid%ny,grid%nz])

  !---- Density Setup
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

  ! Spherical Geometry?
  do k=1,grid%nz
    zz = (grid%zface(k)+grid%zface(k+1))/2.0_wp
    do j=1,grid%ny
      yy = (grid%yface(j)+grid%yface(j+1))/2.0_wp
      do i=1,grid%nx
        xx = (grid%xface(i)+grid%xface(i+1))/2.0_wp
        rr = sqrt(xx**2 + yy**2 + zz**2)
        if (is_finite(par%density_rscale))      grid%rhokap(i,j,k) = exp(-rr / par%density_rscale)
        if (par%rmax > 0.0 .and. rr > par%rmax) grid%rhokap(i,j,k) = 0.0_wp
      enddo
    enddo
  enddo

  ! Normalize density if density is not given in a physical unit (hydrogen number per cm^3).
  ! Note that normalization is done even when the density is not a constant.
  if (par%taumax > 0.0_wp .and. len_trim(par%distance_unit) == 0) then
     dens = (2.0_wp*par%taumax)/sum(grid%rhokap(nxcen,nycen,:)*grid%dz)
     grid%rhokap(:,:,:) = grid%rhokap(:,:,:) * dens
  else if (par%tauhomo > 0.0_wp .and. len_trim(par%distance_unit) == 0) then
     dens = sum(grid%rhokap)/count(grid%rhokap > 0.0)
     do k=1,grid%nz
        grid%rhokap(:,:,k) = grid%rhokap(:,:,k)*(par%tauhomo/grid%zmax)/dens
     enddo
  endif

  !--- Calculate equatorial and polar optical depths
  if (par%xyz_symmetry) then
     taueq        = sum(grid%rhokap(:,nycen,nzcen))*grid%dx
     taupole      = sum(grid%rhokap(nxcen,nycen,:))*grid%dz
     if ((grid%nx/2)*2 /= grid%nx) taueq   = taueq   - grid%rhokap(1,nycen,nzcen)*(grid%dx/2.0_wp)
     if ((grid%nz/2)*2 /= grid%nz) taupole = taupole - grid%rhokap(nxcen,nycen,1)*(grid%dz/2.0_wp)
  else
     taueq        = sum(grid%rhokap(:,nycen,nzcen))*(grid%dx/2.0_wp)
     taupole      = sum(grid%rhokap(nxcen,nycen,:))*(grid%dz/2.0_wp)
  endif

  !--- print out important parameters.
  write(6,'(a,es16.3)') 'tau_equat (dust) : ',taueq
  write(6,'(a,es16.3)') 'tau_pole  (dust) : ',taupole
  end subroutine grid_create
  !------------------------------------
  subroutine grid_destroy(grid)
  use define
  implicit none
  type(grid_type), intent(inout) :: grid
  integer :: ierr

  if (associated(grid%rhokap)) deallocate(grid%rhokap)
  end subroutine grid_destroy
  !=================================================
end module grid_mod