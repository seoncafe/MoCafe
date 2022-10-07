module photon_mod
contains
  subroutine gen_photon(grid,photon)

  use define
  use random
  use peelingoff_mod
  implicit none

  type(grid_type),   intent(inout) :: grid
  type(photon_type), intent(out)   :: photon

  ! local variables
  real(kind=wp) :: sint,cost,phi,sinp,cosp,rp
  real(kind=wp) :: kx,ky,kz,mx,my,mz,nx,ny,nz

  !=== set up photon's position vector.
  if (trim(par%source_geometry) == 'uniform' .and. par%rmax > 0.0_wp) then
     rp   = (rand_number())**(1.0_wp/3.0_wp) * par%rmax
     cost = 2.0_wp*rand_number()-1.0_wp
     sint = sqrt(1.0_wp-cost*cost)
     phi  = twopi*rand_number()
     photon%x = rp*sint*cos(phi)
     photon%y = rp*sint*sin(phi)
     photon%z = rp*cost
  else if (trim(par%source_geometry) == 'uniform' .and. par%rmax <= 0.0_wp) then
     photon%x = (2.0_wp*rand_number()-1.0_wp)*grid%xmax
     photon%y = (2.0_wp*rand_number()-1.0_wp)*grid%ymax
     photon%z = (2.0_wp*rand_number()-1.0_wp)*grid%zmax
  else if (trim(par%source_geometry) == 'uniform_xy' .and. par%rmax <= 0.0_wp) then
     photon%x = grid%xrange*rand_number()+grid%xmin
     photon%y = grid%yrange*rand_number()+grid%ymin
     photon%z = 0.0_wp
  else if (trim(par%source_geometry) == 'gaussian') then
     photon%x = grid%xrange*rand_number()+grid%xmin
     photon%y = grid%yrange*rand_number()+grid%ymin
     photon%z = par%source_zscale/sqrt(2.0_wp)*rand_gauss()
  else if (trim(par%source_geometry) == 'exponential') then
     photon%x = grid%xrange*rand_number()+grid%xmin
     photon%y = grid%yrange*rand_number()+grid%ymin
     photon%z = par%source_zscale*rand_zexp(par%zmax/par%source_zscale)
  else if (trim(par%source_geometry(1:12)) == 'external_sph') then
     !--- random number for incident location
     cost = 2.0_wp*rand_number() - 1.0_wp
     sint = sqrt(1.0_wp-cost*cost)
     phi  = twopi*rand_number()
     cosp = cos(phi)
     sinp = sin(phi)

     !--- reference vectors for isotropically incident ray
     kx = sint*cosp
     ky = sint*sinp
     kz = cost
     mx = cost*cosp
     my = cost*sinp
     mz = -sint
     nx = -sinp
     ny =  cosp
     nz =  0.0_wp

     !--- position of the incident photon
     !--- please note that par%rmax should be positive
     photon%x = par%rmax*kx
     photon%y = par%rmax*ky
     photon%z = par%rmax*kz

     !--- random direction for propagation direction
     !--- note that pi/2 < theta < pi, i.e., -1 < cost < 0.
     cost = -sqrt(rand_number())
     sint = sqrt(1.0_wp-cost*cost)
     phi  = twopi*rand_number()
     cosp = cos(phi)
     sinp = sin(phi)

     photon%kx = sint * (cosp * mx + sinp * nx) + cost * kx
     photon%ky = sint * (cosp * my + sinp * ny) + cost * ky
     photon%kz = sint * (cosp * mz + sinp * nz) + cost * kz
  else if (trim(par%source_geometry(1:12)) == 'external_rec') then
     !--- random direction for propagation direction
     !--- note that pi/2 < theta < pi, i.e., -1 < cost < 0.
     cost = -sqrt(rand_number())
     sint = sqrt(1.0_wp-cost*cost)
     phi  = twopi*rand_number()
     cosp = cos(phi)
     sinp = sin(phi)

     !--- random number for incident location
     select case(floor(6.0*rand_number() + 1.0))
     case (1)
        photon%x = grid%xface(grid%nx+1)
        photon%y = grid%yrange * rand_number() + grid%ymin
        photon%z = grid%zrange * rand_number() + grid%zmin
        kx = 1.0
        ky = 0.0
        kz = 0.0
        mx = 0.0
        my = 1.0
        mz = 0.0
        nx = 0.0
        ny = 0.0
        nz = 1.0
     case (2)
        photon%x = grid%xface(1)
        photon%y = grid%yrange * rand_number() + grid%ymin
        photon%z = grid%zrange * rand_number() + grid%zmin
        kx = -1.0
        ky = 0.0
        kz = 0.0
        mx = 0.0
        my = 1.0
        mz = 0.0
        nx = 0.0
        ny = 0.0
        nz = 1.0
     case (3)
        photon%y = grid%yface(grid%ny+1)
        photon%z = grid%zrange * rand_number() + grid%zmin
        photon%x = grid%xrange * rand_number() + grid%xmin
        kx = 0.0
        ky = 1.0
        kz = 0.0
        mx = 1.0
        my = 0.0
        mz = 0.0
        nx = 0.0
        ny = 0.0
        nz = 1.0
     case (4)
        photon%y = grid%yface(1)
        photon%z = grid%zrange * rand_number() + grid%zmin
        photon%x = grid%xrange * rand_number() + grid%xmin
        kx = 0.0
        ky = -1.0
        kz = 0.0
        mx = 1.0
        my = 0.0
        mz = 0.0
        nx = 0.0
        ny = 0.0
        nz = 1.0
     case (5)
        photon%z = grid%zface(grid%nz+1)
        photon%x = grid%xrange * rand_number() + grid%xmin
        photon%y = grid%yrange * rand_number() + grid%ymin
        kx = 0.0
        ky = 0.0
        kz = 1.0
        mx = 1.0
        my = 0.0
        mz = 0.0
        nx = 0.0
        ny = 1.0
        nz = 0.0
     case (6)
        photon%z = grid%zface(1)
        photon%x = grid%xrange * rand_number() + grid%xmin
        photon%y = grid%yrange * rand_number() + grid%ymin
        kx = 0.0
        ky = 0.0
        kz = -1.0
        mx = 1.0
        my = 0.0
        mz = 0.0
        nx = 0.0
        ny = 1.0
        nz = 0.0
     end select

     photon%kx = sint * (cosp * mx + sinp * nx) + cost * kx
     photon%ky = sint * (cosp * my + sinp * ny) + cost * ky
     photon%kz = sint * (cosp * mz + sinp * nz) + cost * kz
  else
     photon%x = par%xs_point
     photon%y = par%ys_point
     photon%z = par%zs_point
  endif

  !--- Cell Index
  photon%icell = floor((photon%x-grid%xmin)/grid%dx)+1
  photon%jcell = floor((photon%y-grid%ymin)/grid%dy)+1
  photon%kcell = floor((photon%z-grid%zmin)/grid%dz)+1

  !=== set up photon's propagation and direction vetor.
  !--- isotropic emission
  if (trim(par%source_geometry(1:8)) /= 'external') then
     cost = 2.0_wp*rand_number()-1.0_wp
     sint = sqrt(1.0_wp-cost*cost)
     phi  = twopi*rand_number()
     cosp = cos(phi)
     sinp = sin(phi)

     !--- Set propagation direction cosines of the photon
     photon%kx = sint*cosp
     photon%ky = sint*sinp
     photon%kz = cost
  endif

  if (par%use_stokes) then
     if (trim(par%source_geometry(1:8)) == 'external') then
        !--- Find theta and phi in the fixed lab frame.
        cost = photon%kz
        sint = sqrt(1.0_wp - cost*cost)
        cosp = photon%kx/sint
        sinp = photon%ky/sint
     endif

     !--- Set the reference normal perpendicular to the propagation direction.
     photon%mx =  cost * cosp
     photon%my =  cost * sinp
     photon%mz = -sint
     photon%nx = -sinp
     photon%ny =  cosp
     photon%nz =  0.0_wp

     !--- Set the Stokes parameters (assume unpolarized light)
     photon%I = 1.0_wp
     photon%Q = 0.0_wp
     photon%U = 0.0_wp
     photon%V = 0.0_wp
  endif

  photon%nscatt = 0
  photon%wgt    = 1.0_wp
  photon%inside = .true.
  photon%albedo = par%albedo
  photon%hgg    = par%hgg

  !--- comment added, 2021.11.01
  !--- For external isotropic source, the direct ray is a sum of all rays coming from the sightline.
  !--- Then, the direct image cannot be calculated unless the whole system containing the cloud is known.
  if (trim(par%source_geometry(1:8)) /= 'external') then
     call peeling_direct_photon(photon,grid)
  endif

  return
  end subroutine gen_photon
end module photon_mod
