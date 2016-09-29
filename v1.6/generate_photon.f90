module photon_mod
contains
  subroutine generate_photon(photon,grid)

  use define
  use random
  use mathlib
  use sersic
  implicit none

  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(in)    :: grid

! local variables
  real(kind=wp) :: sint,cost,phi
  real(kind=wp) :: sint_v,cost_v,phi_v
  real(kind=wp) :: rp,pp,rs_max
  real(kind=wp) :: BulgeLum
  real(kind=wp) :: z_eq,r_eq
  character(len=16) :: which_source
  integer :: is

!--- Set photon location ------------
  which_source = 'disk'
  if (trim(par%source%bulgename) == 'sersic' .and. par%source%BulgeToDisk > 0.0_wp) then
     BulgeLum = par%source%BulgeToDisk/(1.0_wp + par%source%BulgeToDisk)
     if (rand_number() < BulgeLum) which_source = 'bulge'
  endif

  select case(trim(which_source))
  case('disk')
    select case(trim(par%source%diskname))
    case('exponential')
       photon%z = par%source%zscale*rand_zexp(par%source%zmax/par%source%zscale)
       rp       = par%source%rscale*rand_rexp(par%source%rmax/par%source%rscale)
       phi      = twopi*rand_number()
       photon%x = rp * cos(phi)
       photon%y = rp * sin(phi)
    case('sech')
       photon%z = par%source%zscale*rand_sech2(par%source%zmax/par%source%zscale)
       rp       = par%source%rscale*rand_rexp(par%source%rmax/par%source%rscale)
       phi      = twopi*rand_number()
       photon%x = rp * cos(phi)
       photon%y = rp * sin(phi)
    case default
    end select
  case('bulge')
    select case(trim(par%source%bulgename))
    case('sersic')
       do while(.true.)
          ! rp is a spherical radius
          !rp   = par%source%Reff*rand_sersic(par%source%sersic_index)
          rs_max = sqrt(par%source%zmax**2 + par%source%rmax**2)/par%source%Reff
          rp   = par%source%Reff*rand_sersic(par%source%sersic_index,rs_max)
          cost = 2.0_wp*rand_number()-1.0_wp
          sint = sqrt(1.0_wp-cost*cost)
          phi  = twopi*rand_number()
          photon%x = rp * sint*cos(phi)
          photon%y = rp * sint*sin(phi)
          photon%z = rp * cost * par%source%axial_ratio
          ! Here, rp is a cylindrical radius
          rp       = sqrt(photon%x*photon%x + photon%y*photon%y)
          if (rp <= par%source%rmax .and. abs(photon%z) <= par%source%zmax) exit
       enddo
    case default
    end select
  end select

!=== set up photon's direction vetor. (isotropic emission)
  cost_v = 2.0_wp*rand_number()-1.0_wp
  sint_v = sqrt(1.0_wp-cost_v*cost_v)
  !sint = 1.0_wp-cost*cost
  !if (sint <= 0.0_wp) then
  !   sint = 0.0_wp
  !else
  !   sint = sqrt(sint)
  !endif
  phi_v  = twopi*rand_number()

  photon%vx = sint_v*cos(phi_v)
  photon%vy = sint_v*sin(phi_v)
  photon%vz = cost_v

!--- Cell Index
!--- Note: if photon%x = photon%y = 0,
!          then the phi-cell index should be determined by velocity vector
  rp = sqrt(photon%x**2 + photon%y**2)
  if (photon%x == 0.0_wp .and. photon%y == 0.0_wp) then
     pp = arctan(photon%vy, photon%vx)
  else
     pp = arctan(photon%y, photon%x)
  endif
  photon%pcell = int((pp-grid%pmin)/grid%dp)+1

  !if (grid%pow_idx == 1.0_wp .or. grid%pow_idx == -999.0_wp .or. grid%pow_idx == 0.0_wp) then
  !   photon%rcell = int((rp-grid%rmin)/grid%dr(1))+1
  !   photon%zcell = int((photon%z-grid%zmin)/grid%dz(1))+1
  !else
  !   r_eq         = rp**(1.0d0/grid%pow_idx)
  !   photon%rcell = int((r_eq-grid%rmin_eq)/grid%dr_eq)+1
  !   z_eq         = (abs(photon%z))**(1.0d0/grid%pow_idx) * sign(1.0_wp, photon%z)
  !   photon%zcell = int((z_eq-grid%zmin_eq)/grid%dz_eq)+1
  !endif
  if (grid%r_alpha == 1.0_wp .or. grid%r_alpha == 0.0_wp) then
     photon%rcell = int((rp-grid%rmin)/grid%dr(1))+1
  else
     r_eq         = rp**(1.0d0/grid%r_alpha)
     photon%rcell = int((r_eq-grid%rmin_eq)/grid%dr_eq)+1
  endif
  if (grid%z_alpha == 1.0_wp .or. grid%z_alpha == 0.0_wp) then
     photon%zcell = int((photon%z-grid%zmin)/grid%dz(1))+1
  else
     z_eq         = (abs(photon%z))**(1.0d0/grid%z_alpha) * sign(1.0_wp, photon%z)
     photon%zcell = int((z_eq-grid%zmin_eq)/grid%dz_eq)+1
  endif
  photon%wgt = 1.0_wp

  return
  end subroutine generate_photon
end module photon_mod
