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
  photon%icell = int((photon%x-grid%xmin)/grid%dx)+1
  photon%jcell = int((photon%y-grid%ymin)/grid%dy)+1
  !photon%kcell = int((photon%z-grid%zmin)/grid%dz)+1
  if (grid%z_alpha == 1.0_wp .or. grid%z_alpha == 0.0_wp) then
     photon%kcell = int((photon%z-grid%zmin)/grid%dz(1))+1
  else
     z_eq         = (abs(photon%z))**(1.0d0/grid%z_alpha) * sign(1.0_wp, photon%z)
     photon%kcell = int((z_eq-grid%zmin_eq)/grid%dz_eq)+1
  endif

!---
  if (photon%icell < 1 .or. photon%icell > grid%nx+1 .or. &
      photon%jcell < 1 .or. photon%jcell > grid%ny+1 .or. &
      photon%kcell < 1 .or. photon%kcell > grid%nz+1) then
      write(*,*) 'something wrong in generate_photon.f90'
      stop
  endif
  photon%wgt = 1.0_wp

  return
  end subroutine generate_photon
end module photon_mod
