module density_mod
contains
  subroutine dust_density(model,grid,i,j,k,opacity)
  !subroutine dust_density(model,r,z,opacity)
!-------------------------------------------------
! Dust Extinction Coefficient
!
!  input  : r,z in kpc
!  input  : i,j,k = cell indeces in r, phi, and z
!           tau_faceon (central, face-on, optical depth of the dust disk)
!           rscale, zscale in kpc
!  output : opacity (extinction coefficient in distance_unit^-1).
!
!  The dust is assumed to be distributed isothermally in the direction
!  perpendicular to the galactic plane with scale height zscale, and
!  exponentially in the radial direction with scale length rscale.
!
! History
!    2015-02-04, minor bug fixed. KI Seon, "else opacity = 0.0_wp" were added.
!    2016-03-06, Make the total mass to be conserved regardless of grid spacing scope.
!-------------------------------------------------
  use define, only : wp, dust_type, grid_type, par
  implicit none
  type(dust_type), intent(in)  :: model
  !real(kind=wp),   intent(in)  :: r,z
  type(grid_type), intent(in)  :: grid
  integer,         intent(in)  :: i,j,k
  real(kind=wp),   intent(out) :: opacity
  real(kind=wp) :: mass, volume, rr1,rr2, zz1, zz2

  volume = (grid%rface(i+1)**2 - grid%rface(i)**2)*(grid%zface(k+1) - grid%zface(k))
  rr1    = grid%rface(i)  /model%rscale
  rr2    = grid%rface(i+1)/model%rscale
  zz1    = grid%zface(k)  /model%zscale
  zz2    = grid%zface(k+1)/model%zscale

  select case(trim(model%name))
  case ('sech')
     ! radially exponential, vertically sech2
     if (grid%rface(i+1) <= model%rmax .and. abs(grid%zface(k)) <= model%zmax .and. grid%zface(k+1) <= model%zmax) then
        mass    = model%tau_faceon*model%rscale**2 * ((1.0_wp+rr1)*exp(-rr1) - (1.0_wp+rr2)*exp(-rr2)) * (tanh(zz2) - tanh(zz1))
        opacity = mass/volume
     else
        opacity = 0.0_wp
     endif
  case ('exponential')
     ! radially exponential, vertically exponential
     if (grid%rface(i+1) <= model%rmax .and. abs(grid%zface(k)) <= model%zmax .and. grid%zface(k+1) <= model%zmax) then
        mass    = model%tau_faceon * model%rscale**2 * ((1.0_wp+rr1)*exp(-rr1) - (1.0_wp+rr2)*exp(-rr2))
        if (zz1 < 0.0_wp .and. zz2 >= 0.0_wp) then
           mass = mass * (2.0_wp - exp(-zz2) - exp(-abs(zz1)))
        else
           mass = mass * abs(exp(-abs(zz1)) - exp(-abs(zz2)))
        endif
        opacity = mass/volume
     else
        opacity = 0.0_wp
     endif
  ! To be done.
  !case ('ring')
  !   if (r <= model%rmax .and. abs(z) <= model%zmax) then
  !      !opacity = model%tau_faceon/(2.0_wp*model%zscale) &
  !      !          *exp(-((r-model%rcenter)**2-model%rcenter**2)/(2.0_wp*model%rscale**2)-abs(z)/model%zscale)
  !      opacity = model%tau_faceon*exp(-(r-model%rcenter)**2/(2.0_wp*model%rscale**2))
  !   else
  !      opacity = 0.0_wp
  !   endif
  case default
     opacity = 0.0_wp
  endselect

  return
  end subroutine dust_density
end module density_mod
