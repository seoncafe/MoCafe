module density_mod
contains
  subroutine dust_density(model,r,z,opacity)
!-------------------------------------------------
! Dust Extinction Coefficient
!
!  input  : r,z in kpc
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
!-------------------------------------------------
  use define, only : wp, dust_type, par
  implicit none
  type(dust_type), intent(in) :: model
  real(kind=wp),   intent(in) :: r,z
  real(kind=wp),   intent(out) :: opacity

  select case(trim(model%name))
  case ('sech')
     ! radially exponential, vertically sech2
     if (r <= model%rmax .and. abs(z) <= model%zmax) then
        opacity = model%tau_faceon/(2.0_wp*model%zscale) &
                  *exp(-r/model%rscale)/cosh(z/model%zscale)**2
     else
        opacity = 0.0_wp
     endif
  case ('exponential')
     ! radially exponential, vertically exponential
     if (r <= model%rmax .and. abs(z) <= model%zmax) then
        opacity = model%tau_faceon/(2.0_wp*model%zscale) &
                  *exp(-r/model%rscale-abs(z)/model%zscale)
     else
        opacity = 0.0_wp
     endif
  case default
     opacity = 0.0_wp
  endselect

  return
  end subroutine dust_density
end module density_mod
