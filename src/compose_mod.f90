module compose_mod
!--- Compose an internal source (a single point/extended source, or the
!--- multi-source population of sources_mod) WITH an isotropic external
!--- radiation field in one run.
!---
!--- Activated only when an internal source coexists with an external field
!--- AND par%source_geometry is NOT the 'external_*' external-only shorthand
!--- (that path is untouched and stays bit-identical).  When active, each
!--- photon is drawn internal-or-external in proportion to L_int : L_ext; the
!--- total luminosity par%luminosity is set to L_tot = L_int + L_ext so that
!--- the image normalization already carries the correct absolute scale.
  use define, only : wp, par, pi, twopi, fourpi
  implicit none
  public

  logical       :: compose_ext  = .false.   ! composition mode active
  real(kind=wp) :: int_lum_frac = 0.0_wp    ! internal luminosity fraction L_int/L_tot

contains
  !---------------------------------------------------------------
  !--- surface area [cm^2] of the external-field boundary.
  function ext_surface_area(geom) result(A_surf)
  implicit none
  character(len=*), intent(in) :: geom
  real(kind=wp) :: A_surf, Lx, Ly, Lz
  select case (trim(geom))
  case ('sph')
     A_surf = fourpi*(par%rmax*par%distance2cm)**2
  case ('rec')
     Lx = 2.0_wp*par%xmax*par%distance2cm
     Ly = 2.0_wp*par%ymax*par%distance2cm
     Lz = 2.0_wp*par%zmax*par%distance2cm
     A_surf = 2.0_wp*(Lx*Ly + Ly*Lz + Lx*Lz)
  case ('cyl')
     Lz = 2.0_wp*par%zmax*par%distance2cm
     A_surf = twopi*par%rmax*par%distance2cm*Lz + twopi*(par%rmax*par%distance2cm)**2
  case default
     A_surf = fourpi*(par%rmax*par%distance2cm)**2
  end select
  end function ext_surface_area
end module compose_mod
