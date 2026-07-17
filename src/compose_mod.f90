module compose_mod
!--- Compose an internal source (single point/extended source, or the
!--- multi-source population of sources_mod) WITH an isotropic external
!--- radiation field in one run (MoCafe v2.00, task C2).  Works for both the
!--- monochromatic and the SED path.
!---
!--- Activated only when an internal source coexists with an external field
!--- AND par%source_geometry is NOT the 'external_*' external-only shorthand
!--- (that legacy path is untouched and stays bit-identical).  When active,
!--- each photon is drawn internal-or-external in proportion to L_int : L_ext;
!--- the total luminosity par%luminosity is set to L_tot = L_int + L_ext so the
!--- image normalization already carries the correct absolute scale.
  use define, only : wp
  implicit none
  public

  logical       :: compose_ext  = .false.   ! composition mode active
  real(kind=wp) :: int_lum_frac = 0.0_wp    ! internal luminosity fraction L_int/L_tot
  real(kind=wp) :: Lpacket_tot  = 1.0_wp    ! packet energy L_tot/nphotons [erg/s]
end module compose_mod
