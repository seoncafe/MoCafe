module scatter_mod
  use define
  use random
  use mathlib
  use peelingoff_mod
  public
contains

  subroutine scatter_dust_nostokes(photon,grid)
  implicit none
  ! Calculate new direction cosines
  ! Written by Kwang-il Seon, 2017/09/04
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid

  ! local variables
  real(kind=wp) :: cost,sint,phi,phi2,cosp,sinp
  real(kind=wp) :: kx,ky,kz
  real(kind=wp) :: kx1,ky1,kz1,kr

  photon%nscatt  = photon%nscatt + 1
  par%nscatt_tot = par%nscatt_tot + 1
  photon%wgt     = photon%wgt * par%albedo
  call peeling_scattered_photon(photon,grid)

  cost = rand_henyey_greenstein(par%hgg)
  sint = sqrt(1.0d0-cost**2)

  !--- New phi
  phi  = twopi * rand_number()
  cosp = cos(phi)
  sinp = sin(phi)

  !--- New direction
  kx1 = photon%kx
  ky1 = photon%ky
  kz1 = photon%kz
  kr  = sqrt(kx1*kx1 + ky1*ky1)

  !--- 2020-10-12
  photon%kx = cost*kx1 + sint*(kz1*kx1*cosp - ky1*sinp)/kr
  photon%ky = cost*ky1 + sint*(kz1*ky1*cosp + kx1*sinp)/kr
  photon%kz = cost*kz1 - sint*cosp*kr

  return
end subroutine scatter_dust_nostokes
!=====================================
subroutine scatter_dust_stokes(photon,grid)
  implicit none
  ! Calculate new direction cosines
  ! Written by Kwang-il Seon, 2017/09/04
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid

  ! local variables
  real(kind=wp) :: cost,sint,phi,phi2,cosp,sinp,cos2p,sin2p
  real(kind=wp) :: mr,nr,kr
  real(kind=wp) :: px,py,pz
  real(kind=wp) :: Prand, Pcomp, QoverI, UoverI, S12overS11
  real(kind=wp) :: S11,S12,S33,S34
  real(kind=wp) :: T1,T2,I1,Q1,U1,V1

  photon%nscatt  = photon%nscatt + 1
  par%nscatt_tot = par%nscatt_tot + 1
  photon%wgt     = photon%wgt * par%albedo
  call peeling_scattered_photon(photon,grid)

  !--- Select a new cos(theta) from the numerical table using the alias method. (2021.08.30)
  cost  = rand_alias_linear(scatt_mat%phase_PDF, scatt_mat%alias, scatt_mat%coss, scatt_mat%S11)
  sint  = sqrt(1.0d0-cost**2)

  !--- Calculate Mueller matrix
  call interp_eq(scatt_mat%coss,scatt_mat%S11,cost,S11)
  call interp_eq(scatt_mat%coss,scatt_mat%S12,cost,S12)
  call interp_eq(scatt_mat%coss,scatt_mat%S33,cost,S33)
  call interp_eq(scatt_mat%coss,scatt_mat%S34,cost,S34)
  S12overS11 = S12/S11

  !--- Select a new phi using a rejection method.
  do while(.true.)
    phi    = twopi*rand_number()
    phi2   = 2.0_wp * phi
    QoverI = photon%Q/photon%I
    UoverI = photon%U/photon%I
    Prand  = (1.0_wp + abs(S12overS11) * sqrt(QoverI**2 + UoverI**2))*rand_number()
    Pcomp  =  1.0_wp + S12overS11 * (QoverI*cos(phi2) + UoverI*sin(phi2))
    if (Prand <= Pcomp) exit
  enddo
  cosp = cos(phi)
  sinp = sin(phi)

  ! Which method is faster?
  !cos2p = cos(phi2)
  !sin2p = sin(phi2)
  cos2p = 2.0_wp*cosp*cosp - 1.0_wp
  sin2p = 2.0_wp*sinp*cosp

  T1 =  cos2p*photon%Q + sin2p*photon%U
  T2 = -sin2p*photon%Q + cos2p*photon%U

  !I1 =  S11*photon%I + S12*T1
  !Q1 =  S12*photon%I + S11*T1
  I1 =  S11 + S12*T1
  Q1 =  S12 + S11*T1
  U1 =  S33*T2 + S34*photon%V
  V1 = -S34*T2 + S33*photon%V

  !photon%I  = 1.0_wp
  photon%Q  = Q1/I1
  photon%U  = U1/I1
  photon%V  = V1/I1

  !--- Calculate new reference vectors and propagation vector.
  !--- Do not change the order of the calculation.
  px = cosp * photon%mx + sinp * photon%nx
  py = cosp * photon%my + sinp * photon%ny
  pz = cosp * photon%mz + sinp * photon%nz

  photon%nx = cosp * photon%nx - sinp * photon%mx
  photon%ny = cosp * photon%ny - sinp * photon%my
  photon%nz = cosp * photon%nz - sinp * photon%mz

  photon%mx = cost * px - sint * photon%kx
  photon%my = cost * py - sint * photon%ky
  photon%mz = cost * pz - sint * photon%kz

  photon%kx = sint * px + cost * photon%kx
  photon%ky = sint * py + cost * photon%ky
  photon%kz = sint * pz + cost * photon%kz

  !--- Normalize (Not sure this is really required)
  !kr        = sqrt(photon%kx**2 + photon%ky**2 + photon%kz**2)
  !photon%kx = photon%kx/kr
  !photon%ky = photon%ky/kr
  !photon%kz = photon%kz/kr

  !mr        = sqrt(photon%mx**2 + photon%my**2 + photon%mz**2)
  !photon%mx = photon%mx/mr
  !photon%my = photon%my/mr
  !photon%mz = photon%mz/mr

  !nr        = sqrt(photon%nx**2 + photon%ny**2 + photon%nz**2)
  !photon%nx = photon%nx/nr
  !photon%ny = photon%ny/nr
  !photon%nz = photon%nz/nr
  return
end subroutine scatter_dust_stokes
!-----------------------------------------------------
end module scatter_mod
