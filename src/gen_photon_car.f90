module photon_mod
contains
  subroutine gen_photon(grid,photon)

  use define
  use random
  use peelingoff_mod
  use external_radiation
  use scan_mod,    only : scan_reset_photon
  use clump_mod,   only : active_set_at_point
  use octree_mod,  only : amr_find_leaf
  use sed_mod,     only : sample_sed_lambda, sed_wave, sed_sext, sed_albedo, sed_hgg
  use sources_mod, only : use_sources, gen_source_photon
  implicit none

  type(grid_type),   intent(inout) :: grid
  type(photon_type), intent(out)   :: photon

  ! local variables
  real(kind=wp) :: sint,cost,phi,sinp,cosp,rp

  !--- multi-source (Stage 6): pick a luminosity-weighted component and set
  !--- position/direction/wavelength from it, then run the common tail
  !--- (clump/amr birth cell, scan reset, direct peel-off).
  if (use_sources) then
     call gen_source_photon(grid, photon)
     goto 900
  endif

  !=== set up photon's position vector.
  select case(trim(par%source_geometry))
  case ('uniform')
     if (par%rmax > 0.0_wp) then
        rp   = (rand_number())**(1.0_wp/3.0_wp) * par%rmax
        cost = 2.0_wp*rand_number()-1.0_wp
        sint = sqrt(1.0_wp-cost*cost)
        phi  = twopi*rand_number()
        photon%x = rp*sint*cos(phi)
        photon%y = rp*sint*sin(phi)
        photon%z = rp*cost
        photon%wgt = 1.0_wp
     else if (par%rmax <= 0.0_wp) then
        photon%x = (2.0_wp*rand_number()-1.0_wp)*grid%xmax
        photon%y = (2.0_wp*rand_number()-1.0_wp)*grid%ymax
        photon%z = (2.0_wp*rand_number()-1.0_wp)*grid%zmax
        photon%wgt = 1.0_wp
     endif
  case ('uniform_xy')
     photon%x = grid%xrange*rand_number()+grid%xmin
     photon%y = grid%yrange*rand_number()+grid%ymin
     photon%z = 0.0_wp
     photon%wgt = 1.0_wp
  case ('gaussian')
     photon%x = grid%xrange*rand_number()+grid%xmin
     photon%y = grid%yrange*rand_number()+grid%ymin
     photon%z = par%source_zscale/sqrt(2.0_wp)*rand_gauss()
     photon%wgt = 1.0_wp
  case ('exponential')
     !--- radially exponential disk when source_rscale is set, else plane-uniform.
     if (par%source_rscale > 0.0_wp) then
        rp  = par%source_rscale*rand_r1exp(par%rmax/par%source_rscale)
        phi = twopi*rand_number()
        photon%x = rp*cos(phi);  photon%y = rp*sin(phi)
     else
        photon%x = grid%xrange*rand_number()+grid%xmin
        photon%y = grid%yrange*rand_number()+grid%ymin
     endif
     photon%z = par%source_zscale*rand_zexp(par%zmax/par%source_zscale)
     photon%wgt = 1.0_wp
  case ('external_cyl')
     call external_illumination_cyl(photon,grid)
  case ('external_sph')
     call external_illumination_sph(photon,grid)
  case ('external_rec')
     call external_illumination_rec(photon,grid)
  case default
     photon%x = par%xs_point
     photon%y = par%ys_point
     photon%z = par%zs_point
     photon%wgt = 1.0_wp
  endselect

  !--- Cell Index
  photon%icell = floor((photon%x-grid%xmin)/grid%dx)+1
  photon%jcell = floor((photon%y-grid%ymin)/grid%dy)+1
  photon%kcell = floor((photon%z-grid%zmin)/grid%dz)+1
  if (photon%icell == grid%nx+1 .and. photon%kx < 0.0_wp) photon%icell = grid%nx
  if (photon%jcell == grid%ny+1 .and. photon%ky < 0.0_wp) photon%jcell = grid%ny
  if (photon%kcell == grid%nz+1 .and. photon%kz < 0.0_wp) photon%kcell = grid%nz

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

     if (par%use_stokes) then
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
  endif

  photon%nscatt = 0
  photon%inside = .true.
  photon%albedo = par%albedo
  photon%hgg    = par%hgg

  !--- SED mode: sample the wavelength bin from the source spectrum and set
  !--- the wavelength-dependent dust properties for this photon.
  if (par%use_sed) then
     photon%il     = sample_sed_lambda()
     photon%lambda = sed_wave(photon%il)
     photon%s_ext  = sed_sext(photon%il)
     photon%albedo = sed_albedo(photon%il)
     photon%hgg    = sed_hgg(photon%il)
     !--- physical energy carried by a stellar packet [erg/s].  The Lucy
     !--- driver overrides this per energy pass (star_Lpacket); this default
     !--- covers the standalone save_jlam path (one stellar pass).
     photon%Lpacket = par%luminosity/dble(par%nphotons)
  endif

900 continue
  !--- clump medium: determine the birth clump (0 = vacuum) so the forced
  !--- first scattering and the direct peel-off see the correct starting cell.
  if (trim(par%grid_type) == 'clump') then
     block
       integer(kind=int64) :: active_birth(8)
       integer             :: n_ab
       call active_set_at_point(photon%x, photon%y, photon%z, active_birth, n_ab)
       if (n_ab > 0) photon%icell_clump = int(active_birth(1))
     end block
  else if (trim(par%grid_type) == 'amr') then
     photon%icell_amr = amr_find_leaf(photon%x, photon%y, photon%z)
  endif

  !--- reset each photon's scan accumulators (a,g and/or tau); no-op cost when off
  if (par%use_ag_list .or. par%use_tau_list) call scan_reset_photon()

  !+++ peeled-off.
  call peeling_direct_photon(photon,grid)

  return
  end subroutine gen_photon
end module photon_mod
