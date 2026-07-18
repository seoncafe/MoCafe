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
  use sed_mod,     only : sample_sed_lambda, sample_ext_lambda, &
                          sample_sed_lambda_u, sample_ext_lambda_u, &
                          sed_wave, sed_sext, sed_albedo, sed_hgg
  use sources_mod, only : use_sources, gen_source_photon, gen_source_photon_qmc
  use compose_mod, only : compose_ext, int_lum_frac, Lpacket_tot
  use qmc_mod,     only : qmc_uniforms, QMC_MAXDIM, QMC_NDIM_USED
  implicit none

  type(grid_type),     intent(inout) :: grid
  !--- photon is intent(inout) so that photon%id, set by the caller before
  !--- this call (LaRT convention), survives; every per-photon field that the
  !--- former intent(out) default-initialized is reset explicitly below.
  type(photon_type),   intent(inout) :: photon

  ! local variables
  real(kind=wp) :: sint,cost,phi,sinp,cosp,rp,tanp
  !--- quasi-random (Owen-scrambled Sobol) launch: uq holds the fixed-layout
  !--- launch coordinates for the global photon number photon%id.  use_qmc =
  !--- .false. leaves every original Mersenne Twister draw untouched
  !--- (bit-identical).
  logical       :: use_qmc
  real(kind=wp) :: uq(QMC_MAXDIM), uorigin

  !--- reset the per-photon fields that carry default initializers in
  !--- photon_type (formerly applied by intent(out) on every call).
  photon%is_external = .false.
  photon%icell_clump = 0
  photon%icell_amr   = 0
  photon%il          = 0
  photon%lambda      = 0.0_wp
  photon%s_ext       = 1.0_wp
  photon%Lpacket     = 1.0_wp

  use_qmc = trim(par%launch_sequence) == 'sobol'
  if (use_qmc) call qmc_uniforms(photon%id - 1_int64, uq(1:QMC_NDIM_USED))

  !--- compose (task C2): an internal source and an isotropic external field in
  !--- one run.  Draw the packet internal-or-external in proportion to
  !--- L_int : L_ext (int_lum_frac = L_int/L_tot).  The external branch places
  !--- the photon on the par%ext_geometry boundary and jumps to the shared
  !--- cell-index / tail; the internal branch falls through to the normal
  !--- single-/multi-source generation below.
  if (compose_ext) then
     if (use_qmc) then
        uorigin = uq(1)
     else
        uorigin = rand_number()
     endif
     if (uorigin >= int_lum_frac) then
        photon%is_external = .true.
        select case (trim(par%ext_geometry))
        case ('cyl')
           if (use_qmc) then
              call external_illumination_cyl_qmc(photon,grid,uq(6),uq(7),uq(8),uq(9),uq(4),uq(5))
           else
              call external_illumination_cyl(photon,grid)
           endif
        case ('rec')
           if (use_qmc) then
              call external_illumination_rec_qmc(photon,grid,uq(6),uq(7),uq(8),uq(4),uq(5))
           else
              call external_illumination_rec(photon,grid)
           endif
        case default
           if (use_qmc) then
              call external_illumination_sph1_qmc(photon,grid,uq(6),uq(7),uq(4),uq(5))
           else
              call external_illumination_sph(photon,grid)
           endif
        end select
        goto 100
     endif
  endif

  !--- multi-source (Stage 6): pick a luminosity-weighted component and set
  !--- position/direction/wavelength from it, then run the common tail
  !--- (clump/amr birth cell, scan reset, direct peel-off).
  if (use_sources) then
     if (use_qmc) then
        call gen_source_photon_qmc(grid, photon, uq)
     else
        call gen_source_photon(grid, photon)
     endif
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
  case ('exponential', 'sech', 'exp_spiral')
     !--- radially exponential disk when source_rscale is set, else plane-uniform;
     !--- vertically exp or sech^2; 'exp_spiral' adds log-spiral arms by rejection.
     if (par%source_rscale > 0.0_wp) then
        if (trim(par%source_geometry) == 'exp_spiral' .and. par%spiral_m > 0) then
           tanp = tan(par%spiral_pitch*pi/180.0_wp)
           do
              rp  = par%source_rscale*rand_r1exp(par%rmax/par%source_rscale)
              phi = twopi*rand_number()
              if (rp <= 0.0_wp) cycle
              if ((1.0_wp+par%spiral_amp)*rand_number() <= &
                  1.0_wp + par%spiral_amp*sin(par%spiral_m*(log(rp)/tanp - phi))) exit
           enddo
        else
           rp  = par%source_rscale*rand_r1exp(par%rmax/par%source_rscale)
           phi = twopi*rand_number()
        endif
        photon%x = rp*cos(phi);  photon%y = rp*sin(phi)
     else
        photon%x = grid%xrange*rand_number()+grid%xmin
        photon%y = grid%yrange*rand_number()+grid%ymin
     endif
     if (trim(par%source_geometry) == 'sech') then
        photon%z = par%source_zscale*rand_sech2(par%zmax/par%source_zscale)
     else
        photon%z = par%source_zscale*rand_zexp(par%zmax/par%source_zscale)
     endif
     photon%wgt = 1.0_wp
  case ('external_cyl')
     photon%is_external = .true.
     if (use_qmc) then
        call external_illumination_cyl_qmc(photon,grid,uq(6),uq(7),uq(8),uq(9),uq(4),uq(5))
     else
        call external_illumination_cyl(photon,grid)
     endif
  case ('external_sph')
     photon%is_external = .true.
     if (use_qmc) then
        call external_illumination_sph1_qmc(photon,grid,uq(6),uq(7),uq(4),uq(5))
     else
        call external_illumination_sph(photon,grid)
     endif
  case ('external_rec')
     photon%is_external = .true.
     if (use_qmc) then
        call external_illumination_rec_qmc(photon,grid,uq(6),uq(7),uq(8),uq(4),uq(5))
     else
        call external_illumination_rec(photon,grid)
     endif
  case default
     photon%x = par%xs_point
     photon%y = par%ys_point
     photon%z = par%zs_point
     photon%wgt = 1.0_wp
  endselect

100 continue
  !--- Cell Index
  photon%icell = floor((photon%x-grid%xmin)/grid%dx)+1
  photon%jcell = floor((photon%y-grid%ymin)/grid%dy)+1
  photon%kcell = floor((photon%z-grid%zmin)/grid%dz)+1
  if (photon%icell == grid%nx+1 .and. photon%kx < 0.0_wp) photon%icell = grid%nx
  if (photon%jcell == grid%ny+1 .and. photon%ky < 0.0_wp) photon%jcell = grid%ny
  if (photon%kcell == grid%nz+1 .and. photon%kz < 0.0_wp) photon%kcell = grid%nz

  !=== set up photon's propagation and direction vetor.
  !--- isotropic emission (skipped for external photons, whose direction the
  !--- external_illumination_* routine already set: source_geometry='external_*'
  !--- for the external-only path, photon%is_external for a composed one).
  if (trim(par%source_geometry(1:8)) /= 'external' .and. .not. photon%is_external) then
     if (use_qmc) then
        cost = 2.0_wp*uq(4)-1.0_wp
        phi  = twopi*uq(5)
     else
        cost = 2.0_wp*rand_number()-1.0_wp
        phi  = twopi*rand_number()
     endif
     sint = sqrt(1.0_wp-cost*cost)
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
     !--- external illumination samples the external-field spectrum; internal
     !--- (point/extended) sources sample the stellar source spectrum.
     if (photon%is_external) then
        if (use_qmc) then
           photon%il = sample_ext_lambda_u(uq(3))
        else
           photon%il = sample_ext_lambda()
        endif
     else
        if (use_qmc) then
           photon%il = sample_sed_lambda_u(uq(3))
        else
           photon%il = sample_sed_lambda()
        endif
     endif
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
  !--- compose (task C2): every composed packet carries L_tot/nphotons, so the
  !--- internal and external energy tallies add up regardless of which branch
  !--- (single/multi internal, or external) set Lpacket above.
  if (compose_ext) photon%Lpacket = Lpacket_tot

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

  !+++ peeled-off.  In compose mode an external photon uses the external-boundary
  !--- direct peel (par%ext_geometry); everything else uses the bound
  !--- peeling_direct_photon (internal, or the external-only binding).
  if (compose_ext .and. photon%is_external) then
     select case (trim(par%ext_geometry))
     case ('cyl')
        call peeling_direct_external_cyl(photon,grid)
     case ('rec')
        call peeling_direct_external_rec(photon,grid)
     case default
        call peeling_direct_external_sph(photon,grid)
     end select
  else
     call peeling_direct_photon(photon,grid)
  endif

  return
  end subroutine gen_photon
end module photon_mod
