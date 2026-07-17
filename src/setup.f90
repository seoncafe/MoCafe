module setup_mod
contains
!+++++++++++++++++++++++++++++++++++++++++++
  subroutine read_input
  use define
  use utility
  use iofile_mod, only : io_file_extension
  use scan_mod,   only : scan_setup, scan_na, scan_ng, scan_alist, scan_glist, scan_g0, &
                         scan_nt, scan_tlist, scan_s, scan_taumax_ref
  use sed_mod,     only : setup_sed
  use sources_mod, only : setup_sources
  use mpi
  implicit none

! local variables
  character(len=128) :: model_infile, arg
  character(len=256) :: exepath
  integer :: unit, ierr, status, islash

!--- Read in parameters from params.par using namelist command
  namelist /parameters/ par

! read parameters
  if (command_argument_count() >= 1) then
     call get_command_argument(1, model_infile)
  else
     call get_command_argument(0, arg)
     write(*,*) 'Usage: ',trim(arg),' input_file.'
     stop
  endif

  ! newunit specifier is instrodueced in Fortran 2008.
  open(newunit=unit,file=trim(model_infile),status='old')
  read(unit,parameters)
  close(unit)

  !--- Resolve the SEDust directory relative to the executable when the
  !--- user leaves par%sed_workdir blank, so a fresh checkout runs dust emission
  !--- from any working directory without editing paths.  argv(0) is the path to
  !--- MoCafe.x (e.g. '../../MoCafe.x'); take its directory and append SEDust/sed.
  if (len_trim(par%sed_workdir) == 0) then
     call get_command_argument(0, exepath)
     islash = index(trim(exepath), '/', back=.true.)
     if (islash > 0) then
        par%sed_workdir = exepath(1:islash-1) // '/SEDust/sed'
     else
        par%sed_workdir = 'SEDust/sed'
     endif
  endif

  par%nprint = par%no_print
  if (par%nprint >= par%no_photons) par%nprint = par%no_photons/10
  if (par%no_photons < 10) par%nprint = 1

  par%nphotons   = par%no_photons
  par%nscatt_tot = 0.0_wp
  if (par%nx == 1 .or. par%ny == 1 .or. par%nz == 1) par%xyz_symmetry = .false.

  !--- 2024.05.08/2023.10.14
  if (trim(par%source_geometry(1:12)) == 'external_cyl' .and. par%rmax < 0.0_wp) then
     if (len_trim(par%geometry)==0) par%geometry = 'cylinder'
     par%rmax = minval([par%xmax, par%ymax])
  endif

  if (trim(par%source_geometry(1:12)) == 'external_sph' .and. par%rmax < 0.0_wp) then
     if (len_trim(par%geometry)==0) par%geometry = 'sphere'
     par%rmax = minval([par%xmax, par%ymax, par%zmax])
  endif

  !--- Cylinderical Geometry
  if (par%geometry == 'cylinder') then
     if (par%rmax < 0.0_wp) par%rmax = maxval([par%xmax, par%ymax])
     if (par%nr > 1) then
        par%nx = par%nr
        par%ny = par%nr
     endif
     par%nx = maxval([par%nx,par%ny])
     par%ny = par%nx
     if (trim(par%source_geometry(1:8)) == 'external') par%source_geometry = 'external_cyl'
  endif

  !--- Spherical Geometry
  ! if rmax > 0 and geometry is not a cyliner, then the system is regarded as a sphere.
  if (par%rmax > 0.0_wp .and. trim(par%geometry) /= 'cylinder') then
     if (len_trim(par%geometry)==0) par%geometry = 'sphere'
     par%xmax = par%rmax
     par%ymax = par%rmax
     par%zmax = par%rmax
     if (par%nr > 1) then
        par%nx = par%nr
        par%ny = par%nr
        par%nz = par%nr
     endif
     par%nx   = maxval([par%nx,par%ny,par%nz])
     par%ny   = par%nx
     par%nz   = par%nx
     !--- 2023.10.14
     if (trim(par%source_geometry(1:8)) == 'external') par%source_geometry = 'external_sph'
  endif

  !--- MPI-related parameters
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpar%nproc, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpar%p_rank, ierr)
  !The third argument, key, determines the ordering (rank) within each new communicator.
  !call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, mpar%hostcomm, ierr)
  call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, mpar%p_rank, MPI_INFO_NULL, mpar%hostcomm, ierr)
  call MPI_COMM_RANK(mpar%hostcomm, mpar%h_rank, ierr)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, mpar%h_rank, mpar%p_rank, mpar%SAME_HRANK_COMM, ierr)
  call MPI_COMM_SIZE(mpar%SAME_HRANK_COMM, mpar%SAME_HRANK_NPROC, ierr)

  !--- grid-type selector (LaRT v2.10 style): map legacy booleans and validate.
  if (par%use_clump_medium) par%grid_type = 'clump'
  if (par%use_amr_grid)     par%grid_type = 'amr'
  if (trim(par%grid_type) == 'amr') then
     if (trim(par%amr_type) == 'ramses') then
        if (mpar%p_rank == 0) write(*,'(a)') &
           'ERROR: amr_type = ''ramses'' is not read directly; convert with '// &
           'python/AMR_grid/convert_ramses_to_generic.py first.'
        call MPI_FINALIZE(ierr);  stop
     endif
     if (len_trim(par%amr_file) == 0) then
        if (mpar%p_rank == 0) write(*,'(a)') 'ERROR: grid_type=''amr'' requires par%amr_file.'
        call MPI_FINALIZE(ierr);  stop
     endif
  endif
  if (trim(par%grid_type) == 'clump') then
     if (par%rmax <= 0.0_wp) then
        if (mpar%p_rank == 0) write(*,'(a)') 'ERROR: clump medium requires par%rmax > 0 (sphere radius).'
        call MPI_FINALIZE(ierr);  stop
     endif
     if (par%clump_radius <= 0.0_wp) then
        if (mpar%p_rank == 0) write(*,'(a)') 'ERROR: clump medium requires par%clump_radius > 0.'
        call MPI_FINALIZE(ierr);  stop
     endif
     if (len_trim(par%clump_input_file) == 0 .and. par%clump_N_clumps <= 0 .and. &
         par%clump_f_vol <= 0.0_wp .and. par%clump_f_cov <= 0.0_wp) then
        if (mpar%p_rank == 0) write(*,'(a)') &
           'ERROR: clump medium requires one of clump_N_clumps / clump_f_vol / clump_f_cov (or clump_input_file).'
        call MPI_FINALIZE(ierr);  stop
     endif
  endif

  if (len_trim(par%scatt_mat_file) == 0) par%use_stokes = .false.
  if (par%use_stokes) then
     call setup_scattering_matrix(par%scatt_mat_file)
  endif

  !--- source spectrum shape/absolute mode (par%spectrum_type).  It sets the
  !--- column units of every source spectrum file; a physical type makes the
  !--- file luminosity ABSOLUTE (derived from the file integral when the scale is
  !--- unset), while 'shape' keeps the legacy renormalize-to-scale behavior.
  block
    logical :: is_phys, has_global_spec, has_any_srcspec, src_is_absolute, ext_is_absolute
    integer :: is
    select case (trim(par%spectrum_type))
    case ('shape', 'per_um', 'per_ang', 'per_hz', 'per_ev')
       is_phys = trim(par%spectrum_type) /= 'shape'
    case default
       if (mpar%p_rank == 0) write(*,'(a)') &
          'ERROR: par%spectrum_type must be ''shape'' (default), ''per_um'', '// &
          '''per_ang'', ''per_hz'', or ''per_ev''.'
       call MPI_FINALIZE(ierr);  stop
    end select
    if (is_phys .and. .not. par%use_sed) then
       if (mpar%p_rank == 0) write(*,'(a)') &
          'ERROR: a physical par%spectrum_type requires par%use_sed = .true. (SED mode only).'
       call MPI_FINALIZE(ierr);  stop
    endif
    has_global_spec = len_trim(par%source_spectrum) > 0
    has_any_srcspec = .false.
    if (par%nsource > 1) then
       do is = 1, min(par%nsource, MAX_SRC)
          if (len_trim(par%src_spectrum(is)) > 0) has_any_srcspec = .true.
       enddo
    endif
    src_is_absolute = par%use_sed .and. is_phys .and. (has_global_spec .or. has_any_srcspec)
    !--- external SED illumination also DERIVES the luminosity (pi*J*A_surface)
    !--- when the mean intensity J is known: par%ext_intensity is set, or a
    !--- physical par%spectrum_type + par%ext_spectrum file yields J from its integral.
    ext_is_absolute = par%use_sed .and. (trim(par%source_geometry(1:8)) == 'external') .and. &
       ( par%ext_intensity > -900.0_wp .or. (is_phys .and. len_trim(par%ext_spectrum) > 0) )
    !--- luminosity sentinel (default -999 = unset): map to the legacy 1.0
    !--- unless an absolute file spectrum will DERIVE it from the file integral.
    if (par%luminosity < -900.0_wp .and. .not. src_is_absolute .and. .not. ext_is_absolute) &
       par%luminosity = 1.0_wp
  end block

  !--- resolve the physical length scale (code unit -> cm) before the SED /
  !--- external-field setup, which converts the illuminated boundary area to cm^2
  !--- for the pi*J*A luminosity.
  select case(trim(par%distance_unit))
     case ('kpc')
        par%distance2cm = kpc2cm
     case ('pc')
        par%distance2cm = pc2cm
     case ('au')
        par%distance2cm = au2cm
     case ('')
        par%distance2cm = 1.0_wp
     case default
        par%distance2cm = kpc2cm
  end select

  !--- SED (multi-wavelength) mode: Stage 1 restrictions and setup.
  if (par%use_sed) then
     if (par%use_stokes) then
        if (mpar%p_rank == 0) write(*,'(a)') &
           'ERROR: par%use_sed is incompatible with par%use_stokes (Henyey-Greenstein path only, Stage 1).'
        call MPI_FINALIZE(ierr);  stop
     endif
     if (par%use_ag_list .or. par%use_tau_list) then
        if (mpar%p_rank == 0) write(*,'(a)') &
           'ERROR: par%use_sed is incompatible with the (a,g)/tau scans (use_ag_list / use_tau_list).'
        call MPI_FINALIZE(ierr);  stop
     endif
     if (trim(par%source_geometry(1:8)) == 'external' .and. &
         trim(par%spectrum_type) /= 'shape' .and. len_trim(par%distance_unit) == 0 .and. &
         mpar%p_rank == 0) write(*,'(a)') &
        'WARNING: external SED with a physical par%spectrum_type and blank par%distance_unit: '// &
        'the surface-area conversion uses code units (par%distance2cm = 1).'
     call setup_sed()
     !--- Stage 6: multiple stellar source components (activated when nsource>1).
     if (par%nsource > 1) call setup_sources()
  endif

  !--- monochromatic multiple internal sources (the SED path is handled above).
  !--- setup_sources sets use_sources = .true. so gen_photon dispatches to the
  !--- luminosity-weighted component sampler.
  if (.not. par%use_sed .and. par%nsource > 1) call setup_sources()

  !--- compose (task C2): an internal source coexisting with an external field
  !--- (not the 'external_*' external-only shorthand).  Runs after the internal
  !--- luminosity is fixed (par%luminosity = L_int at this point).
  call setup_compose()

  !--- Stage 2: J_lambda tally restrictions.
  if (par%save_jlam) then
     if (.not. par%use_sed) then
        if (mpar%p_rank == 0) write(*,'(a)') 'ERROR: par%save_jlam requires par%use_sed = .true.'
        call MPI_FINALIZE(ierr);  stop
     endif
     if (.not. (trim(par%grid_type) == 'car' .or. trim(par%grid_type) == 'amr') .or. par%xy_periodic) then
        if (mpar%p_rank == 0) write(*,'(a)') &
           'ERROR: par%save_jlam supports the Cartesian and AMR grids only (no clump / xy_periodic).'
        call MPI_FINALIZE(ierr);  stop
     endif
  endif

  !--- Stage 3/4: dust thermal emission.
  if (par%use_dustemis) then
     if (.not. par%use_sed) then
        if (mpar%p_rank == 0) write(*,'(a)') 'ERROR: par%use_dustemis requires par%use_sed = .true.'
        call MPI_FINALIZE(ierr);  stop
     endif
     select case (trim(par%dust_emission_method))
     case ('lucy')
        if (.not. par%save_jlam) then
           if (mpar%p_rank == 0) write(*,'(a)') &
              'ERROR: dust_emission_method=''lucy'' requires par%save_jlam = .true.'
           call MPI_FINALIZE(ierr);  stop
        endif
        select case (trim(par%dust_model_sed))
        case ('astrodust', 'dl07')
           if (len_trim(par%sed_qtable) == 0 .or. len_trim(par%sed_sizedist) == 0) then
              if (mpar%p_rank == 0) write(*,'(3a)') &
                 'ERROR: dust_model_sed=''', trim(par%dust_model_sed), &
                 ''' requires par%sed_qtable and par%sed_sizedist.'
              call MPI_FINALIZE(ierr);  stop
           endif
        case ('zubko')
           if (len_trim(par%sed_zubko_config) == 0 .or. len_trim(par%sed_zubko_dir) == 0) then
              if (mpar%p_rank == 0) write(*,'(a)') &
                 'ERROR: dust_model_sed=''zubko'' requires par%sed_zubko_config and par%sed_zubko_dir.'
              call MPI_FINALIZE(ierr);  stop
           endif
        case default
           if (mpar%p_rank == 0) write(*,'(3a)') &
              'ERROR: par%dust_model_sed = ''', trim(par%dust_model_sed), &
              ''' unknown (use ''astrodust'', ''dl07'', or ''zubko'').'
           call MPI_FINALIZE(ierr);  stop
        end select
     case ('bw01')
        !--- Bjorkman & Wood needs only the mixture opacity (par%kext_file).
     case default
        if (mpar%p_rank == 0) write(*,'(3a)') &
           'ERROR: par%dust_emission_method = ''', trim(par%dust_emission_method), &
           ''' (use ''lucy'' or ''bw01'').'
        call MPI_FINALIZE(ierr);  stop
     end select
     if (par%luminosity <= 1.0_wp .and. mpar%p_rank == 0) write(*,'(a)') &
        'WARNING: par%luminosity <= 1; set it to the physical stellar luminosity [erg/s] for absolute dust temperatures.'
  endif

  !--- Scan modes: (a,g) -- Seon 2010, PKAS, 25, 177; tau (polychromatic) -- Jonsson
  !--- 2006, MNRAS, 372, 2.  The two compose into scatt(x,y,a,g,tau).
  if (par%use_ag_list .or. par%use_tau_list) then
     if (par%use_stokes) then
        if (mpar%p_rank == 0) write(*,'(a)') &
           'ERROR: scan (use_ag_list / use_tau_list) is incompatible with par%use_stokes (Henyey-Greenstein path only).'
        call MPI_FINALIZE(ierr)
        stop
     endif
     !--- (tau scan + external illumination is now supported: the external
     !--- direct peel writes the tau-axis direct image; see external_radiation.f90.)
     if (.not. par%use_reduced_wgt .and. mpar%p_rank == 0) write(*,'(a)') &
        'WARNING: scan forces a continuous/immortal forward pass (par%use_reduced_wgt ignored).'
     call scan_setup()
     if (mpar%p_rank == 0) then
        if (par%use_ag_list) then
           write(*,'(a,i0,a,i0,a)') 'Scan: (a,g) = ', scan_na, ' albedo x ', scan_ng, ' hgg values'
           write(*,'(a,f7.4)')      'Scan: simulated g0 (par%hgg) = ', scan_g0
           write(*,'(a,10f6.2)')    'Scan: albedo_list = ', scan_alist(1:min(scan_na,10))
           write(*,'(a,10f6.2)')    'Scan: hgg_list    = ', scan_glist(1:min(scan_ng,10))
           if (scan_g0 < 0.4_wp .or. scan_g0 > 0.5_wp) write(*,'(a)') &
              'Scan: note - Seon 2010 finds g0 = 0.4-0.5 most efficient for g = 0.0-0.9.'
        endif
        if (par%use_tau_list) then
           write(*,'(a,i0,a)')   'Scan: tau = ', scan_nt, ' optical-depth values (Jonsson 2006)'
           write(*,'(a,f9.4)')   'Scan: reference taumax = ', scan_taumax_ref
           write(*,'(a,10f8.3)') 'Scan: tau_list = ', scan_tlist(1:min(scan_nt,10))
           write(*,'(a,10f8.3)') 'Scan: tau/tau0 = ', scan_s(1:min(scan_nt,10))
           if (any(scan_s(1:scan_nt) < 0.5_wp) .or. any(scan_s(1:scan_nt) > 3.0_wp)) write(*,'(a)') &
              'Scan: WARNING - some tau/taumax outside [0.5,3]; variance grows away from the reference (Jonsson Fig.1).'
        endif
     endif
  endif

  if (len_trim(par%out_file) == 0) then
     par%base_name = trim(get_base_input_name(model_infile))
     par%out_file  = trim(par%base_name)//'_obs'//trim(io_file_extension(par%file_format))
  else
     par%base_name = trim(get_base_name(trim(par%out_file)))
     !--- par%file_format is authoritative. If par%out_file carries a different
     !    extension than par%file_format implies, rewrite the extension
     !    (preserving the directory prefix) and warn the user. This avoids the
     !    silent mismatch where the main file would land in one format and
     !    derived files (peeloff, sightline, ...) in another, since the
     !    derived names always use par%file_format.
     block
       character(len=128) :: expected_ext, current_ext, dir_part, new_out_file
       integer :: idir, lo
       expected_ext = io_file_extension(par%file_format)
       lo           = len_trim(par%out_file)
       current_ext  = ''
       if (lo >= 8) then
          if (par%out_file(lo-7:lo) == '.fits.gz') current_ext = '.fits.gz'
       endif
       if (trim(current_ext) == '' .and. lo >= 5) then
          if (par%out_file(lo-4:lo) == '.hdf5') current_ext = '.hdf5'
          if (par%out_file(lo-4:lo) == '.fits') current_ext = '.fits'
       endif
       if (trim(current_ext) == '' .and. lo >= 3) then
          if (par%out_file(lo-2:lo) == '.h5') current_ext = '.h5'
       endif
       if (trim(current_ext) /= trim(expected_ext)) then
          idir = index(par%out_file, '/', back=.true.)
          if (idir > 0) then
             dir_part = par%out_file(1:idir)
          else
             dir_part = ''
          endif
          new_out_file = trim(dir_part)//trim(par%base_name)//trim(expected_ext)
          if (mpar%p_rank == 0) then
             write(*,'(a)')  '---'
             write(*,'(4a)') 'WARNING: par%out_file (', trim(par%out_file), &
                  ') does not match par%file_format = ', trim(par%file_format)
             write(*,'(2a)') '         Overriding par%out_file to: ', trim(new_out_file)
             write(*,'(a)')  '         (set par%file_format explicitly to suppress this warning.)'
             write(*,'(a)')  '---'
          endif
          par%out_file = new_out_file
       endif
     end block
  endif

  if (mpar%p_rank == 0) then
     write(*,'(a)')        ''
     write(*,'(3a)')       '+++++ ',trim(model_infile),' +++++'
     write(*,'(2a)')       ' >>> START @ ', get_date_time()
     write(*,'(a,L1)')     'Use_Master_Slave          : ', par%use_master_slave
     write(*,'(a,2f7.4)')  'Dust Parameters(a, g)     : ', par%albedo, par%hgg
     write(*,'(a,es12.3)') 'Dust Extinction per H     : ', par%cext_dust
     write(6,'(a,i14)')    'Total photons             : ', par%nphotons
  endif

  return
  end subroutine read_input
  !---------------------------------------
  subroutine setup_scattering_matrix(scatt_mat_file)
  use define
  use mathlib
  use memory_mod
  use random
  use mpi
  implicit none
  character(len=*), intent(in) :: scatt_mat_file

  character(len=128) :: string_tmp
  real(kind=wp) :: lambda,cext,albedo,hgg
  real(kind=wp) :: S11_norm
  integer       :: i, unit, ierr

  !--- setup Mueller Matrix
  !--- In this code, we will assume spherical dust grains.
  !--- For non-spherical grains, we also need to consider componets (i.e., S22, S21, S13, etc)
  if (len_trim(scatt_mat_file) > 0) then
     if (mpar%h_rank == 0) then
        open(newunit=unit,file=trim(scatt_mat_file),status='old')
        read(unit,*) string_tmp
        read(unit,*) lambda,cext,albedo,hgg,scatt_mat%nPDF
        !read(unit,*) lambda,albedo,hgg,scatt_mat%nPDF
        par%albedo  = albedo
        par%hgg     = hgg
        ! par%lambda0 and lambda is in units of micron.
        par%lambda0 = lambda
        read(unit,*) string_tmp
     endif
     call MPI_BCAST(par%albedo,     1, MPI_DOUBLE_PRECISION, 0, mpar%hostcomm, ierr)
     call MPI_BCAST(par%hgg,        1, MPI_DOUBLE_PRECISION, 0, mpar%hostcomm, ierr)
     call MPI_BCAST(par%lambda0,    1, MPI_DOUBLE_PRECISION, 0, mpar%hostcomm, ierr)
     call MPI_BCAST(scatt_mat%nPDF, 1, MPI_INTEGER,          0, mpar%hostcomm, ierr)

     call create_shared_mem(scatt_mat%coss,      [scatt_mat%nPDF])
     call create_shared_mem(scatt_mat%S11,       [scatt_mat%nPDF])
     call create_shared_mem(scatt_mat%S12,       [scatt_mat%nPDF])
     call create_shared_mem(scatt_mat%S33,       [scatt_mat%nPDF])
     call create_shared_mem(scatt_mat%S34,       [scatt_mat%nPDF])
     call create_shared_mem(scatt_mat%alias,     [scatt_mat%nPDF-1])
     call create_shared_mem(scatt_mat%phase_PDF, [scatt_mat%nPDF-1])

     if (mpar%h_rank == 0) then
        do i=1,scatt_mat%nPDF
           read(unit,*) scatt_mat%coss(i),scatt_mat%S11(i),scatt_mat%S12(i),scatt_mat%S33(i),scatt_mat%S34(i)
        enddo
        close(unit)

        ! Normalize the Mueller matrixes
        S11_norm      = calc_Integral(scatt_mat%coss,scatt_mat%S11)
        scatt_mat%S11 = scatt_mat%S11/S11_norm
        scatt_mat%S12 = scatt_mat%S12/S11_norm
        scatt_mat%S33 = scatt_mat%S33/S11_norm
        scatt_mat%S34 = scatt_mat%S34/S11_norm

        ! setup alias (2021.08.30)
        do i=1,scatt_mat%nPDF-1
           !--- bug-fixed, 2025.03.16
           !scatt_mat%phase_PDF(i) = (scatt_mat%S11(i) + scatt_mat%S12(i))/2.0_wp
           scatt_mat%phase_PDF(i) = (scatt_mat%S11(i) + scatt_mat%S11(i+1))/2.0_wp
        enddo
        scatt_mat%phase_PDF = scatt_mat%phase_PDF / sum(scatt_mat%phase_PDF)
        call random_alias_setup(scatt_mat%phase_PDF, scatt_mat%alias)
     endif
  else
     scatt_mat%nPDF = 0
  endif
  end subroutine setup_scattering_matrix
  !---------------------------------------
  subroutine setup_procedure
  use define
  use raytrace
  use raytrace_clump_mod
  use raytrace_amr_mod
  use peelingoff_mod
  use external_radiation
  use scatter_mod
  use run_simulation_mod
  use write_mod
  implicit none

  !--- Initialize Random Number Generator
  call init_random_seed(par%iseed)

  !--- procedure pointer for raytrace routine
  raytrace_to_edge => raytrace_to_edge_car

  if (par%xyz_symmetry) then
     raytrace_to_tau => raytrace_to_tau_car_xyzsym
  else if (par%xy_periodic) then
     if (par%nx == 1 .and. par%ny == 1) then
        raytrace_to_tau => raytrace_to_tau_car_zonly
     else
        raytrace_to_tau => raytrace_to_tau_car_xyper
     endif
  else
     raytrace_to_tau => raytrace_to_tau_car
  endif

  !--- procedure pointer for peelingoff_mod
  if (par%use_stokes) then
    peeling_direct_photon    => peeling_direct_photon_nostokes
    peeling_scattered_photon => peeling_scattered_photon_stokes
  else
    select case(trim(par%source_geometry))
    case ('external_sph', 'external_sph_radial')
       peeling_direct_photon => peeling_direct_external_sph
    case ('external_cyl')
       peeling_direct_photon => peeling_direct_external_cyl
    case ('external_rec')
       peeling_direct_photon => peeling_direct_external_rec
    case default
       peeling_direct_photon => peeling_direct_photon_nostokes
    end select

    peeling_scattered_photon => peeling_scattered_photon_nostokes
  endif

  !--- procedure pointer for scattering routine
  if (par%use_sed) then
     !--- SED (multi-wavelength) mode: H-G no-Stokes path with
     !--- wavelength-dependent albedo/g and s_ext-rescaled optical depths.
     scattering               => scatter_dust_nostokes_sed
     peeling_scattered_photon => peeling_scattered_photon_nostokes_sed
     !--- external sources keep the (now SED-aware) external direct peel bound
     !--- above; only an internal source uses the SED point-source direct peel.
     if (trim(par%source_geometry(1:8)) /= 'external') &
        peeling_direct_photon => peeling_direct_photon_nostokes_sed
  else if (par%use_tau_list) then
     !--- polychromatic (tau) scan, optionally combined with (a,g).  H-G no-Stokes
     !--- path only; the (a,g) axes collapse to length 1 when use_ag_list = .false.
     scattering               => scatter_dust_nostokes_scan
     peeling_scattered_photon => peeling_scattered_photon_nostokes_scan_tau
     !--- External sources keep the (now tau-aware) external direct peel bound
     !--- above; only an internal source uses the point-source tau direct peel.
     if (trim(par%source_geometry(1:8)) /= 'external') &
        peeling_direct_photon => peeling_direct_photon_nostokes_tau
  else if (par%use_ag_list) then
     !--- (albedo, asymmetry-factor) scan: H-G no-Stokes path only.
     scattering               => scatter_dust_nostokes_scan
     peeling_scattered_photon => peeling_scattered_photon_nostokes_scan
  else if (par%use_stokes) then
     scattering => scatter_dust_stokes
  else
     scattering => scatter_dust_nostokes
  endif

  !--- procedure pointer for simulation run
  if (par%use_master_slave) then
     run_simulation => run_master_slave
  else
     run_simulation => run_equal_number
  endif

  !--- procedure pointer to write output
  write_output => write_output_car

  !--- clump medium: override only the raytrace pointers.  Scattering and
  !--- peeling-off are SHARED with the Cartesian dust path (the physics is
  !--- identical; only the optical-depth integration through the medium
  !--- differs).  This binding runs last so it composes with the (a,g)/tau
  !--- scan and Stokes scattering selections above.
  if (trim(par%grid_type) == 'clump') then
     raytrace_to_tau  => raytrace_to_tau_clump
     raytrace_to_edge => raytrace_to_edge_clump
  else if (trim(par%grid_type) == 'amr') then
     raytrace_to_tau  => raytrace_to_tau_amr
     raytrace_to_edge => raytrace_to_edge_amr
  endif

  end subroutine setup_procedure
  !=================================================
  !--- compose (task C2): decide whether an internal source composes with an
  !--- isotropic external field, and if so build the luminosity split.  Leaves
  !--- compose_ext = .false. (no-op) for every legacy path.
  subroutine setup_compose
  use define
  use compose_mod
  use sed_mod, only : setup_sed_external, ext_surface_area
  use mpi
  implicit none
  real(kind=wp) :: L_int, L_ext, L_tot, J_ext, A_surf
  logical :: ext_on
  integer :: ierr

  compose_ext  = .false.
  int_lum_frac = 0.0_wp
  Lpacket_tot  = par%luminosity/dble(par%nphotons)

  !--- the 'external_*' shorthand is the external-only path: never composed.
  if (trim(par%source_geometry(1:8)) == 'external') return
  if (par%nsource < 1) return

  !--- is an external field switched on?
  if (par%use_sed) then
     ext_on = (par%ext_intensity > -900.0_wp) .or. (len_trim(par%ext_spectrum) > 0) &
              .or. (par%ext_tstar > 0.0_wp)
  else
     ext_on = (par%ext_intensity > 0.0_wp)
  endif
  if (.not. ext_on) return

  !--- restrictions outside the scope of this task.
  if (par%use_stokes) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: composing an internal source with an external field is not supported with par%use_stokes.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (par%use_ag_list .or. par%use_tau_list) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: composing an internal source with an external field is not supported with the '// &
        '(a,g)/tau scans (use_ag_list / use_tau_list).'
     call MPI_FINALIZE(ierr);  stop
  endif

  !--- resolve the external-field boundary geometry (default from geometry/rmax).
  if (len_trim(par%ext_geometry) == 0) then
     if (trim(par%geometry) == 'cylinder') then
        par%ext_geometry = 'cyl'
     else if (par%rmax > 0.0_wp) then
        par%ext_geometry = 'sph'
     else
        par%ext_geometry = 'rec'
     endif
  endif
  select case (trim(par%ext_geometry))
  case ('sph')
     if (.not. (par%rmax > 0.0_wp)) then
        if (mpar%p_rank == 0) write(*,'(a)') &
           'ERROR: par%ext_geometry = ''sph'' requires par%rmax > 0.'
        call MPI_FINALIZE(ierr);  stop
     endif
  case ('cyl', 'rec')
     !--- ok
  case default
     if (mpar%p_rank == 0) write(*,'(3a)') &
        'ERROR: par%ext_geometry = ''', trim(par%ext_geometry), ''' (use ''sph'', ''cyl'', or ''rec'').'
     call MPI_FINALIZE(ierr);  stop
  end select

  !--- internal luminosity (already fixed by setup_sed / setup_sources).
  L_int = par%luminosity

  !--- external luminosity L_ext = pi * J * A_surface.
  if (par%use_sed) then
     !--- also builds the external-field spectrum alias (sample_ext_lambda).
     call setup_sed_external(compose_mode=.true., lum_out=L_ext)
  else
     J_ext  = par%ext_intensity
     A_surf = ext_surface_area(par%ext_geometry)
     L_ext  = pi * J_ext * A_surf
  endif

  if (.not. (L_ext > 0.0_wp)) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: composed external luminosity L_ext = pi*J*A is not positive (check par%ext_intensity / ext_spectrum).'
     call MPI_FINALIZE(ierr);  stop
  endif

  L_tot          = L_int + L_ext
  int_lum_frac   = L_int/L_tot
  Lpacket_tot    = L_tot/dble(par%nphotons)
  par%luminosity = L_tot
  compose_ext    = .true.

  if (mpar%p_rank == 0) then
     write(*,'(a)')        '--- Compose: internal source + external field (task C2) ---'
     write(*,'(2a)')       'external boundary geometry: ', trim(par%ext_geometry)
     write(*,'(a,es12.4)') 'L_internal (erg/s)        : ', L_int
     write(*,'(a,es12.4)') 'L_external (erg/s)        : ', L_ext
     write(*,'(a,es12.4)') 'L_total    (erg/s)        : ', L_tot
     write(*,'(a,f9.5)')   'internal photon fraction  : ', int_lum_frac
  endif
  end subroutine setup_compose
  !=================================================
end module setup_mod
