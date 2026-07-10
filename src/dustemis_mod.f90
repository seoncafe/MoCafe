module dustemis_mod
!--- Dust thermal emission via the SEDust library (MoCafe v2.00, Stage 3).
!--- Mode 1 (Lucy 1999): from the per-cell mean intensity J_lambda tallied by
!--- jtally_mod, compute each cell's dust emission spectrum with SEDust
!--- (equilibrium + stochastically heated grains + PAHs), then produce the
!--- emergent dust thermal SED.
!---
!--- This first version is NON-ITERATIVE: a single stellar pass sets J_lambda,
!--- and the dust emission is computed once.  It is exact in the limit that the
!--- dust is optically thin to its own (re-emitted) radiation.  Dust
!--- self-absorption (the Lucy iteration proper) is a later refinement.
!---
!--- Units.  MoCafe J_lambda is assembled in CGS [erg s^-1 cm^-2 sr^-1 um^-1]
!--- from the raw pathlength tally jt_sum = Sum(wgt*dl):
!---   J_cgs(il) = E_p * jt_sum(il,cell) / (4 pi V dlam) / dist_cm^2 ,
!--- with E_p = par%luminosity/nphotons [erg/s].  SEDust's Planck (bbody) is SI
!--- [W m^-2 sr^-1 m^-1], so J_SI = 1.0e3 * J_cgs.  dust_emission returns
!--- lamI_total(lambda) = lambda * (emission per H) whose SHAPE gives the
!--- spectral PDF; the absolute per-cell luminosity is set to the locally
!--- absorbed power (radiative equilibrium), which is exact by construction and
!--- sidesteps the SEDust absolute-normalization convention.
  use define
  use dust_lib, only : dust_model_t, build_astrodust, build_dl07, build_zubko, &
                       dust_emission, dust_nlam, dust_lambda
  implicit none
  private
  public :: setup_dustemis, compute_dustemis, write_dustemis

  type(dust_model_t)         :: dmodel
  integer                    :: nl_sed = 0            ! SEDust wavelength grid length
  real(kind=wp), allocatable :: lam_sed(:)           ! SEDust lambda grid [um]
  real(kind=wp), allocatable :: cell_Lemit(:)        ! per-cell emitted (=absorbed) L [erg/s]
  real(kind=wp), allocatable :: cell_pdf(:,:)        ! (nl_mc, ncell) emission energy fraction per bin
  real(kind=wp), allocatable :: cell_Teq(:)          ! per-cell luminosity-weighted equilibrium-ish T [K] (diagnostic)
  real(kind=wp) :: Labs_total = 0.0_wp               ! total absorbed L over the grid [erg/s]
  integer :: ncell_tot = 0

contains
  !---------------------------------------------------------------
  subroutine setup_dustemis(grid)
  use mpi
  use ifport, only : chdir, getcwd
  implicit none
  type(grid_type), intent(in) :: grid
  integer :: ierr, cstat
  character(len=512) :: cwd_save

  !--- SEDust reads its dielectric tables via paths hard-coded relative to its
  !--- sed/ directory ('../data/dielectric/...'), so build the model from
  !--- par%sed_workdir (default SEDust/sed) and restore the working directory.
  !--- par%sed_qtable / par%sed_sizedist are given as absolute paths.
  !--- chdir/getcwd are the Intel IFPORT integer functions (return 0 on success).
  cstat = getcwd(cwd_save)
  if (len_trim(par%sed_workdir) > 0) then
     cstat = chdir(trim(par%sed_workdir))
     if (cstat /= 0 .and. mpar%p_rank == 0) write(*,'(3a)') &
        'WARNING: could not chdir to par%sed_workdir = ''', trim(par%sed_workdir), &
        ''' (SEDust dielectric files may not be found).'
  endif

  select case (trim(par%dust_model_sed))
  case ('astrodust')
     call build_astrodust(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
                          par%sed_NT, par%sed_Tlo, par%sed_Thi)
  case default
     cstat = chdir(trim(cwd_save))
     if (mpar%p_rank == 0) write(*,'(3a)') &
        'ERROR: par%dust_model_sed = ''', trim(par%dust_model_sed), &
        ''' not wired yet (use ''astrodust'').'
     call MPI_FINALIZE(ierr);  stop
  end select

  cstat = chdir(trim(cwd_save))

  nl_sed = dust_nlam(dmodel)
  lam_sed = dust_lambda(dmodel)          ! allocatable assignment
  ncell_tot = grid%nx*grid%ny*grid%nz

  if (mpar%p_rank == 0) then
     write(*,'(a)')      '--- Dust thermal emission (SEDust, Mode 1 Lucy, non-iterative) ---'
     write(*,'(2a)')     'dust model                : ', trim(dmodel%name)
     write(*,'(a,i6)')   'SEDust wavelength grid    : ', nl_sed
     write(*,'(a,2es11.3)') 'SEDust lambda range (um)  : ', lam_sed(1), lam_sed(nl_sed)
  endif
  end subroutine setup_dustemis

  !---------------------------------------------------------------
  !--- compute per-cell emission from the reduced J tally.
  subroutine compute_dustemis(grid)
  use mpi
  use sed_mod,    only : sed_nlam, sed_wave, sed_dwave, sed_sext, sed_albedo
  use jtally_mod, only : jt_sum
  implicit none
  type(grid_type), intent(in) :: grid
  real(kind=wp), allocatable :: Jcgs(:), Jsi(:), Jsed(:), lamI_sed(:), lamI_mc(:), emis(:)
  real(kind=wp) :: E_p, vol, dist2, jnorm, Labs, esum, lam_mean
  integer :: nl_mc, i, j, k, ic, il, ierr, ndone, nmine

  nl_mc = sed_nlam
  ndone = 0
  nmine = (ncell_tot - mpar%p_rank + mpar%nproc - 1)/mpar%nproc
  E_p   = par%luminosity/dble(par%nphotons)
  vol   = grid%dx*grid%dy*grid%dz
  dist2 = par%distance2cm**2

  if (.not. allocated(cell_Lemit)) allocate(cell_Lemit(ncell_tot))
  if (.not. allocated(cell_pdf))   allocate(cell_pdf(nl_mc, ncell_tot))
  if (.not. allocated(cell_Teq))   allocate(cell_Teq(ncell_tot))
  cell_Lemit(:)  = 0.0_wp
  cell_pdf(:,:)  = 0.0_wp
  cell_Teq(:)    = 0.0_wp

  allocate(Jcgs(nl_mc), Jsi(nl_mc), Jsed(nl_sed), lamI_sed(nl_sed), lamI_mc(nl_mc), emis(nl_mc))

  if (mpar%p_rank == 0) write(*,'(a,i0,a)') 'computing dust emission for ~', nmine, ' cells/rank ...'

  !--- distribute cells across MPI ranks (round-robin on the linear index).
  do ic = mpar%p_rank+1, ncell_tot, mpar%nproc
     ndone = ndone + 1
     if (mpar%p_rank == 0 .and. nmine >= 500 .and. mod(ndone, nmine/5) == 0) &
        write(*,'(a,i0,a,i0,a)') '   dust emission: rank0 ', ndone, '/', nmine, ' cells'
     k = (ic-1)/(grid%nx*grid%ny) + 1
     j = mod((ic-1)/grid%nx, grid%ny) + 1
     i = mod(ic-1, grid%nx) + 1
     if (grid%rhokap(i,j,k) <= 0.0_wp) cycle
     if (all(jt_sum(:,i,j,k) <= 0.0_wp)) cycle

     !--- physical mean intensity in this cell (SI units for SEDust).
     jnorm = E_p/(fourpi*vol*dist2)
     do il = 1, nl_mc
        Jcgs(il) = jt_sum(il,i,j,k)*jnorm/sed_dwave(il)     ! erg/s/cm^2/sr/um
     enddo
     Jsi(:) = Jcgs(:)*1.0e3_wp                              ! -> W/m^2/sr/m

     !--- absorbed power in this cell [erg/s] (radiative equilibrium => emitted).
     Labs = E_p*grid%rhokap(i,j,k) * &
            sum(jt_sum(:,i,j,k)*sed_sext(:)*(1.0_wp - sed_albedo(:)))
     if (Labs <= 0.0_wp) cycle

     !--- SEDust emission spectrum (shape): resample J onto the SEDust grid,
     !--- solve, resample lamI back onto the MoCafe grid.
     call resample_to_sed(sed_wave, Jsi, lam_sed, Jsed)
     call dust_emission(dmodel, Jsed, lamI_sed)
     call resample_from_sed(lam_sed, lamI_sed, sed_wave, lamI_mc)

     !--- emission energy fraction per MoCafe bin: j_lam ~ lamI/lam, energy ~ j_lam*dlam.
     do il = 1, nl_mc
        emis(il) = max(lamI_mc(il), 0.0_wp)/sed_wave(il) * sed_dwave(il)
     enddo
     esum = sum(emis)
     if (esum <= 0.0_wp) cycle

     cell_Lemit(ic)  = Labs
     cell_pdf(:,ic)  = emis(:)/esum
     !--- diagnostic "colour temperature": energy-weighted mean wavelength -> Wien.
     lam_mean        = sum(sed_wave(:)*cell_pdf(:,ic))
     cell_Teq(ic)    = 2897.77_wp/max(lam_mean, 1.0e-30_wp)   ! Wien [K], lam in um
  enddo

  !--- combine partial per-cell results across ranks.
  call MPI_ALLREDUCE(MPI_IN_PLACE, cell_Lemit, ncell_tot,       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, cell_pdf,   nl_mc*ncell_tot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, cell_Teq,   ncell_tot,       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  Labs_total = sum(cell_Lemit)

  if (mpar%p_rank == 0) then
     write(*,'(a,es14.6)') 'total absorbed = emitted L : ', Labs_total
     write(*,'(a,es14.6)') 'input stellar luminosity   : ', par%luminosity
     write(*,'(a,f10.5)')  'reprocessed fraction       : ', Labs_total/par%luminosity
  endif

  deallocate(Jcgs, Jsi, Jsed, lamI_sed, lamI_mc, emis)
  end subroutine compute_dustemis

  !---------------------------------------------------------------
  !--- emergent dust thermal SED: raytrace each emitting cell to each
  !--- observer, attenuate at each wavelength, accumulate a per-observer
  !--- 1-D emergent SED and a 3-D (x,y,lambda) image.  Also writes the
  !--- intrinsic (unattenuated) SED and the per-cell T/L maps.
  subroutine write_dustemis(grid)
  use mpi
  use sed_mod,    only : sed_nlam, sed_wave, sed_dwave, sed_sext
  use iofile_mod
  use utility,    only : get_base_name
  implicit none
  type(grid_type), intent(in) :: grid
  type(io_file_type) :: file
  character(len=192) :: filename
  type(photon_type)  :: p
  real(kind=wp), allocatable :: sed_emergent(:,:), sed_intrinsic(:)
  real(kind=wp), allocatable :: Tmap(:,:,:), Lmap(:,:,:)
  real(kind=wp) :: tau, r2, atten, wcell
  integer :: nl_mc, ic, i, j, k, il, kobs, ierr, status

  nl_mc = sed_nlam
  allocate(sed_emergent(nl_mc, par%nobs), sed_intrinsic(nl_mc))
  sed_emergent(:,:) = 0.0_wp
  sed_intrinsic(:)  = 0.0_wp

  !--- intrinsic (whole-grid) dust SED: sum of cell emission, no attenuation.
  do ic = 1, ncell_tot
     if (cell_Lemit(ic) > 0.0_wp) sed_intrinsic(:) = sed_intrinsic(:) + cell_Lemit(ic)*cell_pdf(:,ic)
  enddo

  !--- emergent SED: attenuate each cell's emission to each observer.
  !--- distribute cells across ranks; raytrace at the reference wavelength and
  !--- scale the optical depth per bin by s_ext (grey-rescaling, as in transport).
  do ic = mpar%p_rank+1, ncell_tot, mpar%nproc
     if (cell_Lemit(ic) <= 0.0_wp) cycle
     k = (ic-1)/(grid%nx*grid%ny) + 1
     j = mod((ic-1)/grid%nx, grid%ny) + 1
     i = mod(ic-1, grid%nx) + 1
     p%x = (grid%xface(i)+grid%xface(i+1))*0.5_wp
     p%y = (grid%yface(j)+grid%yface(j+1))*0.5_wp
     p%z = (grid%zface(k)+grid%zface(k+1))*0.5_wp
     p%icell = i;  p%jcell = j;  p%kcell = k
     do kobs = 1, par%nobs
        p%kx = observer(kobs)%x - p%x
        p%ky = observer(kobs)%y - p%y
        p%kz = observer(kobs)%z - p%z
        r2   = p%kx**2 + p%ky**2 + p%kz**2
        p%kx = p%kx/sqrt(r2);  p%ky = p%ky/sqrt(r2);  p%kz = p%kz/sqrt(r2)
        call raytrace_to_edge(p, grid, tau)
        wcell = cell_Lemit(ic)/(fourpi*r2*par%distance2cm**2)
        do il = 1, nl_mc
           atten = exp(-sed_sext(il)*tau)
           sed_emergent(il,kobs) = sed_emergent(il,kobs) + wcell*cell_pdf(il,ic)*atten
        enddo
     enddo
  enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE, sed_emergent, nl_mc*par%nobs, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

  if (mpar%p_rank /= 0) then
     deallocate(sed_emergent, sed_intrinsic)
     return
  endif

  !--- per-cell T and L maps (3-D).
  allocate(Tmap(grid%nx,grid%ny,grid%nz), Lmap(grid%nx,grid%ny,grid%nz))
  do ic = 1, ncell_tot
     k = (ic-1)/(grid%nx*grid%ny) + 1
     j = mod((ic-1)/grid%nx, grid%ny) + 1
     i = mod(ic-1, grid%nx) + 1
     Tmap(i,j,k) = cell_Teq(ic)
     Lmap(i,j,k) = cell_Lemit(ic)
  enddo

  status = 0
  filename = trim(get_base_name(par%out_file))//'_dustsed'//trim(io_file_extension(par%file_format))
  call io_open_new(file, trim(filename), status)
  call io_append_image(file, sed_emergent, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','SED_emergent','emergent dust SED L_lam-like (nlam,nobs)',status)
  call io_put_keyword(file,'SED_UNIT','luminosity/(4pi dist_cm^2)/bin','emergent SED unit',status)
  call io_put_keyword(file,'TOT_LUM', par%luminosity, 'input stellar luminosity',status)
  call io_put_keyword(file,'L_ABS',   Labs_total,     'total absorbed = emitted L',status)
  call io_put_keyword(file,'DIST_CM', par%distance2cm,'distance unit (cm)',status)
  call io_append_image(file, sed_intrinsic, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','SED_intrinsic','intrinsic (unattenuated) dust SED [erg/s per bin]',status)
  call io_append_image(file, sed_wave, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','Wavelength','bin centers [um]',status)
  call io_append_image(file, sed_dwave, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','Dwavelength','bin widths [um]',status)
  call io_append_image(file, Tmap, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','Tdust','per-cell dust colour temperature [K] (Wien, diagnostic)',status)
  call io_append_image(file, Lmap, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','Ldust','per-cell emitted luminosity [erg/s]',status)
  call io_close(file, status)
  write(*,'(2a)') 'dust SED written to: ', trim(filename)

  deallocate(sed_emergent, sed_intrinsic, Tmap, Lmap)
  end subroutine write_dustemis

  !---------------------------------------------------------------
  !--- log-log resample of y(xin) onto xout (both ascending), clamped.
  subroutine resample_to_sed(xin, yin, xout, yout)
  implicit none
  real(kind=wp), intent(in)  :: xin(:), yin(:), xout(:)
  real(kind=wp), intent(out) :: yout(:)
  call loglog_interp(xin, yin, xout, yout)
  end subroutine resample_to_sed

  subroutine resample_from_sed(xin, yin, xout, yout)
  implicit none
  real(kind=wp), intent(in)  :: xin(:), yin(:), xout(:)
  real(kind=wp), intent(out) :: yout(:)
  call loglog_interp(xin, yin, xout, yout)
  end subroutine resample_from_sed

  subroutine loglog_interp(x, y, xnew, ynew)
  implicit none
  real(kind=wp), intent(in)  :: x(:), y(:), xnew(:)
  real(kind=wp), intent(out) :: ynew(:)
  integer :: n, m, i, iq
  real(kind=wp) :: lx, t, ly
  n = size(x);  m = size(xnew)
  do iq = 1, m
     lx = log(xnew(iq))
     if (xnew(iq) <= x(1)) then
        ynew(iq) = max(y(1), 0.0_wp)
     else if (xnew(iq) >= x(n)) then
        ynew(iq) = max(y(n), 0.0_wp)
     else
        do i = 2, n
           if (xnew(iq) <= x(i)) exit
        enddo
        t = (lx - log(x(i-1)))/(log(x(i)) - log(x(i-1)))
        if (y(i-1) > 0.0_wp .and. y(i) > 0.0_wp) then
           ly = log(y(i-1)) + t*(log(y(i)) - log(y(i-1)))
           ynew(iq) = exp(ly)
        else
           ynew(iq) = max((1.0_wp-t)*y(i-1) + t*y(i), 0.0_wp)
        endif
     endif
  enddo
  end subroutine loglog_interp

end module dustemis_mod
