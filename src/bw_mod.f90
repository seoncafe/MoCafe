module bw_mod
!--- Bjorkman & Wood (2001, ApJ 554, 615) immediate-reemission dust emission
!--- (MoCafe v2.00, Stage 4).  Approximate, equilibrium, large-grain (mixture
!--- mean-opacity) treatment with NO global iteration: a photon packet keeps a
!--- fixed energy and random-walks through scatter and absorb+immediate-reemit
!--- events until it escapes.  Each absorption raises the local dust
!--- temperature (from the running absorbed energy) and the packet is re-emitted
!--- with a wavelength drawn from the temperature-correction spectrum
!--- (proportional to kappa_abs * dB_lambda/dT), so the time-averaged emitted
!--- spectrum equals kappa_abs*B_lambda(T_final) without iterating.
!---
!--- Uses the size-integrated mixture absorption opacity kappa_abs(lambda) =
!--- C_ext*(1-albedo) from the SED extinction table (par%kext_file); PAH /
!--- stochastic-heating features are NOT captured (that is Mode 1, Lucy+SEDust).
  use define
  use utility, only : time_stamp
  use sed_mod, only : sed_nlam, sed_wave, sed_dwave, sed_sext, sed_albedo, sed_cext, sed_cext_ref
  implicit none
  private
  public :: setup_bw, bw_reemit, bw_finalize, bw_Tmap, bw_write, run_bw

  integer :: nl = 0, NT = 0
  real(kind=wp), allocatable :: kabs(:)        ! (nl) C_abs/H [cm^2/H] on the RT grid
  real(kind=wp), allocatable :: Tgrid(:)       ! (NT) temperature grid [K]
  real(kind=wp), allocatable :: kappB(:)       ! (NT) integral C_abs*B_lam(T) dlam [erg/s/sr/H]
  real(kind=wp), allocatable :: cdf_emit(:,:)  ! (nl, NT) cumulative emission spectrum vs T
  real(kind=wp), allocatable :: cell_A(:)      ! (ncell) absorbed power (rank-local) [erg/s]
  real(kind=wp), allocatable :: cell_T(:)      ! (ncell) current dust temperature [K]
  integer :: nx=0, ny=0, nz=0, ncell_tot=0
  real(kind=wp) :: cellfac = 0.0_wp            ! Vcode*distance2cm^2/Cext_ref (N_H = rhokap*cellfac)

  real(kind=wp), parameter :: hc2_cgs = 1.1910429723971884e-5_wp  ! 2 h c^2 [erg cm^2 / s]
  real(kind=wp), parameter :: hck_cgs = 1.4387768775039337_wp     ! h c / k_B [cm K]

contains
  !---------------------------------------------------------------
  subroutine setup_bw(grid)
  use mpi
  implicit none
  type(grid_type), intent(in) :: grid
  real(kind=wp) :: Tlo, Thi, dlnT, T, lam_cm, x, ex, psum
  integer :: il, it

  nl = sed_nlam
  NT = par%sed_NT
  Tlo = par%sed_Tlo;  Thi = par%sed_Thi
  nx = grid%nx;  ny = grid%ny;  nz = grid%nz
  ncell_tot = nx*ny*nz
  cellfac = grid%dx*grid%dy*grid%dz * par%distance2cm**2 / sed_cext_ref

  allocate(kabs(nl), Tgrid(NT), kappB(NT), cdf_emit(nl, NT))
  allocate(cell_A(ncell_tot), cell_T(ncell_tot))
  cell_A(:) = 0.0_wp;  cell_T(:) = Tlo

  !--- mixture absorption cross section per H [cm^2/H].
  kabs(:) = sed_cext(:)*(1.0_wp - sed_albedo(:))

  !--- temperature grid and the kappB(T) = int C_abs B_lam dlam table, plus the
  !--- temperature-correction emission CDF (proportional to C_abs*dB/dT).
  dlnT = log(Thi/Tlo)/dble(NT-1)
  do it = 1, NT
     T = Tlo*exp(dble(it-1)*dlnT)
     Tgrid(it) = T
     kappB(it) = 0.0_wp
     psum = 0.0_wp
     do il = 1, nl
        lam_cm = sed_wave(il)*1.0e-4_wp
        x = hck_cgs/(lam_cm*T)
        if (x < 700.0_wp) then
           ex = exp(x)
           !--- C_abs * B_lam(T) * dlam  [erg/s/sr/H]; B_lam in erg/s/cm^2/sr/cm,
           !--- dlam converted um->cm.
           kappB(it) = kappB(it) + kabs(il)*(hc2_cgs/lam_cm**5/(ex-1.0_wp))*(sed_dwave(il)*1.0e-4_wp)
           !--- emission PDF shape ~ C_abs * dB_lam/dT * dlam
           !--- dB/dT ~ B * x*ex/((ex-1)T); keep only the shape (T-const drops).
           cdf_emit(il,it) = kabs(il)*(hc2_cgs/lam_cm**5/(ex-1.0_wp))*(x*ex/(ex-1.0_wp))*(sed_dwave(il)*1.0e-4_wp)
        else
           cdf_emit(il,it) = 0.0_wp
        endif
        psum = psum + cdf_emit(il,it)
     enddo
     !--- cumulative (normalized) emission distribution at this T.
     if (psum > 0.0_wp) then
        cdf_emit(1,it) = cdf_emit(1,it)/psum
        do il = 2, nl
           cdf_emit(il,it) = cdf_emit(il-1,it) + cdf_emit(il,it)/psum
        enddo
     endif
  enddo

  if (mpar%p_rank == 0) then
     write(*,'(a)')      '--- Dust thermal emission (Bjorkman & Wood 2001, immediate reemission) ---'
     write(*,'(a,i0,a,2f8.1,a)') 'temperature grid: NT=', NT, ', T=[', Tlo, Thi, '] K'
     write(*,'(a)')      'mixture mean opacity (no PAH/stochastic features)'
  endif
  end subroutine setup_bw

  !---------------------------------------------------------------
  !--- absorption + immediate reemission: deposit the packet energy into the
  !--- cell, update the dust temperature from the running absorbed power, and
  !--- resample the packet wavelength from the temperature-correction spectrum.
  !--- The MPI ranks each transport 1/nproc of the photons, so the rank-local
  !--- absorbed power is scaled by nproc to estimate the cell total (standard
  !--- MPI B&W approximation).
  subroutine bw_reemit(photon, grid)
  use random, only : rand_number
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(in)    :: grid
  integer :: ic, it, lo, hi, mid, il
  real(kind=wp) :: N_H, target_kappB, u, cost, sint, phi

  ic = (photon%kcell-1)*nx*ny + (photon%jcell-1)*nx + photon%icell
  cell_A(ic) = cell_A(ic) + photon%Lpacket

  !--- equilibrium temperature: 4 pi N_H kappB(T) = absorbed power.
  N_H = grid%rhokap(photon%icell,photon%jcell,photon%kcell)*cellfac
  if (N_H <= 0.0_wp) then
     it = 1
  else
     target_kappB = cell_A(ic)*dble(mpar%nproc)/(fourpi*N_H)
     !--- invert the monotonic kappB(T) by binary search.
     if (target_kappB <= kappB(1)) then
        it = 1
     else if (target_kappB >= kappB(NT)) then
        it = NT
     else
        lo = 1;  hi = NT
        do while (hi - lo > 1)
           mid = (lo+hi)/2
           if (kappB(mid) < target_kappB) then
              lo = mid
           else
              hi = mid
           endif
        enddo
        it = hi
     endif
  endif
  cell_T(ic) = Tgrid(it)

  !--- sample the reemission wavelength from the temperature-correction CDF.
  u = rand_number()
  lo = 1;  hi = nl
  do while (lo < hi)
     mid = (lo+hi)/2
     if (u <= cdf_emit(mid,it)) then
        hi = mid
     else
        lo = mid+1
     endif
  enddo
  il = lo
  photon%il     = il
  photon%lambda = sed_wave(il)
  photon%s_ext  = sed_sext(il)
  photon%albedo = sed_albedo(il)
  photon%hgg    = 0.0_wp     ! isotropic scattering after reemission is not tracked; g reset

  !--- isotropic reemission direction.
  cost = 2.0_wp*rand_number() - 1.0_wp
  sint = sqrt(1.0_wp - cost*cost)
  phi  = twopi*rand_number()
  photon%kx = sint*cos(phi)
  photon%ky = sint*sin(phi)
  photon%kz = cost
  end subroutine bw_reemit

  !---------------------------------------------------------------
  !--- combine the per-rank absorbed power and finalize cell temperatures.
  subroutine bw_finalize()
  use mpi
  implicit none
  integer :: ic, it, lo, hi, mid, ierr
  real(kind=wp) :: N_H, target_kappB
  !--- (kept for a global temperature map; transport used the nproc-scaled
  !--- rank-local estimate.)  Here we ALLREDUCE the true total absorbed power.
  call MPI_ALLREDUCE(MPI_IN_PLACE, cell_A, ncell_tot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  end subroutine bw_finalize

  !---------------------------------------------------------------
  !--- per-cell dust temperature map from the final total absorbed power.
  subroutine bw_Tmap(grid)
  implicit none
  type(grid_type), intent(in) :: grid
  integer :: ic, i, j, k, it, lo, hi, mid
  real(kind=wp) :: N_H, target_kappB
  do ic = 1, ncell_tot
     k = (ic-1)/(nx*ny) + 1
     j = mod((ic-1)/nx, ny) + 1
     i = mod(ic-1, nx) + 1
     N_H = grid%rhokap(i,j,k)*cellfac
     if (N_H <= 0.0_wp .or. cell_A(ic) <= 0.0_wp) then
        cell_T(ic) = Tgrid(1);  cycle
     endif
     target_kappB = cell_A(ic)/(fourpi*N_H)
     if (target_kappB <= kappB(1)) then
        it = 1
     else if (target_kappB >= kappB(NT)) then
        it = NT
     else
        lo = 1;  hi = NT
        do while (hi-lo > 1)
           mid = (lo+hi)/2
           if (kappB(mid) < target_kappB) then; lo = mid; else; hi = mid; endif
        enddo
        it = hi
     endif
     cell_T(ic) = Tgrid(it)
  enddo
  end subroutine bw_Tmap

  !---------------------------------------------------------------
  !--- Bjorkman & Wood transport: fixed-energy packets random-walk through
  !--- scatter and absorb+immediate-reemit events until they escape, peeling
  !--- to the observers at birth (direct), each scattering (HG), and each
  !--- reemission (isotropic).  No forced first scattering / reduced weight
  !--- (B&W needs discrete interactions).  Stellar packets are split across
  !--- ranks (equal share).
  subroutine run_bw(grid)
  use photon_mod,   only : gen_photon
  use peelingoff_mod, only : peeling_scattered_photon_nostokes_sed, peeling_reemit_sed
  use random,       only : rand_number, rand_henyey_greenstein
  use mpi
  implicit none
  type(grid_type), intent(inout) :: grid
  type(photon_type) :: photon
  integer(kind=int64) :: ip
  integer :: ierr
  real(kind=wp) :: tau, xi, cost, sint, phi, cosp, sinp, kx1, ky1, kz1, kr
  real(kind=wp) :: dtime

  do ip = mpar%p_rank+1, par%nphotons, mpar%nproc
     call gen_photon(grid, photon)     ! sets position/dir/wavelength, Lpacket, direct peel
     do while (photon%inside)
        tau = -log(rand_number())
        call raytrace_to_tau(photon, grid, tau/photon%s_ext)
        if (.not. photon%inside) exit
        xi = rand_number()
        if (xi < photon%albedo) then
           !--- scattering: peel (HG at current wavelength), then new direction.
           call peeling_scattered_photon_nostokes_sed(photon, grid)
           cost = rand_henyey_greenstein(photon%hgg)
           sint = sqrt(1.0_wp - cost*cost)
           phi  = twopi*rand_number();  cosp = cos(phi);  sinp = sin(phi)
           kx1 = photon%kx;  ky1 = photon%ky;  kz1 = photon%kz
           kr  = sqrt(kx1*kx1 + ky1*ky1)
           photon%kx = cost*kx1 + sint*(kz1*kx1*cosp - ky1*sinp)/kr
           photon%ky = cost*ky1 + sint*(kz1*ky1*cosp + kx1*sinp)/kr
           photon%kz = cost*kz1 - sint*cosp*kr
        else
           !--- absorption + immediate reemission (new wavelength, isotropic dir),
           !--- then peel the reemission.
           call bw_reemit(photon, grid)
           call peeling_reemit_sed(photon, grid)
        endif
        photon%nscatt = photon%nscatt + 1
     enddo
     if (mpar%p_rank == 0 .and. mod(ip, int(par%nprint,int64)) == 0) then
        call time_stamp(dtime)
        write(6,'(es14.5,a,f8.3,a)') dble(ip),' photons calculated: ',dtime/60.0_wp,' mins'
     endif
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end subroutine run_bw

  !---------------------------------------------------------------
  subroutine bw_write(grid)
  use iofile_mod
  use utility, only : get_base_name
  implicit none
  type(grid_type), intent(in) :: grid
  type(io_file_type) :: file
  character(len=192) :: filename
  real(kind=wp), allocatable :: Tmap(:,:,:), Amap(:,:,:)
  integer :: ic, i, j, k, status
  if (mpar%p_rank /= 0) return
  allocate(Tmap(nx,ny,nz), Amap(nx,ny,nz))
  do ic = 1, ncell_tot
     k = (ic-1)/(nx*ny) + 1;  j = mod((ic-1)/nx, ny) + 1;  i = mod(ic-1, nx) + 1
     Tmap(i,j,k) = cell_T(ic);  Amap(i,j,k) = cell_A(ic)
  enddo
  status = 0
  filename = trim(get_base_name(par%out_file))//'_bwdust'//trim(io_file_extension(par%file_format))
  call io_open_new(file, trim(filename), status)
  call io_append_image(file, Tmap, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','Tdust','equilibrium dust temperature [K] (B&W 2001)',status)
  call io_append_image(file, Amap, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','Labs','per-cell absorbed = emitted power [erg/s]',status)
  call io_put_keyword(file,'TOT_LUM', par%luminosity, 'input stellar luminosity',status)
  call io_put_keyword(file,'L_ABS',   sum(cell_A),    'total absorbed power',status)
  call io_close(file, status)
  write(*,'(2a)') 'B&W dust temperature/absorption written to: ', trim(filename)
  write(*,'(a,es14.6,a,f8.4)') 'total absorbed L = ', sum(cell_A), '  (fraction ', sum(cell_A)/par%luminosity, ')'
  deallocate(Tmap, Amap)
  end subroutine bw_write

end module bw_mod
