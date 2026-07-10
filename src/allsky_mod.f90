module allsky_mod
!--- HEALPix all-sky map for an interior observer (MoCafe v2.00): the Milky-Way
!--- case, where the observer sits inside the dusty medium and sees light
!--- arriving from every direction.  Instead of a tangent-plane image, each
!--- scatter / dust-reemission / direct event peels toward the observer point
!--- (par%allsky_x/y/z) and is binned into the HEALPix pixel (RING scheme,
!--- nside = par%allsky_nside) of the sky direction it arrives from -- a
!--- panchromatic all-sky intensity map (e.g. the diffuse Galactic light and
!--- the FIR cirrus as seen from the Sun's position).  HEALPix (equal-area
!--- pixels, npix = 12 nside^2) follows LaRT's observer_heal implementation.
!---
!--- Activated by par%allsky = .true.  Surface brightness per wavelength bin.
  use define
  use utility, only : get_base_name
  use healpix, only : nside2npix, vec2pix, pix2vec
  implicit none
  private
  public :: allsky_on, setup_allsky, allsky_peel, allsky_reduce, allsky_write

  logical :: allsky_on = .false.
  integer :: nside = 0, npix = 0, nl = 0
  real(kind=wp) :: obsx, obsy, obsz, steradian_pix
  real(kind=wp), pointer :: sky(:,:) => null()     ! (npix, nlambda) surface brightness

contains
  !---------------------------------------------------------------
  subroutine setup_allsky()
  use sed_mod,    only : sed_nlam
  use memory_mod, only : create_mem
  use mpi
  implicit none
  integer :: ierr
  nl    = sed_nlam
  nside = par%allsky_nside
  npix  = nside2npix(nside)
  if (npix <= 0) then
     if (mpar%p_rank == 0) write(*,'(a,i0,a)') &
        'ERROR: par%allsky_nside = ', nside, ' is not a valid HEALPix nside (power of 2).'
     call MPI_FINALIZE(ierr);  stop
  endif
  obsx = par%allsky_x;  obsy = par%allsky_y;  obsz = par%allsky_z
  steradian_pix = fourpi/dble(npix)      ! HEALPix: equal-area pixels
  call create_mem(sky, [npix, nl])
  sky(:,:) = 0.0_wp
  allsky_on = .true.
  if (mpar%p_rank == 0) then
     write(*,'(a)')      '--- HEALPix all-sky interior observer ---'
     write(*,'(a,3es11.3)') 'observer position     : ', obsx, obsy, obsz
     write(*,'(a,i0,a,i0,a)') 'HEALPix map: nside=', nside, ', npix=', npix, ' (RING)'
  endif
  end subroutine setup_allsky

  !---------------------------------------------------------------
  !--- peel an event toward the interior observer.  kind: 'd' direct (isotropic
  !--- point-source birth), 's' scattered (HG at photon%hgg), 'e' reemission
  !--- (isotropic).  The contribution is a surface brightness binned into the
  !--- HEALPix pixel of the sky direction the light comes from (observer->event).
  subroutine allsky_peel(photon, grid, kind)
  use scan_mod, only : hg_kernel
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  character(len=1),  intent(in) :: kind
  type(photon_type) :: pobs
  real(kind=wp) :: dx,dy,dz,r2,r,tau,cosa,peel,wgt
  integer :: ipix

  if (.not. allsky_on) return
  dx = obsx - photon%x;  dy = obsy - photon%y;  dz = obsz - photon%z
  r2 = dx*dx + dy*dy + dz*dz
  if (r2 <= 0.0_wp) return
  r  = sqrt(r2)
  pobs = photon
  pobs%kx = dx/r;  pobs%ky = dy/r;  pobs%kz = dz/r

  !--- optical depth from the event to the interior observer: integrate only
  !--- over the segment of length r (photon-wavelength units via s_ext).
  tau = raytrace_tau_segment(photon%x, photon%y, photon%z, pobs%kx, pobs%ky, pobs%kz, r, grid) * photon%s_ext

  select case (kind)
  case ('s')
     cosa = photon%kx*pobs%kx + photon%ky*pobs%ky + photon%kz*pobs%kz
     peel = hg_kernel(cosa, photon%hgg)/fourpi
  case default    ! 'd' or 'e': isotropic
     peel = 1.0_wp/fourpi
  end select

  !--- HEALPix pixel of the incoming sky direction (observer -> event = -pobs%k).
  call vec2pix(nside, -pobs%kx, -pobs%ky, -pobs%kz, ipix)
  if (ipix < 1 .or. ipix > npix) return

  !--- accumulate the peel flux (same per-event quantity as the SED observer
  !--- images); allsky_write applies the surface-brightness normalization
  !--- no_photons * steradian_pix * dist_cm^2 / luminosity (as output_normalize).
  wgt = peel/r2 * exp(-tau) * photon%wgt
  sky(ipix,photon%il) = sky(ipix,photon%il) + wgt
  end subroutine allsky_peel

  !---------------------------------------------------------------
  !--- reference-wavelength optical depth along a straight segment of length
  !--- seg from (x0,y0,z0) in direction (kx,ky,kz), on the Cartesian grid
  !--- (Amanatides & Woo cell walk, capped at the segment length).
  function raytrace_tau_segment(x0, y0, z0, kx, ky, kz, seg, grid) result(tau)
  use octree_mod, only : amr_grid, amr_find_leaf, amr_cell_exit, amr_next_leaf
  implicit none
  real(kind=wp), intent(in) :: x0, y0, z0, kx, ky, kz, seg
  type(grid_type), intent(in) :: grid
  real(kind=wp) :: tau, xp, yp, zp, d, tnext, tx, ty, tz, delx, dely, delz, step
  integer :: icell, jcell, kcell, istep, jstep, kstep, idx

  !--- AMR: walk the octree from (x0,y0,z0) accumulating rhokap*step, capped at
  !--- the segment length seg (event -> interior observer).
  if (trim(par%grid_type) == 'amr') then
     block
       integer  :: il, ilcell, iface, il_new
       real(kind=wp) :: t_exit
       tau = 0.0_wp;  xp = x0;  yp = y0;  zp = z0;  d = 0.0_wp
       il = amr_find_leaf(xp, yp, zp)
       do while (il > 0)
          ilcell = amr_grid%icell_of_leaf(il)
          call amr_cell_exit(ilcell, xp, yp, zp, kx, ky, kz, t_exit, iface)
          step = min(d+t_exit, seg) - d
          if (step > 0.0_wp) tau = tau + step*amr_grid%rhokap(il)
          if (d+t_exit >= seg) exit
          d = d + t_exit
          xp = xp + t_exit*kx;  yp = yp + t_exit*ky;  zp = zp + t_exit*kz
          il_new = amr_next_leaf(ilcell, iface, xp, yp, zp)
          if (il_new <= 0) exit
          il = il_new
       enddo
     end block
     return
  endif

  tau = 0.0_wp
  xp = x0;  yp = y0;  zp = z0
  icell = min(max(floor((xp-grid%xmin)/grid%dx)+1, 1), grid%nx)
  jcell = min(max(floor((yp-grid%ymin)/grid%dy)+1, 1), grid%ny)
  kcell = min(max(floor((zp-grid%zmin)/grid%dz)+1, 1), grid%nz)
  if (kx > 0.0_wp) then; istep=1; tx=(grid%xface(icell+1)-xp)/kx; delx=grid%dx/kx
  else if (kx < 0.0_wp) then; istep=-1; tx=(grid%xface(icell)-xp)/kx; delx=-grid%dx/kx
  else; istep=0; tx=hugest; delx=hugest; endif
  if (ky > 0.0_wp) then; jstep=1; ty=(grid%yface(jcell+1)-yp)/ky; dely=grid%dy/ky
  else if (ky < 0.0_wp) then; jstep=-1; ty=(grid%yface(jcell)-yp)/ky; dely=-grid%dy/ky
  else; jstep=0; ty=hugest; dely=hugest; endif
  if (kz > 0.0_wp) then; kstep=1; tz=(grid%zface(kcell+1)-zp)/kz; delz=grid%dz/kz
  else if (kz < 0.0_wp) then; kstep=-1; tz=(grid%zface(kcell)-zp)/kz; delz=-grid%dz/kz
  else; kstep=0; tz=hugest; delz=hugest; endif

  d = 0.0_wp
  do
     idx = minloc([tx,ty,tz], dim=1)
     if (idx == 1) then; tnext = tx; else if (idx == 2) then; tnext = ty; else; tnext = tz; endif
     step = min(tnext, seg) - d
     if (step > 0.0_wp) tau = tau + step*grid%rhokap(icell,jcell,kcell)
     if (tnext >= seg) exit
     d = tnext
     if (idx == 1) then; icell = icell+istep; tx = tx+delx
     else if (idx == 2) then; jcell = jcell+jstep; ty = ty+dely
     else; kcell = kcell+kstep; tz = tz+delz; endif
     if (icell < 1 .or. icell > grid%nx .or. jcell < 1 .or. jcell > grid%ny .or. &
         kcell < 1 .or. kcell > grid%nz) exit
  enddo
  end function raytrace_tau_segment

  !---------------------------------------------------------------
  subroutine allsky_reduce()
  use memory_mod, only : reduce_mem
  implicit none
  if (.not. allsky_on) return
  call reduce_mem(sky)
  end subroutine allsky_reduce

  !---------------------------------------------------------------
  subroutine allsky_write()
  use sed_mod, only : sed_wave, sed_dwave
  use iofile_mod
  implicit none
  type(io_file_type) :: file
  character(len=192) :: filename
  real(kind=wp), allocatable :: pixdir(:,:)
  real(kind=wp) :: vx, vy, vz
  integer :: status, ip
  if (.not. allsky_on .or. mpar%p_rank /= 0) return

  !--- surface-brightness normalization, matching the SED observer images
  !--- (output_normalize): divide by no_photons * steradian_pix * dist_cm^2 /
  !--- luminosity.  Result is luminosity/(dist*dist_cm)^2/sr per wavelength bin.
  block
    real(kind=wp) :: scale_factor
    scale_factor = par%no_photons*steradian_pix*par%distance2cm**2/par%luminosity
    if (scale_factor > 0.0_wp) sky(:,:) = sky(:,:)/scale_factor
  end block

  status = 0
  filename = trim(get_base_name(par%out_file))//'_allsky'//trim(io_file_extension(par%file_format))
  call io_open_new(file, trim(filename), status)
  !--- HEALPix surface-brightness map: AllSky(npix, nlambda), RING ordering.
  call io_append_image(file, sky, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','AllSky','HEALPix all-sky intensity (npix,lambda)',status)
  call io_put_keyword(file,'PIXTYPE','HEALPIX','HEALPix pixelization',status)
  call io_put_keyword(file,'ORDERING','RING','HEALPix ordering scheme',status)
  call io_put_keyword(file,'NSIDE', nside, 'HEALPix resolution parameter',status)
  call io_put_keyword(file,'NPIX',  npix,  'number of pixels (12 nside^2)',status)
  call io_put_keyword(file,'INDXSCHM','IMPLICIT','pixel index scheme',status)
  call io_put_keyword(file,'FIRSTPIX', 1, 'first pixel index (1-based)',status)
  call io_put_keyword(file,'OMEGAPIX', steradian_pix, 'pixel solid angle [sr]',status)
  call io_put_keyword(file,'OBSX', obsx, 'observer X',status)
  call io_put_keyword(file,'OBSY', obsy, 'observer Y',status)
  call io_put_keyword(file,'OBSZ', obsz, 'observer Z',status)
  call io_put_keyword(file,'SB_UNIT','luminosity/dist_cm^2/sr/bin','surface brightness',status)
  !--- unit direction vector of each pixel center (npix x 3), for plotting.
  allocate(pixdir(npix,3))
  do ip = 1, npix
     call pix2vec(nside, ip, vx, vy, vz)
     pixdir(ip,1) = vx;  pixdir(ip,2) = vy;  pixdir(ip,3) = vz
  enddo
  call io_append_image(file, pixdir, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','PixelDir','pixel-center unit vectors (npix,3)',status)
  call io_append_image(file, sed_wave, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','Wavelength','bin centers [um]',status)
  call io_append_image(file, sed_dwave, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','Dwavelength','bin widths [um]',status)
  call io_close(file, status)
  write(*,'(2a)') 'HEALPix all-sky map written to: ', trim(filename)
  deallocate(pixdir)
  end subroutine allsky_write

end module allsky_mod
