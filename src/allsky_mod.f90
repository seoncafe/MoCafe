module allsky_mod
!--- All-sky map for an interior observer (MoCafe v2.00, Stage 7): the
!--- Milky-Way case, where the observer sits inside the dusty medium and sees
!--- light arriving from every direction.  Instead of a tangent-plane image,
!--- each scatter / dust-reemission / direct event peels toward the observer
!--- point (par%obsx/y/z) and is binned by the sky direction it arrives from
!--- into an equirectangular (longitude, latitude) x wavelength map -- a
!--- panchromatic all-sky intensity cube (e.g. the diffuse Galactic light and
!--- the FIR cirrus as seen from the Sun's position).
!---
!--- Activated by par%allsky = .true. with the observer position given by
!--- par%obsx/y/z (interior).  Surface brightness per wavelength bin.
  use define
  use utility, only : get_base_name
  implicit none
  private
  public :: allsky_on, setup_allsky, allsky_peel, allsky_reduce, allsky_write

  logical :: allsky_on = .false.
  integer :: nlon = 0, nlat = 0, nl = 0
  real(kind=wp) :: obsx, obsy, obsz, dlon, dlat
  real(kind=wp), pointer :: sky(:,:,:) => null()   ! (nlon, nlat, nlambda) intensity

contains
  !---------------------------------------------------------------
  subroutine setup_allsky()
  use sed_mod,    only : sed_nlam
  use memory_mod, only : create_mem
  implicit none
  nl   = sed_nlam
  nlon = par%allsky_nlon
  nlat = par%allsky_nlat
  obsx = par%allsky_x;  obsy = par%allsky_y;  obsz = par%allsky_z
  dlon = twopi/dble(nlon)         ! [0,2pi)
  dlat = pi/dble(nlat)            ! [-pi/2, pi/2]
  call create_mem(sky, [nlon, nlat, nl])
  sky(:,:,:) = 0.0_wp
  allsky_on = .true.
  if (mpar%p_rank == 0) then
     write(*,'(a)')      '--- All-sky interior observer (Stage 7) ---'
     write(*,'(a,3es11.3)') 'observer position     : ', obsx, obsy, obsz
     write(*,'(a,i0,a,i0)') 'sky map (lon x lat)   : ', nlon, ' x ', nlat
  endif
  end subroutine setup_allsky

  !---------------------------------------------------------------
  !--- peel an event toward the interior observer.  kind: 'd' direct (isotropic
  !--- point-source birth), 's' scattered (HG at photon%hgg), 'e' reemission
  !--- (isotropic).  The contribution is a surface brightness binned by the sky
  !--- direction the light comes from (observer -> event).
  subroutine allsky_peel(photon, grid, kind)
  use scan_mod, only : hg_kernel
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  character(len=1),  intent(in) :: kind
  type(photon_type) :: pobs
  real(kind=wp) :: dx,dy,dz,r2,r,tau,cosa,peel,wgt,lon,lat,domega
  integer :: ilon, ilat

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

  !--- sky pixel from the direction observer -> event = -pobs%k.
  lon = atan2(-pobs%ky, -pobs%kx);  if (lon < 0.0_wp) lon = lon + twopi
  lat = asin(max(-1.0_wp, min(1.0_wp, -pobs%kz)))
  ilon = floor(lon/dlon) + 1
  ilat = floor((lat + halfpi)/dlat) + 1
  if (ilon < 1) ilon = 1;  if (ilon > nlon) ilon = nlon
  if (ilat < 1) ilat = 1;  if (ilat > nlat) ilat = nlat

  !--- surface brightness: flux / pixel solid angle.
  domega = dlon*dlat*cos(lat)
  if (domega <= 0.0_wp) domega = dlon*dlat
  wgt = peel/r2 * exp(-tau) * photon%wgt / domega
  sky(ilon,ilat,photon%il) = sky(ilon,ilat,photon%il) + wgt
  end subroutine allsky_peel

  !---------------------------------------------------------------
  !--- reference-wavelength optical depth along a straight segment of length
  !--- seg from (x0,y0,z0) in direction (kx,ky,kz), on the Cartesian grid
  !--- (Amanatides & Woo cell walk, capped at the segment length).
  function raytrace_tau_segment(x0, y0, z0, kx, ky, kz, seg, grid) result(tau)
  implicit none
  real(kind=wp), intent(in) :: x0, y0, z0, kx, ky, kz, seg
  type(grid_type), intent(in) :: grid
  real(kind=wp) :: tau, xp, yp, zp, d, tnext, tx, ty, tz, delx, dely, delz, step
  integer :: icell, jcell, kcell, istep, jstep, kstep, idx

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
  integer :: status
  if (.not. allsky_on .or. mpar%p_rank /= 0) return
  status = 0
  filename = trim(get_base_name(par%out_file))//'_allsky'//trim(io_file_extension(par%file_format))
  call io_open_new(file, trim(filename), status)
  call io_append_image(file, sky, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','AllSky','all-sky intensity (lon,lat,lambda)',status)
  call io_put_keyword(file,'CTYPE1','GLON-CAR','longitude axis',status)
  call io_put_keyword(file,'CTYPE2','GLAT-CAR','latitude axis',status)
  call io_put_keyword(file,'CDELT1', dlon*rad2deg, 'deg/pixel (lon)',status)
  call io_put_keyword(file,'CDELT2', dlat*rad2deg, 'deg/pixel (lat)',status)
  call io_put_keyword(file,'OBSX', obsx, 'observer X',status)
  call io_put_keyword(file,'OBSY', obsy, 'observer Y',status)
  call io_put_keyword(file,'OBSZ', obsz, 'observer Z',status)
  call io_put_keyword(file,'SB_UNIT','luminosity/dist_cm^2/sr/bin','surface brightness',status)
  call io_append_image(file, sed_wave, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','Wavelength','bin centers [um]',status)
  call io_append_image(file, sed_dwave, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','Dwavelength','bin widths [um]',status)
  call io_close(file, status)
  write(*,'(2a)') 'all-sky map written to: ', trim(filename)
  end subroutine allsky_write

end module allsky_mod
