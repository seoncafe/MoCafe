module write_mod
  use define
  use iofile_mod
  use utility
  use scan_mod, only : scan_na, scan_ng, scan_alist, scan_glist, scan_g0, &
                       scan_nt, scan_tlist, scan_s, scan_taumax_ref
  use sed_mod,  only : sed_nlam, sed_wave, sed_dwave, sed_lum, sed_cext, &
                       sed_albedo, sed_hgg, sed_cext_ref
  implicit none
  !---- this should be accesible within this module.
  character(len=128) :: fname_backup
  private :: fname_backup
  private :: write_output_peeling
  public  :: write_output_car
contains
!------------------
!-- 2023-06-19,
!-- 2017-06-28, Kwang-il Seon
!------------------
  subroutine write_output_car(filename,grid)
  use define
  implicit none
  character(len=*),    intent(in) :: filename
  type(grid_type),     intent(in) :: grid
  integer          :: k
  character(len=4) :: filename_end

  do k = 1, par%nobs
     if (par%nobs == 1) then
        filename_end = ''
     else
        write(filename_end,'(a,i3.3)') '_',k
     endif
     call write_output_peeling(trim(filename),grid,observer(k),suffix=trim(filename_end))
  enddo
  end subroutine write_output_car
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine write_output_peeling(filename,grid,obs,suffix)
  use define
  implicit none
  character(len=*),    intent(in) :: filename
  type(grid_type),     intent(in) :: grid
  type(observer_type), intent(inout) :: obs
  character(len=*), optional, intent(in) :: suffix
  !--------------------
  type(io_file_type) :: file
  integer            :: status=0
  character(len=128) :: filename1, filename_end
  real(real64)       :: min_value, min_value_tmp
  character(len=8)   :: ext

  if (present(suffix)) then
     filename_end = trim(suffix)
  else
     filename_end = ''
  endif

  ext = io_file_extension(par%file_format)

  !--- Initialize output file name.
  filename1 = trim(get_base_name(filename))//trim(filename_end)//trim(ext)

  !--- updated on 2024.04.04
  if (par%out_bitpix == -32 .and. .not. par%out_bitpix_force) then
     if (par%use_tau_list) then
        min_value     = minval(obs%scatt_agt, obs%scatt_agt > 0.0)
        min_value_tmp = minval(obs%direc_t,   obs%direc_t   > 0.0)
     else if (par%use_ag_list) then
        min_value     = minval(obs%scatt_ag, obs%scatt_ag > 0.0)
        min_value_tmp = minval(obs%direc,    obs%direc    > 0.0)
     else if (par%use_sed) then
        min_value     = minval(obs%scatt_sed, obs%scatt_sed > 0.0)
        min_value_tmp = minval(obs%direc_sed, obs%direc_sed > 0.0)
     else
        min_value     = minval(obs%scatt, obs%scatt > 0.0)
        min_value_tmp = minval(obs%direc, obs%direc > 0.0)
     endif
     if (min_value_tmp < min_value) min_value = min_value_tmp

     if (min_value < tiny(1.0_real32)) par%out_bitpix = -64
  endif

  !--- open output file for peel-off.
  call io_open_new(file,trim(filename1),status)

  !--- write scattered data
  if (par%use_tau_list .and. par%use_ag_list) then
     !--- 5-D scattered image: (x, y, albedo, hgg, tau)
     call io_append_image(file,obs%scatt_agt,status,bitpix=par%out_bitpix)
     call io_put_keyword(file,'EXTNAME','Scattered','J(x,y,albedo,hgg,tau)',status)
     call write_agt_axes(file)
     call write_common_header(file,obs)
  else if (par%use_tau_list) then
     !--- a,g are singletons (na = ng = 1): collapse to a 3-D image (x, y, tau)
     call io_append_image(file,obs%scatt_agt(:,:,1,1,:),status,bitpix=par%out_bitpix)
     call io_put_keyword(file,'EXTNAME','Scattered','J(x,y,tau) (intensity)',status)
     call write_tau_axis_keys(file, 3)
     call write_common_header(file,obs)
  else if (par%use_ag_list) then
     !--- 4-D scattered image: (x, y, albedo, hgg)
     call io_append_image(file,obs%scatt_ag,status,bitpix=par%out_bitpix)
     call io_put_keyword(file,'EXTNAME','Scattered','J(x,y,albedo,hgg)',status)
     call write_ag_axes(file)
     call write_common_header(file,obs)
  else if (par%use_sed) then
     !--- 3-D wavelength-resolved scattered image: (x, y, lambda)
     call io_append_image(file,obs%scatt_sed,status,bitpix=par%out_bitpix)
     call io_put_keyword(file,'EXTNAME','Scattered','J(x,y,lambda) (intensity)',status)
     call write_wave_axis_keys(file, 3)
     call write_common_header(file,obs)
  else
     call io_append_image(file,obs%scatt,status,bitpix=par%out_bitpix)
     call io_put_keyword(file,'EXTNAME','Scattered','J(x,y) (intensity)',status)
     call write_common_header(file,obs)
  endif

  !--- write direct data
  if (par%use_tau_list) then
     !--- 3-D direct image: (x, y, tau)
     call io_append_image(file,obs%direc_t,status,bitpix=par%out_bitpix)
     call io_put_keyword(file,'EXTNAME','Direct','J(x,y,tau) (intensity)',status)
     call write_tau_axis_keys(file, 3)
     call write_common_header(file,obs)
  else if (par%use_sed) then
     !--- 3-D wavelength-resolved direct image: (x, y, lambda)
     call io_append_image(file,obs%direc_sed,status,bitpix=par%out_bitpix)
     call io_put_keyword(file,'EXTNAME','Direct','J(x,y,lambda) (intensity)',status)
     call write_wave_axis_keys(file, 3)
     call write_common_header(file,obs)
  else
     call io_append_image(file,obs%direc,status,bitpix=par%out_bitpix)
     call io_put_keyword(file,'EXTNAME','Direct','J(x,y) (intensity)',status)
     call write_common_header(file,obs)
  endif

  if (par%save_direc0) then
     if (par%use_sed) then
        call io_append_image(file,obs%direc0_sed,status,bitpix=par%out_bitpix)
        call io_put_keyword(file,'EXTNAME','Direct0','J(x,y,lambda) (intensity)',status)
        call write_wave_axis_keys(file, 3)
        call write_common_header(file,obs)
     else
        call io_append_image(file,obs%direc0,status,bitpix=par%out_bitpix)
        call io_put_keyword(file,'EXTNAME','Direct0','J(x,y) (intensity)',status)
        call write_common_header(file,obs)
     endif
  endif

  !--- SED mode: write the wavelength grid and per-bin dust/source tables as
  !--- 1-D datasets, so the analysis side does not depend on header keywords.
  if (par%use_sed) then
     call io_append_image(file,sed_wave,status,bitpix=-64)
     call io_put_keyword(file,'EXTNAME','Wavelength','bin centers [um]',status)
     call io_append_image(file,sed_dwave,status,bitpix=-64)
     call io_put_keyword(file,'EXTNAME','Dwavelength','bin widths [um]',status)
     call io_append_image(file,sed_cext,status,bitpix=-64)
     call io_put_keyword(file,'EXTNAME','Cext','C_ext/H [cm^2/H] per bin',status)
     call io_append_image(file,sed_albedo,status,bitpix=-64)
     call io_put_keyword(file,'EXTNAME','Albedo','dust albedo per bin',status)
     call io_append_image(file,sed_hgg,status,bitpix=-64)
     call io_put_keyword(file,'EXTNAME','Hgg','asymmetry factor g per bin',status)
     call io_append_image(file,sed_lum,status,bitpix=-64)
     call io_put_keyword(file,'EXTNAME','Source_lum','source luminosity fraction per bin',status)
  endif

  !--- close output file
  call io_close(file,status)

  if (par%use_stokes) then
     !--- open new output file for Stokes parameters.
     filename1 = trim(get_base_name(filename))//'_stokes'//trim(filename_end)//trim(ext)

     !--- updated on 2024.04.04
     if (par%out_bitpix == -32) then
        min_value     = minval(obs%I, obs%I > 0.0)
        min_value_tmp = minval(obs%Q, obs%Q > 0.0)
        if (min_value_tmp < min_value) min_value = min_value_tmp
        min_value_tmp = minval(obs%U, obs%U > 0.0)
        if (min_value_tmp < min_value) min_value = min_value_tmp

        if (min_value < tiny(1.0_real32)) par%out_bitpix = -64
     endif

     call io_open_new(file,trim(filename1),status)

     !--- write Stokes I, Q, U, V data
     call io_append_image(file,obs%I,status,bitpix=par%out_bitpix)
     call io_put_keyword(file,'EXTNAME','Stokes_I','Stokes I image',status)
     call write_common_header(file,obs)
     call io_append_image(file,obs%Q,status,bitpix=par%out_bitpix)
     call io_put_keyword(file,'EXTNAME','Stokes_Q','Stokes Q image',status)
     call write_common_header(file,obs)
     call io_append_image(file,obs%U,status,bitpix=par%out_bitpix)
     call io_put_keyword(file,'EXTNAME','Stokes_U','Stokes U image',status)
     call write_common_header(file,obs)
     call io_append_image(file,obs%V,status,bitpix=par%out_bitpix)
     call io_put_keyword(file,'EXTNAME','Stokes_V','Stokes V image',status)
     call write_common_header(file,obs)

     call io_close(file,status)
  endif
  end subroutine write_output_peeling
  !-----------------------------------------------------------------
  subroutine write_common_header(file,obs)
     implicit none
     type(io_file_type),  intent(inout) :: file
     type(observer_type), intent(in)    :: obs
     real(real64)        :: cd1_1, cd1_2, cd2_1, cd2_2
     real(real64)        :: crpix1, crpix2, crval1, crval2
     integer             :: equinox = 2000
     integer             :: status=0
     cd1_1  = par%dxim
     cd1_2  = 0.0_wp
     cd2_1  = 0.0_wp
     cd2_2  = par%dyim
     crpix1 = (par%nxim+1)/2.0_wp
     crpix2 = (par%nyim+1)/2.0_wp
     crval1 = 0.0_wp
     crval2 = 0.0_wp
     call io_put_keyword(file,'EQUINOX',equinox,   'Equinox of Ref. Coord.' ,status)
     call io_put_keyword(file,'CD1_1'  ,cd1_1  ,   'Degree / Pixel',status)
     call io_put_keyword(file,'CD2_1'  ,cd2_1  ,   'Degree / Pixel',status)
     call io_put_keyword(file,'CD1_2'  ,cd1_2  ,   'Degree / Pixel',status)
     call io_put_keyword(file,'CD2_2'  ,cd2_2  ,   'Degree / Pixel' ,status)
     call io_put_keyword(file,'CRPIX1' ,crpix1 ,   'Reference Pixel in X',status)
     call io_put_keyword(file,'CRPIX2' ,crpix2 ,   'Reference Pixel in Y',status)
     call io_put_keyword(file,'CRVAL1' ,crval1 ,   'R.A. (Degree)',status)
     call io_put_keyword(file,'CRVAL2' ,crval2 ,   'Dec  (Degree)',status)
     call io_put_keyword(file,'CTYPE1' ,'RA--TAN', 'Coordinate Type',status)
     call io_put_keyword(file,'CTYPE2' ,'DEC-TAN', 'Coordinate Type',status)
     call io_put_keyword(file,'XMAX',     par%xmax,         'xmax',status)
     call io_put_keyword(file,'YMAX',     par%ymax,         'ymax',status)
     call io_put_keyword(file,'ZMAX',     par%zmax,         'zmax',status)
     call io_put_keyword(file,'DISTANCE', par%distance,     'Distance',status)
     call io_put_keyword(file,'DISTUNIT', par%distance_unit,'Distance Unit',status)
     call io_put_keyword(file,'DIST_CM',  par%distance2cm,  'Distance Unit (cm)',status)
     call io_put_keyword(file,'nphotons', par%no_photons,   'number of photons',status)
     call io_put_keyword(file,'alpha',    obs%alpha, 'Observer alpha (degree)',status)
     call io_put_keyword(file,'beta',     obs%beta,  'Observer beta (degree)',status)
     call io_put_keyword(file,'gamma',    obs%gamma, 'Observer gamma (degree)',status)
     call io_put_keyword(file,'obsx',     obs%x,     'Observer X coordinate',status)
     call io_put_keyword(file,'obsy',     obs%y,     'Observer Y coordinate',status)
     call io_put_keyword(file,'obsz',     obs%z,     'Observer Z coordinate',status)
     call io_put_keyword(file,'taumax',   par%taumax, 'tau_max (z-direction)', status)
     call io_put_keyword(file,'tauhomo',  par%tauhomo,'tau_homo (z-direction)',status)
     call io_put_keyword(file,'TOT_LUM',  par%luminosity, 'TOTAL luminosity' ,status)
     call io_put_keyword(file,'src_geom', trim(par%source_geometry),'source geometry',status)
     call io_put_keyword(file,'SB_value', 'luminosity/(distance*dist_cm)^2/sr', 'surface brightness value',status)
     call io_put_keyword(file,'SB_unit',  'luminosity unit/distance unit^2/sr', 'surface brightness unit' ,status)
     call io_put_keyword(file,'out_norm', trim(par%output_normalization), 'output normalization' ,status)
  end subroutine write_common_header
  !-----------------------------------------------------------------
  subroutine write_ag_axes(file)
  !--- (albedo, asymmetry-factor) axis metadata for the 4-D Scattered image.
  !--- Axis 3 = albedo, axis 4 = hgg.  WCS keywords (CRVAL/CDELT) describe the
  !--- grid when it is uniform; the explicit AVALnnn / GVALnnn lists are
  !--- authoritative and recover an arbitrary (possibly non-uniform) grid.
  implicit none
  type(io_file_type), intent(inout) :: file
  integer          :: status, ia, ig
  character(len=8) :: key
  real(real64)     :: da, dg

  status = 0
  call io_put_keyword(file,'AG_NA', scan_na, 'number of albedo values (axis 3)', status)
  call io_put_keyword(file,'AG_NG', scan_ng, 'number of hgg values (axis 4)',    status)
  call io_put_keyword(file,'AG_G0', scan_g0, 'simulated asymmetry factor g0',    status)

  da = 0.0_real64
  dg = 0.0_real64
  if (scan_na > 1) da = scan_alist(2) - scan_alist(1)
  if (scan_ng > 1) dg = scan_glist(2) - scan_glist(1)
  call io_put_keyword(file,'CTYPE3','ALBEDO',   'axis 3 = dust albedo',   status)
  call io_put_keyword(file,'CRPIX3', 1.0_real64,'reference pixel (axis 3)',status)
  call io_put_keyword(file,'CRVAL3', scan_alist(1),'albedo at CRPIX3',     status)
  call io_put_keyword(file,'CDELT3', da,         'albedo step',            status)
  call io_put_keyword(file,'CTYPE4','HGG',       'axis 4 = asymmetry factor g', status)
  call io_put_keyword(file,'CRPIX4', 1.0_real64,'reference pixel (axis 4)',status)
  call io_put_keyword(file,'CRVAL4', scan_glist(1),'hgg at CRPIX4',        status)
  call io_put_keyword(file,'CDELT4', dg,         'hgg step',               status)

  do ia = 1, scan_na
     write(key,'(a,i3.3)') 'AVAL', ia
     call io_put_keyword(file, key, scan_alist(ia), 'albedo value', status)
  enddo
  do ig = 1, scan_ng
     write(key,'(a,i3.3)') 'GVAL', ig
     call io_put_keyword(file, key, scan_glist(ig), 'hgg value', status)
  enddo
  end subroutine write_ag_axes
  !-----------------------------------------------------------------
  subroutine write_agt_axes(file)
  !--- Axis metadata for the 5-D Scattered image: axis 3 = albedo, axis 4 = hgg
  !--- (written by write_ag_axes), axis 5 = tau (target taumax).
  implicit none
  type(io_file_type), intent(inout) :: file
  call write_ag_axes(file)             ! albedo (axis 3), hgg (axis 4)
  call write_tau_axis_keys(file, 5)    ! tau (axis 5)
  end subroutine write_agt_axes
  !-----------------------------------------------------------------
  subroutine write_tau_axis_keys(file, naxis)
  !--- Polychromatic (tau) axis metadata on FITS/HDF5 axis `naxis`.  The explicit
  !--- TVALnnn list is authoritative (recovers a non-uniform tau grid); CRVAL/CDELT
  !--- describe a uniform grid.  AG_TAU0 is the simulated reference taumax.
  implicit none
  type(io_file_type), intent(inout) :: file
  integer, intent(in) :: naxis
  integer          :: status, it
  character(len=8) :: key
  character(len=1) :: ax
  real(real64)     :: dt

  status = 0
  write(ax,'(i1)') naxis
  call io_put_keyword(file,'AG_NT',   scan_nt,         'number of tau values',         status)
  call io_put_keyword(file,'AG_TAU0', scan_taumax_ref, 'reference (simulated) taumax', status)

  dt = 0.0_real64
  if (scan_nt > 1) dt = scan_tlist(2) - scan_tlist(1)
  call io_put_keyword(file,'CTYPE'//ax,'TAU',          'tau axis (target taumax)',  status)
  call io_put_keyword(file,'CRPIX'//ax, 1.0_real64,    'reference pixel (tau axis)',status)
  call io_put_keyword(file,'CRVAL'//ax, scan_tlist(1), 'tau at reference pixel',    status)
  call io_put_keyword(file,'CDELT'//ax, dt,            'tau step',                  status)

  do it = 1, scan_nt
     write(key,'(a,i3.3)') 'TVAL', it
     call io_put_keyword(file, key, scan_tlist(it), 'tau (target taumax) value', status)
  enddo
  end subroutine write_tau_axis_keys
  !-----------------------------------------------------------------
  subroutine write_wave_axis_keys(file, naxis)
  !--- Wavelength axis metadata on FITS/HDF5 axis `naxis` (SED mode).  The
  !--- grid is log-spaced; CRVAL/CDELT describe it in log10(lambda/um).  The
  !--- explicit per-bin values are stored in the 1-D 'Wavelength' dataset.
  implicit none
  type(io_file_type), intent(inout) :: file
  integer, intent(in) :: naxis
  integer          :: status
  character(len=1) :: ax

  status = 0
  write(ax,'(i1)') naxis
  call io_put_keyword(file,'SED_NLAM', sed_nlam,       'number of wavelength bins',      status)
  call io_put_keyword(file,'SED_LMIN', par%lambda_min, 'min wavelength [um] (bin edge)', status)
  call io_put_keyword(file,'SED_LMAX', par%lambda_max, 'max wavelength [um] (bin edge)', status)
  call io_put_keyword(file,'SED_LREF', par%lambda_ref, 'reference wavelength [um]',      status)
  call io_put_keyword(file,'SED_CREF', sed_cext_ref,   'C_ext/H at lambda_ref [cm^2/H]', status)
  call io_put_keyword(file,'CTYPE'//ax,'WAVE-LOG',     'wavelength axis (log-spaced)',   status)
  call io_put_keyword(file,'CRPIX'//ax, 1.0_real64,    'reference pixel (wave axis)',    status)
  call io_put_keyword(file,'CRVAL'//ax, log10(sed_wave(1)), 'log10(lambda/um) at CRPIX', status)
  call io_put_keyword(file,'CDELT'//ax, log10(sed_wave(2)/sed_wave(1)), 'log10 step',    status)
  end subroutine write_wave_axis_keys
  !-----------------------------------------------------------------
end module write_mod
