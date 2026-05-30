!-- Modification History
!   2021-06-11, now the routines for sightline_tau are gathered in a separate module.
!--
module sightline_tau_mod
  use define
  use utility
  use iofile_mod
  use memory_mod
  use raytrace, only : raytrace_to_edge_car
contains
  !--------------------------------------------------
  subroutine make_sightline_tau(grid)
  !--- calculate the optical depth along sight lines that are projected to a detector-plane pixel (2020/09/20).
  use mpi
  !--- now, write subrotines are defined in define.f90 (2023.01.16).
  !use write_mod
  implicit none
  type(grid_type),  intent(in) :: grid
  !--- local variables
  type(photon_type) :: pobs
  real(kind=wp) :: delt(6),dist
  real(kind=wp) :: tau_dust
  real(kind=wp) :: kx,ky,kz,kr,u1
  integer       :: ix,iy
  integer       :: loop,loop1,loop2,nsize
  integer       :: ierr
  integer       :: i,jj,kk,j0
  character(len=4) :: filename_end

  !--- memories allocation
  do i=1,par%nobs
     call create_shared_mem(observer(i)%tau, [par%nxim,par%nyim])
  enddo

  !--- calculate tau, N(gas) maps.
  do i=1,par%nobs
    nsize = observer(i)%nxim * observer(i)%nyim
    call loop_divide(nsize,mpar%nproc,mpar%p_rank,loop1,loop2)
    do loop=loop1,loop2
       call array_2D_indices(observer(i)%nxim,observer(i)%nyim,loop,ix,iy)
       kx = tan((ix - (observer(i)%nxim+1.0_wp)/2.0_wp) * observer(i)%dxim/rad2deg)
       ky = tan((iy - (observer(i)%nyim+1.0_wp)/2.0_wp) * observer(i)%dyim/rad2deg)
       kz = -1.0_wp
       kr = sqrt(kx*kx + ky*ky + kz*kz)
       kx = kx/kr
       ky = ky/kr
       kz = kz/kr
       pobs%kx = observer(i)%rmatrix(1,1)*kx + observer(i)%rmatrix(2,1)*ky + observer(i)%rmatrix(3,1)*kz
       pobs%ky = observer(i)%rmatrix(1,2)*kx + observer(i)%rmatrix(2,2)*ky + observer(i)%rmatrix(3,2)*kz
       pobs%kz = observer(i)%rmatrix(1,3)*kx + observer(i)%rmatrix(2,3)*ky + observer(i)%rmatrix(3,3)*kz
       if (pobs%kx == 0.0_wp) then
          delt(1) = hugest
          delt(2) = hugest
       else
          delt(1) = (grid%xmax-observer(i)%x)/pobs%kx
          delt(2) = (grid%xmin-observer(i)%x)/pobs%kx
       endif
       if (pobs%ky == 0.0_wp) then
          delt(3) = hugest
          delt(4) = hugest
       else
          delt(3) = (grid%ymax-observer(i)%y)/pobs%ky
          delt(4) = (grid%ymin-observer(i)%y)/pobs%ky
       endif
       if (pobs%kz == 0.0_wp) then
          delt(5) = hugest
          delt(6) = hugest
       else
          delt(5) = (grid%zmax-observer(i)%z)/pobs%kz
          delt(6) = (grid%zmin-observer(i)%z)/pobs%kz
       endif

       !-- Find the farthest boundary where the ray touches the grid system.
       !-- We measure optical depth from the distant universe toward Earth. (2020.10.20)
       dist = -999.9_wp
       do jj=1,6
          if (delt(jj) > 0.0_wp .and. delt(jj) < hugest) then
             !-- Do not delete this part (2023.01.20).
             pobs%x     = observer(i)%x + pobs%kx * delt(jj)
             pobs%y     = observer(i)%y + pobs%ky * delt(jj)
             pobs%z     = observer(i)%z + pobs%kz * delt(jj)
             pobs%icell = floor((pobs%x-grid%xmin)/grid%dx)+1
             pobs%jcell = floor((pobs%y-grid%ymin)/grid%dy)+1
             pobs%kcell = floor((pobs%z-grid%zmin)/grid%dz)+1
             !--- to reduce numerical errors (2021.08.05)
             if (jj == 1) pobs%icell = grid%nx+1
             if (jj == 2) pobs%icell = 1
             if (jj == 3) pobs%jcell = grid%ny+1
             if (jj == 4) pobs%jcell = 1
             if (jj == 5) pobs%kcell = grid%nz+1
             if (jj == 6) pobs%kcell = 1
             if (pobs%icell >=1 .and. pobs%icell <= grid%nx+1 .and. &
                 pobs%jcell >=1 .and. pobs%jcell <= grid%ny+1 .and. &
                 pobs%kcell >=1 .and. pobs%kcell <= grid%nz+1) then
                if (delt(jj) > dist) then
                   dist = delt(jj)
                   j0   = jj
                endif
             endif
          endif
       enddo

       if (dist > 0.0_wp .and. dist < hugest) then
          pobs%x     = observer(i)%x + pobs%kx * dist
          pobs%y     = observer(i)%y + pobs%ky * dist
          pobs%z     = observer(i)%z + pobs%kz * dist
          pobs%icell = floor((pobs%x-grid%xmin)/grid%dx)+1
          pobs%jcell = floor((pobs%y-grid%ymin)/grid%dy)+1
          pobs%kcell = floor((pobs%z-grid%zmin)/grid%dz)+1

          !--- After finding the starting position of the ray, the direction vector should be revsersed. (2020.10.20)
          pobs%kx = -pobs%kx
          pobs%ky = -pobs%ky
          pobs%kz = -pobs%kz

          !--- to reduce numerical errors (2021.08.05)
          if (j0 == 1) then
             pobs%icell = grid%nx+1
             pobs%x     = grid%xface(grid%nx+1)
          else if (j0 == 2) then
             pobs%icell = 1
             pobs%x     = grid%xface(1)
          else if (j0 == 3) then
             pobs%jcell = grid%ny+1
             pobs%y     = grid%yface(grid%ny+1)
          else if (j0 == 4) then
             pobs%jcell = 1
             pobs%y     = grid%yface(1)
          else if (j0 == 5) then
             pobs%kcell = grid%nz+1
             pobs%z     = grid%zface(grid%nz+1)
          else if (j0 == 6) then
             pobs%kcell = 1
             pobs%z     = grid%zface(1)
          endif

          if (pobs%icell == grid%nx+1 .and. pobs%kx < 0.0_wp) pobs%icell = grid%nx
          if (pobs%jcell == grid%ny+1 .and. pobs%ky < 0.0_wp) pobs%jcell = grid%ny
          if (pobs%kcell == grid%nz+1 .and. pobs%kz < 0.0_wp) pobs%kcell = grid%nz

          if (pobs%icell >=1 .and. pobs%icell <= grid%nx .and. &
              pobs%jcell >=1 .and. pobs%jcell <= grid%ny .and. &
              pobs%kcell >=1 .and. pobs%kcell <= grid%nz) then
             !--- note that there should be no limitation on tau (tau_huge) in this routine (2020.11.08).
             call raytrace_to_edge_car(pobs,grid,tau_dust)
             observer(i)%tau(ix,iy) = tau_dust
          endif
       endif
    enddo

    call reduce_mem(observer(i)%tau, shared_memory=.true.)
  enddo

  if (mpar%p_rank == 0) then
     !--- write tau and N(gas) maps
     do i = 1, par%nobs
        if (par%nobs == 1) then
           filename_end = ''
        else
           write(filename_end,'(a,i3.3)') '_',i
        endif
        call write_sightline_tau(trim(par%out_file),grid,observer(i), suffix=trim(filename_end))
     enddo
  endif

  !--- memories deallocation
  do i=1, par%nobs
     call destroy_mem(observer(i)%tau)
  enddo
  end subroutine make_sightline_tau
  !-------------------------------------------------------
  subroutine write_sightline_tau(filename,grid,obs,suffix)
  implicit none
  character(len=*),    intent(in) :: filename
  type(grid_type),     intent(in) :: grid
  type(observer_type), intent(in) :: obs
  character(len=*), optional, intent(in) :: suffix
  !--------------------
  type(io_file_type) :: file
  integer            :: status=0
  character(len=128) :: filename1, filename_end
  real(real64)       :: cd1_1, cd1_2, cd2_1, cd2_2
  real(real64)       :: crpix1, crpix2, crval1, crval2
  integer            :: equinox = 2000

  if (present(suffix)) then
     filename_end = trim(suffix)
  else
     filename_end = ''
  endif

  !--- Initialize output file name.
  filename1 = trim(get_base_name(filename))//'_tau'//trim(filename_end)//trim(io_file_extension(par%file_format))

  call io_open_new(file,trim(filename1),status)

  !--- write Images for dust optical depth
  call io_append_image(file,obs%tau,status,bitpix=par%out_bitpix)
  call io_put_keyword(file,'EXTNAME','TAU_dust','dust optical depth',status)

  !--- write keywords
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
  call io_put_keyword(file,'nphotons', par%no_photons,   'number of photons',status)
  call io_put_keyword(file,'alpha',    obs%alpha, 'Observer alpha (degree)',status)
  call io_put_keyword(file,'beta',     obs%beta,  'Observer beta (degree)',status)
  call io_put_keyword(file,'gamma',    obs%gamma, 'Observer gamma (degree)',status)
  call io_put_keyword(file,'obsx',     obs%x,     'Observer X coordinate',status)
  call io_put_keyword(file,'obsy',     obs%y,     'Observer Y coordinate',status)
  call io_put_keyword(file,'obsz',     obs%z,     'Observer Z coordinate',status)
  call io_put_keyword(file,'taumax',   par%taumax, 'tau_max (z-direction)', status)
  call io_put_keyword(file,'tauhomo',  par%tauhomo,'tau_homo (z-direction)',status)
  call io_put_keyword(file,'lum',      par%luminosity,'luminosity',status)
  call io_put_keyword(file,'src_geom', trim(par%source_geometry),'source geometry',status)

  call io_close(file,status)
  end subroutine write_sightline_tau
  !-------------------------------------------------------
end module sightline_tau_mod
