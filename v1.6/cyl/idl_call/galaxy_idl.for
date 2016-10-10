c234567
c  This fortran 77 file is to make interface between fortran 90 routines and IDL.
c
      subroutine galaxy_idl(no_photons,nprint,hgg,albedo,luminosity,
     &  dust1,tauface1,dust_rscale1,dust_zscale1,dust_rmax1,dust_zmax1,
     &  dust2,tauface2,dust_rscale2,dust_zscale2,dust_rmax2,dust_zmax2,
     &  rmax,pmax,zmax,nr,np,nz,
     &  disk_name,disk_rscale,disk_zscale,disk_rmax,disk_zmax,
     &  bulge_name,sersic_index,Reff,axial_ratio,BulgeToDisk,
     &  inclination_angle,position_angle,phase_angle,distance,
     *  nxim,nyim,dxim,dyim,im_scatt,im_direc,
     &  im_scatt_sig,im_direc_sig,
     &  psf_file,left_right_fold)

c----- Parameter declarations ----------------------------------------
      implicit none
c input/output variables
      real*8  no_photons,hgg,albedo,luminosity
      integer nprint
      character(len=128) dust1,dust2,disk_name,bulge_name
      real*8  tauface1,dust_rscale1,dust_zscale1,dust_rmax1,dust_zmax1,
     &        tauface2,dust_rscale2,dust_zscale2,dust_rmax2,dust_zmax2,
     &        disk_rscale,disk_zscale,disk_rmax,disk_zmax,
     &        sersic_index,Reff,axial_ratio,BulgeToDisk
      real*8  rmax,pmax,zmax
      integer nr,np,nz,nxim,nyim
      real    im_scatt(nxim,nyim),im_direc(nxim,nyim),
     &        im_scatt_sig(nxim,nyim),im_direc_sig(nxim,nyim)
      real*8  inclination_angle,position_angle,phase_angle,
     &        distance,dxim,dyim
      character(len=128) psf_file
      integer left_right_fold

      print*,'test'
      end
