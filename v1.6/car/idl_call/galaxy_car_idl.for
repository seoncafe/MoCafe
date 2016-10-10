c234567
      subroutine galaxy_car_idl(no_photons,nprint,hgg,albedo,luminosity,
     &           iseed4,
     &  dust1,tauface1,dust_rscale1,dust_zscale1,dust_rmax1,dust_zmax1,
     &  dust2,tauface2,dust_rscale2,dust_zscale2,dust_rmax2,dust_zmax2,
     &  xmax,ymax,zmax,nx,ny,nz,
     &  source_name,source_rscale,source_zscale,source_rmax,source_zmax,
     &  obs_angles,obs_distance,nxim,nyim,dxim,dyim,im_scatt,im_direc)

c----- Parameter declarations ----------------------------------------
      implicit none
c input/output variables
      real*8 no_photons,hgg,albedo,luminosity
      integer nprint
      character(len=100) dust1,dust2,source_name
      real*8  tauface1,dust_rscale1,dust_zscale1,dust_rmax1,dust_zmax1,
     &        tauface2,dust_rscale2,dust_zscale2,dust_rmax2,dust_zmax2,
     &        source_rscale,source_zscale,source_rmax,source_zmax
      real*8  xmax,ymax,zmax
      integer  nx,ny,nz,nxim,nyim
!      real*8  im_scatt(nxim,nyim),im_direc(nxim,nyim)
      real    im_scatt(nxim,nyim),im_direc(nxim,nyim)
      real*8  obs_angles(3),obs_distance,dxim,dyim
      integer iseed4

      print*,'test'
      end
