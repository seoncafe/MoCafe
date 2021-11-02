--------------------------------------------------------------------------------
Author:
   Kwang-il Seon, KASI

Version:
   1.02 - 2021/11/01

Purpose:
   Radiative Transfer in cartesian grid system (x, y, z) \
   Basic Monte-Carlo Simulation for a reflection nebula with a point source embedded
   in a uniform constant-density cloud.

   Unit of the output image = (luminosity unit) cm^-2 sr^-1

Howto run:

   1. Edit model.in
      no_photons : number of photons
      hgg : asymmetry factor for the Henyey-Greenstein function
      albedo : dust albedo
      taumax : optical depth from the center to the outer edge
      luminosity : luminosity
      outfile : fit file name for the output image
      xmax, ymax, zmax : dust cloud size (-xmax < x < xmax, etc)
      nx, ny, nz : number of cells in x, y, z coordinates
      iseed : initial random number seed (iseed = 0 to generate a random "random seed")
      nprint : how often do you wan to print
      obs_anges : rotation angle about z, y, and x-axes (see docs/radiative_trasnfer.pdf)
      obs_distance : distance to the observer
      nxim, nyim : number of pixels of the output image

   2. make
      You need the CFITSIO library installed. (http://heasarc.gsfc.nasa.gov/fitsio/)

   3. MoCafe.x

   4. change directory to idl ("cd idl"), run idl ("idl" or "spear")
       and excute plot_out (".r plot_out")
--------------------------------------------------------------------------------
Note:
   1. If one wants to run the code with a single thread,
      comment out (with #) the line starting with OMP in Makefile.
      The number of threads "NTHREADS" should be specified in Makefile. (NTHREADS=8 for quad-core cpu.)

   2. a) Use main.f90 to run serial code and if you don't care about parallel code.
      b) Use main_omp.f90 and comment out the line starting with OMP in Makefile to run parallel (openmp) code.
      c) Use main_omp.f90, comment out the line starting with OMP in Makefile and set NTHREADS to thread number
          to run serial code and compare with parallel code.
      You should obtain exactly the same results from b) and c) if you assume the same parameters and iseed (/= 0).
      However, you shouldn't expect the same results from a) and b).
      a) and b) would give statistically-consistent results, but not the same.

   3. Use ./docs/isotherm.f90 to generate isothermal density profile.

   4.
         iwp (working precision for integer) is introduced in define.f90 and other routines.
         The random_kiss 32-bit code modified to return 0.0 < u < 1.0
         The original 32-bit code returned 0.0 < u <= 1.0, because of numerial round-off.
