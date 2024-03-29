## MoCafe (Monte Carlo Radiative Transfer)

### Version:
   1.10 - 2023.10.14 \
   1.02 - 2021.11.01

### Purpose:
   Dust Radiative Transfer in cartesian grid system (x, y, z) \
   Basic Monte-Carlo Simulation for a Dust RT in a cartesian grid system. \
   This version is capable of dealing with the polarization of continuum radiation. \
   Unit of the output image = (luminosity unit) cm^-2 sr^-1

### Howto run:

   1. Edit model.in \
      no_photons : number of photons \
      hgg : asymmetry factor for the Henyey-Greenstein function \
      albedo : dust albedo \
      taumax : optical depth from the center to the outer edge \
      outfile : fit file name for the output image \
      xmax, ymax, zmax : dust cloud size (-xmax < x < xmax, etc) \
      nx, ny, nz : number of cells in x, y, z coordinates \
      iseed : initial random number seed (iseed = 0 to generate a random "random seed") \
      nprint : how often do you wan to print \
      alpha, beta, gamma : rotation angle about z, y, and x-axes (see docs/radiative_trasnfer.pdf) \
      distance : distance to the observer \
      nxim, nyim : number of pixels of the output image

   2. make
      You need the CFITSIO library installed. (http://heasarc.gsfc.nasa.gov/fitsio/)

   3. MoCafe.x

### References
If you use MoCafe, please acknowledge the following papers.
   - [Seon & Draine (2016, ApJ, 833, 201)](https://ui.adsabs.harvard.edu/abs/2016ApJ...833..201S/abstract)
   - [Seon (2018, ApJL, 862, 87)](https://ui.adsabs.harvard.edu/abs/2018ApJ...862...87S/abstract)

### Author:
   Kwang-il Seon, KASI
