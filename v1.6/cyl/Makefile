#******************************** INTEL Fortran **********************
MPI	= 1
#OMP    = 1

ifdef MPI
        MAIN            = main_mpi
        FC              = mpif90
        OMP_FLAGS       =
        MPI_FLAGS       = -DMPI
        exec            = galaxy_mpi.x
else ifdef OMP
        MAIN            = main_omp
        FC              = ifort
        OMP_FLAGS       = -openmp
        MPI_FLAGS       =
        exec            = galaxy.x
else
        MAIN            = main_serial
        FC              = ifort
        OMP_FLAGS       =
        MPI_FLAGS       =
        exec            = galaxy_serial.x
endif
FFLAGS  =	-xhost -fast -cpp $(OMP_FLAGS) $(MPI_FLAGS)
LDFLAGS =	-xhost -fast $(OMP_FLAGS) $(MPI_FLAGS) -lcfitsio -lfftw3f -lfftw3

HOST = $(shell hostname)
ifeq ($(HOST), polaris.kasi.re.kr)
    ifdef MPI
	FC	= mpiifort
    endif
	FFLAGS  = -xhost -ipo -O3 -no-prec-div -fp-model fast=2 -cpp $(OMP_FLAGS) $(MPI_FLAGS)
	LDFLAGS = -xhost -ipo -O3 -no-prec-div -fp-model fast=2 $(OMP_FLAGS) $(MPI_FLAGS) \
		  -lcfitsio -lfftw3f -lfftw3 -L/home/kiseon/local/lib
endif

#FFLAGS  =	-g -xhost -ipo -no-prec-div -O3
#LDFLAGS =	-g -xhost -ipo -no-prec-div -O3 -lcfitsio
#FFLAGS  =	-xhost -fast -DFOLDING
#FFLAGS =       -gen-interface -warn interface -check all -traceback -warn all -fstack-protector
#******************************** GNU Fortran **********************
# Note: FFTW3 in /opt/local/lib was compiled with gcc
#       FFTW3 in /usr/local/lib was compiled with icc
#       They are not compatible.
#OMP    =       -fopenmp -cpp
#FC      =      gfortran
#FFLAGS  =      -O3 $(OMP)
#LDFLAGS =      -O3 -xhost -lcfitsio -lfftw3 -lfftw3f $(OMP) -L/opt/local/lib
#exec	=	galaxy_omp.x
#
#*********************************************************************
.SUFFIXES: .f .f90 .o

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

.f90.o:
	$(FC) $(FFLAGS) -c -o $@ $<

OBJSB	= \
	define.o \
	random.o \
	mathlib.o \
	sersic.o \
	fits_module.o \
	dust_density.o \
	grid_mod.o \
	setup.o \
	output.o \
	raytrace.o \
	generate_photon.o \
	detect.o \
	scatter.o

default:clean
	make main
	make clean

main:$(OBJSB) $(MAIN).o
	$(FC) $(MAIN).o $(OBJSB) $(LDFLAGS) -o $(exec)

clean:
	/bin/rm -f *.o *.mod *.so

distclean: cleanall

cleanall:
	/bin/rm -f *.o *.mod galaxy*.x *_genmod.f90 *.i90
	/bin/rm -f *.fits.gz *.log
