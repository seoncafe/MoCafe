#******************************** INTEL Fortran **********************
MAIN		= main
ext_FLAGS	= -cpp -DMPI

DEBUG		= 0

HOST  = $(shell hostname -s)
UNAME = $(shell uname)
ifeq ($(UNAME), Linux)
    CC      = mpiicc
    FC      = mpiifort
    #FC      = mpif90
    CFLAGS  = -ipo -O3 -no-prec-div -fp-model fast=2
    FFLAGS  = -ipo -O3 -no-prec-div -fp-model fast=2 $(ext_FLAGS)
    #FFLAGS  = -xhost -ipo -O3 -no-prec-div -fp-model fast=2 $(ext_FLAGS)
    #FFLAGS  = -O3 $(ext_FLAGS)
    LDFLAGS = $(FFLAGS) -lcfitsio
    ifeq ($(HOST), polaris)
       LDFLAGS = $(FFLAGS) -lcfitsio -L$(HOME)/local/lib
    endif
else ifeq ($(UNAME), Darwin)
    # Disable the warnings for a while (Big Sur gives the warnings. 2020-11-15).
    # this is not required anymore (2021.01.12)
    #ifeq ($(shell uname -r), 20.1.0)
    #ifeq ($(shell uname -r), 20.2.0)
    #   extra = -diag-disable 11012,11013,11021
    #endif
    CC      = mpicc
    FC      = mpif90
    # 10440 warning: Use of IPO with macOS* SDK 11 or higher is not supported, disabling the warning. (2021.01.12)
    CFLAGS  = -xhost -fast -diag-disable 10440,10441
    FFLAGS  = -xhost -fast -diag-disable 10440,10441 $(FLAGS)
    #FFLAGS  = -g -O0 $(ext_FLAGS)
    LDFLAGS = $(extra) $(FFLAGS) -lcfitsio
endif

ifeq ($(DEBUG), 1)
   FFLAGS  = -check all,noarg_temp_created -fpe0 -debug all -traceback -g -O0 $(ext_FLAGS)
   #FFLAGS  = -check all,noarg_temp_created -fpe0 -debug all -traceback -g -assume ieee_fpe_flags -O0 -xhost $(ext_FLAGS)
   LDFLAGS = $(FFLAGS) -lcfitsio
   ifeq ($(HOST), polaris)
      LDFLAGS = $(FFLAGS) -lcfitsio -L$(HOME)/local/lib
   endif
endif

#*********************************************************************
.SUFFIXES: .c .f .f90 .o

.c.o:
	$(CC) $(CFLAGS) $(OMP_FLAGS) -c -o $@ $<

.f.o:
	$(FC) $(FFLAGS) $(OMP_FLAGS) -c -o $@ $<

.f90.o:
	$(FC) $(FFLAGS) $(OMP_FLAGS) -c -o $@ $<

OBJSB	= \
	define.o \
	utility_c.o \
	utility.o \
	random_mt_v2.52.o \
	memory_mod_mpi_v2.o \
	fitsio_mod.o \
        read_mod.o \
	mathlib.o \
	raytrace_car.o \
	observer_mod.o \
	grid_mod_car.o \
	write_mod_car.o \
	peelingoff_mod.o \
	gen_photon_car.o \
	run_simulation_mod_v2.o \
	scattering_car.o \
	setup.o \
	output_sum_car.o \

default:clean
	make main
	make clean

main:$(OBJSB) $(MAIN).o
	$(FC) $(MAIN).o $(OBJSB) $(LDFLAGS) $(OMP_FLAGS) -o MoCafe.x

clean:
	/bin/rm -f *.o *.mod

distclean: cleanall

cleanall:
	/bin/rm -f *.o *.mod nebula.x
	/bin/rm -f *.fits.gz *.log
