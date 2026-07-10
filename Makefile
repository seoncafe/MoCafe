#************* Choose Compiler *************
HOST  = $(shell hostname -s)
UNAME = $(shell uname)

ifeq ($(UNAME), Linux)
    FC      = mpiifort
    ifeq (, $(shell which ifort))
        FC  = mpiifx
    endif
    ifeq (, $(shell which mpiifx))
        FC  = mpif90
    endif
endif
ifeq ($(F90), mpif90)
    FC      = mpif90
endif
ifeq ($(F90), mpiifx)
    FC      = mpiifx
endif
ifeq ($(F90), mpiifort)
    FC      = mpiifort
endif

#---------------------------
# Fortran sources, objects (*.o) and modules (*.mod) all live in src/;
# only the final MoCafe.x is written to the top-level directory.
SRCDIR		= src

MAIN		= main
FLAGS		= -cpp -DMPI
DEBUG		= 0
HDF5		= 1

# HDF5 installation prefix (set when HDF5=1). Override on command line if needed:
#   make HDF5=1 HDF5_PREFIX=/usr/local/hdf5
HDF5_PREFIX	?= /data/opt/hdf5_intel

ifneq ($(HDF5), 0)
   FLAGS    += -DHDF5 -I$(HDF5_PREFIX)/include
   HDF5_LIBS = $(HDF5_PREFIX)/lib/libhdf5_fortran.a $(HDF5_PREFIX)/lib/libhdf5.a -lsz -ldl -lz -lm
endif

ifeq ($(FC), $(filter $(FC), mpiifort mpiifx))
    FFLAGS  = -ipo -O3 -no-prec-div -fp-model fast=2 $(FLAGS) $(Fextra)
    MODFLAG = -module $(SRCDIR)
    LDFLAGS = $(extra) $(FFLAGS) -lcfitsio
    ifeq ($(HOST), polaris)
       LDFLAGS = $(FFLAGS) -lcfitsio -L$(HOME)/local/lib
    endif
else ifeq ($(FC), mpif90)
    #--- GNU compilers
    FFLAGS = -O3 -ffpe-summary=none $(FLAGS) $(EXTRAFLAG)
    MODFLAG = -J$(SRCDIR) -I$(SRCDIR)
else ifeq ($(UNAME), Darwin)
    FC      = mpif90
    FFLAGS  = -O3 $(ext_FLAGS) $(FLAGS)
    MODFLAG = -J$(SRCDIR) -I$(SRCDIR)
    LDFLAGS = $(extra) $(FFLAGS) -lcfitsio
endif

ifeq ($(DEBUG), 1)
   ifeq ($(FC), $(filter $(FC), mpiifort mpiifx))
      FFLAGS  = -check all,noarg_temp_created -fpe0 -debug all -traceback -g -O0 $(FLAGS)
      #FFLAGS  = -check all,arg_temp_created -fpe0 -debug all -traceback -g -O0 $(FLAGS)
      #FFLAGS  = -check all,noarg_temp_created -fpe0 -debug all -traceback -g -assume ieee_fpe_flags -O0 -xhost $(FLAGS)
   else
      FFLAGS  = -O0 -fimplicit-none  -Wall -Wextra -fcheck=all -std=f2008  -pedantic -fall-intrinsics -fbacktrace $(FLAGS)
   endif
endif

LDFLAGS = $(extra) $(FFLAGS) -lcfitsio -L/usr/local/lib $(HDF5_LIBS)
#*********************************************************************
.SUFFIXES: .f .f90 .o

$(SRCDIR)/%.o: $(SRCDIR)/%.f
	$(FC) $(FFLAGS) $(MODFLAG) $(OMP_FLAGS) -c -o $@ $<

$(SRCDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $(MODFLAG) $(OMP_FLAGS) -c -o $@ $<

OBJS	= \
	define.o \
	utility.o \
	scan_mod.o \
	random_mt.o \
	memory_mod_mpi.o \
	fitsio_mod.o \
	hdf5io_mod.o \
	iofile_mod.o \
        read_mod.o \
	clump_mod.o \
	physics_amr_mod.o \
	octree_mod.o \
	read_generic_amr.o \
	mathlib.o \
	sed_mod.o \
	jtally_mod.o \
	raytrace_car.o \
	raytrace_clump.o \
	raytrace_amr.o \
	observer_mod.o \
	grid_mod_car.o \
	grid_mod_clump.o \
	grid_mod_amr.o \
	write_mod_car.o \
	peelingoff_mod.o \
	external_radiation.o \
	gen_photon_car.o \
	run_simulation_mod.o \
	sightline_tau_mod.o \
	scattering_car.o \
	setup.o \
	output_sum_car.o

OBJSB	= $(addprefix $(SRCDIR)/, $(OBJS))

default:clean
	$(MAKE) main
	$(MAKE) clean

main:$(OBJSB) $(SRCDIR)/$(MAIN).o
	$(FC) $(SRCDIR)/$(MAIN).o $(OBJSB) $(LDFLAGS) $(OMP_FLAGS) -o MoCafe.x

clean:
	/bin/rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod *.o *.mod

distclean: cleanall

cleanall:
	/bin/rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod *.o *.mod nebula.x
	/bin/rm -f *.fits.gz *.log
