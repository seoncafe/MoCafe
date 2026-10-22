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

# Serial compiler under the MPI wrapper; the SEDust archive and the HDF5
# install are chosen for it, because Fortran module files are compiler-specific.
ifeq ($(FC), mpif90)
    SERIAL_FC     = gfortran
    HDF5_DEFAULT  = /opt/hdf5_gfortran
else ifeq ($(FC), mpiifx)
    SERIAL_FC     = ifx
    HDF5_DEFAULT  = /opt/hdf5_intel
else
    SERIAL_FC     = ifort
    HDF5_DEFAULT  = /opt/hdf5_intel
endif

# SEDust dust-emission library (Stage 3), self-contained under SEDust/ and made
# by SEDust/sed/build_lib.sh with $(SERIAL_FC).  The archive is compiled against
# the HDF5 Fortran modules with that compiler, so it has to come from the same
# HDF5 and compiler as MoCafe; `make main` therefore rebuilds it whenever HDF5,
# HDF5_PREFIX or SERIAL_FC differs from the values recorded in $(SEDUST_STAMP),
# or a source under SEDust/sed/src/ is newer.  Running build_lib.sh by hand uses
# its own defaults instead.  SEDust's
# 'sed_mathlib' module is named so to avoid clashing with MoCafe's own
# 'mathlib' module.
SEDUST_LIBDIR	?= SEDust/sed/lib
FLAGS		+= -I$(SEDUST_LIBDIR)
SEDUST_LIB	= $(SEDUST_LIBDIR)/libsedust.a
SEDUST_STAMP	= $(SEDUST_LIBDIR)/hdf5_build.txt
SEDUST_BUILD	= HDF5=$(HDF5) HDF5_PREFIX=$(HDF5_PREFIX) FC=$(SERIAL_FC)

# HDF5 installation prefix (set when HDF5=1), by default the install for the
# compiler: /opt/hdf5_intel for ifx/ifort, /opt/hdf5_gfortran for gfortran.
# Override on command line if needed:
#   make HDF5=1 HDF5_PREFIX=/usr/local/hdf5
# Each default is a local rebuild of HDF5 2.2.0 (source and patch in
# /nfs/mocafe/kiseon/MoCafe/new/hdf5-2.2.0/LOCAL_CHANGES) whose static Fortran
# archives contain their own C stubs, so the link line is the 1.14 one.  An
# unmodified 2.x install keeps the C side of the Fortran interface (the h5*_c
# routines) in a separate libhdf5_f90cstub.a, which is linked when present;
# without it such an install leaves 352 of them undefined.  The libaec and
# pthread that libhdf5.settings also lists need no flag here: the shared
# libsz.so pulls in libaec.so, and the MPI wrapper links -lpthread.
HDF5_PREFIX	?= $(HDF5_DEFAULT)

ifneq ($(HDF5), 0)
   FLAGS    += -DHDF5 -I$(HDF5_PREFIX)/include
   HDF5_LIBS = $(HDF5_PREFIX)/lib/libhdf5_fortran.a $(wildcard $(HDF5_PREFIX)/lib/libhdf5_f90cstub.a) \
               $(HDF5_PREFIX)/lib/libhdf5.a -lsz -ldl -lz -lm
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

# libsedust.a reads the SEDust optics products as HDF5, so it must come
# BEFORE the HDF5 libraries: a static archive resolves left to right, and
# it is libsedust that needs those symbols.  The other order leaves ~100
# undefined references.  A text-only archive (HDF5=0 ./build_lib.sh) needs
# none of them and links either way.
# libsedust.a is compiled with OpenMP (build_lib.sh), so the link pulls in the
# OpenMP runtime of the same compiler.
ifeq ($(SERIAL_FC), gfortran)
    OMP_LINK = -fopenmp
else
    OMP_LINK = -qopenmp
endif
LDFLAGS = $(extra) $(FFLAGS) -lcfitsio -L/usr/local/lib $(SEDUST_LIB) $(HDF5_LIBS) $(OMP_LINK)
#*********************************************************************
.SUFFIXES: .f .f90 .o

$(SRCDIR)/%.o: $(SRCDIR)/%.f
	$(FC) $(FFLAGS) $(MODFLAG) $(OMP_FLAGS) -c -o $@ $<

$(SRCDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $(MODFLAG) $(OMP_FLAGS) -c -o $@ $<

OBJS	= \
	define.o \
	qmc_mod.o \
	compose_mod.o \
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
	healpix.o \
	spectrum_sampler_mod.o \
	grain_model_mod.o \
	sed_mod.o \
	random_bulge.o \
	sources_mod.o \
	cellinfo_mod.o \
	allsky_mod.o \
	jtally_mod.o \
	mrw_mod.o \
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
	dustemis_mod.o \
	lucy_mod.o \
	bw_mod.o \
	setup.o \
	output_sum_car.o

OBJSB	= $(addprefix $(SRCDIR)/, $(OBJS))


.PHONY: default main sedust clean distclean cleanall FORCE

default:clean
	$(MAKE) main
	$(MAKE) clean

main:$(SEDUST_LIB) $(OBJSB) $(SRCDIR)/$(MAIN).o
	$(FC) $(SRCDIR)/$(MAIN).o $(OBJSB) $(LDFLAGS) $(OMP_FLAGS) -o MoCafe.x

sedust: $(SEDUST_STAMP)
	cd SEDust/sed && $(SEDUST_BUILD) ./build_lib.sh

# The stamp is rewritten only when the HDF5 or compiler selection changes, so its time
# stamp moves, and the archive is rebuilt, only then.  It is written before
# the build, so the archive it describes is always the newer file.
$(SEDUST_STAMP): FORCE
	@mkdir -p $(SEDUST_LIBDIR)
	@echo '$(SEDUST_BUILD)' | cmp -s - $@ || echo '$(SEDUST_BUILD)' > $@

$(SEDUST_LIB): $(SEDUST_STAMP) $(wildcard SEDust/sed/src/*.f90) SEDust/sed/build_lib.sh
	$(MAKE) sedust

FORCE:

clean:
	/bin/rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod *.o *.mod

distclean: cleanall

cleanall:
	/bin/rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod *.o *.mod MoCafe.x
	/bin/rm -f *.fits.gz *.log
