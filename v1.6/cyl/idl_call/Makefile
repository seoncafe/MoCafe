#-------------------------------- INTEL Fortran --------------------------
SRC_DIR	=	..
FC      =	ifort
CC      =	icc
LD	=	ifort

FFLAGS  =	-xHost -ipo -O3 -no-prec-div -fp-model fast=2 -openmp -cpp -I$(SRC_DIR)/
LDFLAGS =	-xHost -ipo -O3 -no-prec-div -fp-model fast=2 -openmp -lcfitsio -lfftw3 -lfftw3f
CFLAGS  =	-xHost -ipo -O3 -no-prec-div -fp-model fast=2 -I/usr/local/exelis/idl82/external/include

HOST = $(shell hostname)
ifeq ($(HOST), fims.kasi.re.kr)
   SOFLAGS =	-fPIC -shared
else
   SOFLAGS =	-fPIC -shared -bundle -flat_namespace -undefined suppress
endif

# for debug
#FFLAGS  =	-xhost -fast -fpp -DDEBUG
#LDFLAGS =	-xhost -fast -fpp -DDEBUG
#FFLAGS  =	-g -fpp
#LDFLAGS =	-g -fpp
#FFLAGS  =	-g -fopenmp -fpp -warn all -traceback
#LDFLAGS =	-g -fopenmp -fpp -warn all -traceback
#*********************************************************************
.SUFFIXES: .c .f90 .o

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

.f90.o:
	$(FC) $(FFLAGS) -c -o $@ $<

OBJS2	= \
	$(SRC_DIR)/define.o \
	$(SRC_DIR)/random.o \
	$(SRC_DIR)/mathlib.o \
	$(SRC_DIR)/sersic.o \
	$(SRC_DIR)/fits_module.o \
	$(SRC_DIR)/generate_photon.o \
	$(SRC_DIR)/dust_density.o \
	$(SRC_DIR)/grid_mod.o \
	$(SRC_DIR)/setup.o \
	./setup_idl.o \
	$(SRC_DIR)/output.o \
	$(SRC_DIR)/raytrace.o \
	$(SRC_DIR)/detect.o \
	$(SRC_DIR)/scatter.o

galaxy.so: clean $(OBJS2) galaxy_idl.o galaxy_idl_load.o
	@/bin/rm -f *.so
	$(LD) galaxy_idl_load.o galaxy_idl.o $(OBJS2) $(LDFLAGS) -o galaxy_idl.so $(SOFLAGS)
	#cp galaxy_idl.dlm $(expldust)/Galaxy
	#cp galaxy_idl.so $(expldust)/Galaxy
	@make clean

clean:
	@/bin/rm -f $(SRC_DIR)/*.mod $(SRC_DIR)/*.o *.mod *.o *.dcl

cleanall:
	@/bin/rm -f $(SRC_DIR)/*.mod $(SRC_DIR)/*.o *.mod *.o
	@/bin/rm -f galaxy
	@/bin/rm -f *.fits.gz out*.log *_genmod.f90
	@/bin/rm -f *.so
