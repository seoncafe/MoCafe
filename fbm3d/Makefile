#******************************** INTEL Fortran **********************
#OMP	=	-fopenmp
FC      =	ifort
UNAME = $(shell uname)
ifeq ($(UNAME), Linux)
	FFLAGS	= -xhost -ipo -O3 -no-prec-div -fp-model fast=2 -fpp $(OMP)
else
	FFLAGS	= -xhost -fast -fpp $(OMP)
endif
ifdef OMP
	LDFLAGS =	${FFLAGS} -lfftw3f_omp -lfftw3f -lcfitsio
else
	LDFLAGS =	${FFLAGS} -lfftw3f -lcfitsio
endif

#********************************
.SUFFIXES: .f .f90 .o

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

.f90.o:
	$(FC) $(FFLAGS) -c -o $@ $<

OBJS	= \
	define.o \
	random_light.o \
	fbm3d.o \
	write_output.o

default: make_fbm3d.x
	/bin/rm -rf *.o *.mod

make_fbm3d.x:$(OBJS) make_fbm3d.o
	$(FC) make_fbm3d.o $(OBJS) $(LDFLAGS) -o make_fbm3d.x
	/bin/rm -rf *.o *.mod

clean:
	/bin/rm -f *.o *.mod

cleanall:
	/bin/rm -f *.o *.mod make_fbm3d.x
