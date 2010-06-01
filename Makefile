#now using v11+ intel ifort

#Intel. Note heap-arrays needs v10+ of compiler, and avoids Seg Faults for large arrays
F90C     = ifort

healpix = $(HEALPIX)
LAPACKL = -mkl=sequential -lmkl_lapack -lmpi -lhealpix

FFLAGS = -O3 -ip -fpp -error-limit 5 -DMPIPIX -DMPI -DTHREEJ -heap-arrays

ifndef CFITSIO
cfitsio = /usr/local/cfitsio/intel10/64/3.040/lib
else
cfitsio = $(CFITSIO)
endif

#cosmos seems to have only openmp healpix installed
ifeq ($(COSMOHOST),cosmos)
LINKFLAGS = -openmp
endif

F90FLAGS = $(FFLAGS) -I$(INCLUDE) -I$(healpix)/include -L$(cfitsio) -L$(healpix)/lib $(LAPACKL) -lcfitsio

OBJFILES= toms760.o inifile.o utils.o spin_alm_tools.o \
   HealpixObj.o HealpixVis.o SimLens.o

default: simlens
all: simlens

spin_alm_tools.o:  utils.o toms760.o
HealpixObj.o: spin_alm_tools.o
HealpixVis.o: HealpixObj.o
SimLens.o: HealpixVis.o inifile.o

.f.o:
	f77 $(F90FLAGS) -c $<

%.o: %.f90
	$(F90C) $(F90FLAGS) -c $*.f90

%.o: %.F90
	$(F90C) $(F90FLAGS) -c $*.F90


simlens: $(OBJFILES) 	
	$(F90C) -o simlens $(OBJFILES) $(F90FLAGS) $(LINKFLAGS)

clean:
	rm -f *.o* *.e* *.mod *.d *.pc *.obj core* *.il
