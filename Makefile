#now using v11+ intel ifort

#Intel. Note heap-arrays needs v10+ of compiler, and avoids Seg Faults for large arrays
F90C     = mpif90
#F90C    = ifort

healpix = $(HEALPIX)
LAPACKL = -mkl=sequential -lmpi -lhealpix -openmp

#remove -xHost if cluster is not homogeneous
#add -DHEALPIXI4B if using older healpix and get errors about arguments not matching
FFLAGS = -O3 -xHost -ip -fpp -error-limit 5 -DMPIPIX -DMPI -heap-arrays

#cfitsio = /usr/local/Cluster-Apps/cfitsio/intel/3.300
cfitsio ?= $(CFITSIO)


F90FLAGS = $(FFLAGS) -I$(INCLUDE) -I$(healpix)/include -L$(cfitsio)/lib -L$(healpix)/lib $(LAPACKL) -lcfitsio

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
