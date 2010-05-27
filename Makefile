#You will need to edit the LAPACKL paths, and edit the options if you are
#not using v11+ intel ifort

#Intel. Note heap-arrays needs v10+ of compiler, and avoids Seg Faults for large arrays
F90C     = mpif90

FFLAGS = -ip -O3 -fpp -error-limit 5 -DMPIPIX -DMPI -heap-arrays -mkl=parallel
#use these lines instead for non-MPI runs
#F90C     = ifort 
#FFLAGS =  -ip -O3 -fpp -error-limit 5 -heap-arrays -mkl=parallel
LAPACKL = -L/usr/local/cfitsio/intel10/64/3.040/lib \
	 -L/usr/local/healpix/intel10/64/2.01/serial/lib \
	 -I/usr/local/healpix/intel10/64/2.01/serial/include \
	 -lhealpix -lcfitsio -lmkl_lapack -lguide -lpthread 



F90FLAGS = $(FFLAGS) -I$(INCLUDE) $(LAPACKL)

OBJFILES= toms760.o inifile.o utils.o spin_alm_tools.o \
   HealpixObj.o HealpixVis.o SimLens.o

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
	$(F90C) -o simlens $(OBJFILES) $(F90FLAGS)

clean:
	rm -f *.o* *.e* *.mod *.d *.pc *.obj core* *.il
