#You will need to edit the LAPACKL paths, and edit the options if you are
#not using intel ifc

#Intel ifort 8 and mkl 6.1:
F90C     = mpif77
FFLAGS =  -O2 -Vaxlib -W0 -WB -fpp -tpp7 -xW -ip -DMPIPIX
#use these commented lines instead for non-MPI runs
#F90C     = ifort
#FFLAGS =  -O2 -Vaxlib -W0 -WB -fpp
LAPACKL =  -L/home/antlewis/Healpix_2.00/lib -L/home/antlewis/cfitsio \
       -I/home/antlewis/Healpix_2.00/include -lhealpix -lcfitsio \
	-L/opt/intel/mkl/lib/32 -lmkl_lapack -lmkl_ia32 -lguide -lpthread 


#Digital/Compaq fortran; run with e.g. dmpirun -pf orca.procdist where orca.procdist contains node info
#F90C    = f90
#FFLAGS  = -fpp -lfmpi -lmpi -math_library fast -fpe1 -DMPIPIX  -O4 -tune host -arch host
#LAPACKL = -lcxmlp -L/cita/h/home-3/antlewis/healpix-1.2/lib \
#        -L/cita/h/home-3/antlewis/cfitsio/orca/ \
#        -I/cita/h/home-3/antlewis/healpix-1.2/include -lhealpix -lcfitsio



F90FLAGS = $(FFLAGS) $(INCLUDE) $(LAPACKL)

OBJFILES= inifile.o utils.o spin_alm_tools.o \
   HealpixObj.o HealpixVis.o SimLens.o

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
