############# Makefile for ramses2gadget
# RF 2010

LIBS=

# --- No MPI, ifort ----------------------------------
#F90 = ifort -g -O0 -assume byterecl 
#FFLAGS = -cpp

# --- Vanilla, gfortran
#F90 = gfortran -O3 -frecord-marker=4 
#FFLAGS = -x f95-cpp-input 

# --- Openmp, gfortran ----------------------------------
F90 = gfortran -O3 -fopenmp -frecord-marker=4
FFLAGS = -x f95-cpp-input -DOPENMP

# --- MPI, pgf90, FULLA ----------------------------------
#source /opt/pgi/setup.sh
#export MPIDIR=/usr/local/mpich2-pgi
#F90 = pgf90 
#FFLAGS = -O3 -Mmpi=mpich2 -Mpreprocess -DMPI

# rules
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@

# linking
ramses2gadget: ramses2gadget.o
	$(F90) ramses2gadget.o -o ramses2gadget $(LIBS)

ramses2gadget_openmp: ramses2gadget_openmp.o
	$(F90) ramses2gadget_openmp.o -o ramses2gadget_openmp $(LIBS)

helloworld: helloworld.o
	$(F90) helloworld.o -o helloworld $(LIBS)

# clearance
clean :
	rm *.o 
