############# Makefile for ramses2gadget
# PS 2012
#############################################################################

#F90 = ifort -O2 -assume byterecl -DMPI
F90 = mpif90.openmpi -DMPI -DNPRE=4
FFLAGS = -cpp

#F90 = gfortran -O3 -fopenmp -frecord-marker=4
#FFLAGS = -x f95-cpp-input -DOPENMP


all: r2g checkr2g Makefile

r2g: r2g.f90 Makefile
	$(F90) $(FFLAGS) r2g.f90 -o r2g

checkr2g: checkr2g.f90 Makefile
	$(F90) $(FFLAGS) checkr2g.f90 -o checkr2g

# clearance
clean :
	rm *.o
	rm *~
tidy: clean
	rm r2g checkr2g
