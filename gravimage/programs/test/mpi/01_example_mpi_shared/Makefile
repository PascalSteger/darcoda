#mpif90 -fPIC fd.f90
#mpif90 -fPIC itimes.f90
all:
	mpif90 -shared -fPIC fd.f90 itimes.f90 -o it.so

clean:
	rm -rf *.o *.so *.out *.mod
