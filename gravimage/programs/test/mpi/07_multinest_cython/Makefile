export F90=mpif90 -ffree-line-length-none -cpp -Wall -Werror -DMPI
export INC=-I/usr/include/python3.4/ -I/usr/lib64/python3.4/site-packages/numpy/core/include
export OBJ=nested_wrapper.o nested.o posterior.o priors.o xmeans_clstr.o kmeans_clstr.o utils1.o randomns.o
export LIB=-L. -L/lib -L/usr -L/usr/lib64 -L/usr/lib64/gcc/x86_64-pc-linux-gnu/4.8.3
export LINK=-llapack -lblas -lpython3.4 -lc -lgfortran

randomns.o: randomns.f90
	$(F90) -c randomns.f90 -o randomns.o

utils1.o: randomns.o utils1.f90
	$(F90) -c utils1.f90 -o utils1.o

priors.o: utils1.o priors.f90
	$(F90) -c priors.f90 -o priors.o

kmeans_clstr.o: randomns.o utils1.o kmeans_clstr.f90
	$(F90) -c kmeans_clstr.f90 -o kmeans_clstr.o

xmeans_clstr.o: utils1.o kmeans_clstr.o xmeans_clstr.f90
	$(F90) -c xmeans_clstr.f90 -o xmeans_clstr.o

posterior.o: utils1.o xmeans_clstr.o posterior.f90
	$(F90) -c posterior.f90 -o posterior.o

nested.o: utils1.o priors.o kmeans_clstr.o xmeans_clstr.o posterior.o nested.f90
	$(F90) -c nested.f90 -o nested.o

nested_wrapper.o: nested_wrapper.f90
	$(F90) -c nested_wrapper.f90 -o nested_wrapper.o

pywrap: nested.o nested_wrapper.o randomns.o utils1.o priors.o kmeans_clstr.o xmeans_clstr.o posterior.o
	cython --embed multinest_demo.pyx
	mpicc $(INC) -c multinest_demo.c -o multinest_demo.o
	mpif90 $(INC) $(LIB) multinest_demo.o $(OBJ) $(LINK) -o multinest_demo

clean:
	rm -rf *.o *.mod a.out *.so

all: clean nested.o
	$(F90) -c nested_wrapper.f90 -o nested_wrapper.o
