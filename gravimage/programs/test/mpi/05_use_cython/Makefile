fmesh.o: fmesh.f90
	gfortran -fPIC -c fmesh.f90 -o fmesh.o

fmesh_wrapper.o: fmesh_wrapper.f90
	gfortran -c fmesh_wrapper.f90 -o fmesh_wrapper.o

cwrap: fmesh.o fmesh_wrapper.o
	gcc -c wrap_mesh.c
	gcc wrap_mesh.o fmesh_wrapper.o fmesh.o -lgfortran  -o wrap_mesh

pywrap: fmesh.o fmesh_wrapper.o
	cython --embed mesh_use.pyx
	gcc -I/usr/include/python3.4/ -I/usr/lib64/python3.4/site-packages/numpy/core/include mesh_use.c fmesh_wrapper.o fmesh.o -lpython3.4 -lgfortran -o mesh_use

test: testf.f90 testc.c Makefile
	gfortran -c testf.f90
	gcc -c testc.c
	gfortran -o test testf.o testc.o
	./test

clean:
	rm -rf *.o *.mod a.out *.so


all: pywrap
