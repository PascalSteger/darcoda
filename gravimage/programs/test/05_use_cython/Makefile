fmesh.o: fmesh.f90
	gfortran -fPIC -c fmesh.f90 -o fmesh.o

fmesh_wrapper.o: fmesh_wrapper.f90
	gfortran -c fmesh_wrapper.f90 -o fmesh_wrapper.o

cwrap: fmesh.o fmesh_wrapper.o
	gcc -c wrap_mesh.c
	gcc wrap_mesh.o fmesh_wrapper.o fmesh.o -o wrap_mesh

pywrap: fmesh.o fmesh_wrapper.o
	cython --embed mesh_use.pyx
	gcc -I/usr/include/python3.4/ -I/usr/lib64/python3.4/site-packages/numpy/core/include mesh_use.c fmesh_wrapper.o fmesh.o -lpython3.4 -o mesh_use

clean:
	rm -rf *.o *.mod *.c a.out *.so

all: comp
	python3 setup.py build_ext --inplace
	python3 launch.py
