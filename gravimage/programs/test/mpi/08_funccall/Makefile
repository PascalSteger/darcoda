DEFS=
#-DNPY_NO_DEPRECATED_API -DNPY_1_7_API_VERSION=9 -Wcpp
WARN=
# -Wall #-Werror
CYWARN= #-Werror -Wextra
FOPTS=-fno-underscoring

default: all

grav: grav.py
	cython $(CYWARN) --embed-positions -f --fast-fail -3 grav.py
# rename __pyx_mdef_4grav_3* to pyx_mdef_4grav_3* in grav.c
	sed -i 's/__pyx_mdef/pyx_mdef/g' grav.c
	gcc $(WARN) -I/usr/include/python3.4/ -I/usr/lib64/python3.4/site-packages/numpy/core/include -lpython3.4 -lgfortran -c grav.c


#gcc -Wall -Werror -I/usr/include/python3.4/ -I/usr/lib64/python3.4/site-packages/numpy/core/include -lpython3.4 -lgfortran -c usegrav.c
#gcc -Wall -Werror -I/usr/include/python3.4/ -I/usr/lib64/python3.4/site-packages/numpy/core/include -lpython3.4 -lgfortran usegrav.c grav.c -o usegrav

conn: conn.pyx
	cython $(CYWARN) --embed-positions -f --fast-fail -3 --embed conn.pyx
	sed -i 's/__pyx_mdef/pyx_mdef/g' conn.c
	gcc $(WARN) $(DEFS) -I/usr/include/python3.4/ -I/usr/lib64/python3.4/site-packages/numpy/core/include -lpython3.4 -c conn.c


mesh: fmesh.f90 fmesh_wrapper.f90
	gfortran $(WARN) $(FOPTS) -c fmesh.f90 -o fmesh.o
	gfortran $(WARN) $(FOPTS) -c fmesh_wrapper.f90 -o fmesh_wrapper.o

all: mesh grav conn
	gcc $(WARN) -I/usr/include/python3.4/ -I/usr/lib64/python3.4/site-packages/numpy/core/include conn.o fmesh_wrapper.o fmesh.o grav.o -lpython3.4 -lgfortran -o conn
#	gcc $(WARN) -I/usr/include/python3.4/ -I/usr/lib64/python3.4/site-packages/numpy/core/include conn.o grav.o -lpython3.4 -lgfortran -o conn
