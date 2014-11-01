# dmonly 256: not OLD, DMONLY, 8byte
# full hydro: OLD, not DMONLY, 8byte

# take a choice: 3.09 output format or OLD ones (e.g. sim_aboley)
OLD = 
#OLD = -DOLD

# take a choice: full hydro sim or DM only
DMO = 
#DMO = -DDM

# take a choice: real as 8 byte or 4 byte
NPRE = 8
#NPRE = 4

# take a choice: compiler
FC = gfortran -frecord-marker=4 -O2 -x f95-cpp-input
#FC = f95 

# take a choice: optimized or debugging
DEBUG = -O2
#DEBUG = -O0 -g

# full compiler command:
COPT = -frecord-marker=4 -x f95-cpp-input -Wall
FCC = $(FC) $(COPT) $(DEBUG) -DNPRE=$(NPRE) $(OLD) $(DMO)

#output directory: tools/bin or tools/twiddle/act
OUTDIR = twiddle/act
#OUTDIR = bin/

all: get_temp get_sphere_dm get_sphere_stars amr2map part2map metal2map sod get_particles octreef r2g Makefile

# TODO: check
sod: onsim/sod.f90 Makefile
	$(FCC) onsim/sod.f90 -o $(OUTDIR)/sod

# TODO: check
r2g: onsim/r2g.f90 Makefile
	$(FCC) onsim/r2g.f90 -o $(OUTDIR)/r2g

# TODO: check
octreef: onhalo/octreef.f90 Makefile
	$(FCC) onhalo/octreef.f90 -o $(OUTDIR)/octreef

rho_cell: onsim/rho_cell.f90 Makefile
	$(FCC) onsim/rho_cell.f90 -o $(OUTDIR)/rho_cell

get_sphere_dm: onsim/get_sphere_dm.f90 Makefile
	$(FCC) onsim/get_sphere_dm.f90 -o $(OUTDIR)/get_sphere_dm

get_sphere_stars: onsim/get_sphere_stars.f90 Makefile
	$(FCC) onsim/get_sphere_stars.f90 -o $(OUTDIR)/get_sphere_stars

get_sphere_gas: onsim/get_sphere_gas.f90 Makefile
	$(FCC) onsim/get_sphere_gas.f90 -o $(OUTDIR)/get_sphere_gas

get_temp: onsim/get_temp.f90 Makefile
	$(FCC) onsim/get_temp.f90 -o $(OUTDIR)/get_temp

count_stars: onsim/count_stars.f90 Makefile
	$(FCC) onsim/count_stars.f90 -o $(OUTDIR)/count_stars

get_particles: onsim/get_particles.f90 Makefile
	$(FCC) onsim/get_particles.f90 -o $(OUTDIR)/get_particles

part2map: onsim/part2map.f90 Makefile
	$(FCC) onsim/part2map.f90 -o $(OUTDIR)/part2map

metal2map: onsim/metal2map.f90 Makefile
	$(FCC) onsim/metal2map.f90 -o $(OUTDIR)/metal2map

amr2map: onsim/amr2map.f90
	$(FCC) onsim/amr2map.f90 -o $(OUTDIR)/amr2map

amr2map_dmo: onsim/amr2map.f90
	$(FCC) -DDM onsim/amr2map.f90 -o $(OUTDIR)/amr2map_dmo

test: helpf90/testout.f90 src/testin.f90
	$(FCC) helpf90/testout.f90 -o testout
	$(FCC) helpf90/testin.f90 -o testin
	rm test.dat
	./testout
	./testin

exe:
	chmod ug+x $(OUTDIR)/*.py || chmod ug+x *.py
	chmod ug+x $(OUTDIR)/*.sh
	chmod ug+x *help/*.sh
	chmod ug+x *help/*.py

clean:
	rm -f *~ */*~ a.out */a.out *mod