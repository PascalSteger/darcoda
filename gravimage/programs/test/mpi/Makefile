FC = gfortran
FFLAGS = -cpp -L. -L/usr/lib64/gcc/x86_64-pc-linux-gnu/4.8.3 -L/usr -L/usr/lib -I. -ffree-line-length-none -O2
LIBS = -L. -L/lib -L/usr -L/usr/lib -L/usr/lib64/gcc/x86_64-pc-linux-gnu/4.8.3
INCS = -I. -I/usr/include

LINKLIBS = -lm -lc -llapack -lblas -lgfortran

NSOBJECTS = randomns.o utils1.o priors.o kmeans_clstr.o xmeans_clstr.o posterior.o nested.o

%.o: %.f90
	$(FC) -fPIC $(INCS) $(LIBS) $(FFLAGS) $(LINKLIBS) -c -o $@ $^

default: libnest3.a libnest3.so Makefile
	cp libnest3.* /usr/local/lib64/

all: libnest3.a libnest3.so eggboxC

libnest3.so: $(NSOBJECTS)
	ld -shared $(LIBS) -o $@ $^  $(LINKLIBS)

libnest3.a: $(NSOBJECTS)
	ar r $@ $^

eggboxC:
	make -C example_eggbox_C

clean:
	-rm libnest3.*  *.o *.mod

cleanall: clean_example_eggbox_C

clean_exec:
	-rm eggboxC

clean_example_eggbox_C:
	make -C example_eggbox_C clean
