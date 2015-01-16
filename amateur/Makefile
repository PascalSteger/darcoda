CC	= g++
# for optimized execution of program
COFLAGS = -O3 -Wall -Wfatal-errors
# for debug purposes
CDFLAGS = -O0 -Wall -Wfatal-errors -g -fpermissive

#-msse2
#-DSINGLE_PRECISION -DGAUSSIAN_KERNEL

# run on 64bit machine achtland
AIPATH	= -I. -I/data/achtland1/psteger/local/include/
ALPATH	= -L/usr/local/lib/ -L/data/achtland1/psteger/local/lib/
ALFLAGS	= -lhdf5 -lsrfftw -lsfftw -lgsl -lgslcblas -lm -lANN
# -lcblas

# run on 64bit machine plumpy
PIPATH  = -I. -I/usr/include/ -I/data/achtland1/psteger/localp/include/
PLPATH  = -L/usr/local/lib/ -L/data/achtland1/psteger/localp/lib/
PLFLAGS = -lhdf5 -lfftw3 -lgsl -lgslcblas -lm -lANN

# run on 64bit machine igloo (AMD Phenom 9750, 8GB RAM)
SIPATH  = -I. -I/usr/include/ -I/usr/local/include/
SLPATH  = -L/usr/lib -L/usr/local/lib/
SLFLAGS = -lhdf5 -lfftw3 -lgsl -lgslcblas -lm -lANN

# run on 32bit machine morrigan
MIPATH  = -I. -I/usr/include/ -I/usr/local/include/
MLPATH	= -L/usr/local/lib/ -L/home/hahn/local/lib
MLFLAGS = -lhdf5 -lfftw3 -lgsl -lgslcblas -lm -lANN

# run on 32bit machine snowball
NIPATH = -I. -I/home/ps/amd/local/include
NLPATH = -L/usr/local/lib -L/home/ps/amd/local/lib
NLFLAGS	= -lhdf5 -lfftw3 -lgsl -lgslcblas -lm -lANN

# choose previously introduced compilation setups
IPATH  = $(PIPATH)
LPATH  = $(PLPATH)
LFLAGS = $(PLFLAGS)
# choose: COFLAGS for fast execution or CDFLAGS for debug purposes
MYC     = $(CC) $(CDFLAGS) $(IPATH)
#MYC     = $(CC) $(COFLAGS) $(IPATH)

FILES = Global.cpp ReadGadget.c ExtractSnap.cpp TidalField.cpp	   \
AnalyzeSnap.cpp AHFrun.cpp MergerTreeRun.cpp SubhaloMatch.cpp	   \
HaloTrackerRun.cpp ScalarField.cpp GreensFunction.cpp Vector3.cpp  \
Matrix33.cpp Array3.cpp Eigensystem3.cpp Halo.cpp Distance.cpp	   \
Cosmology.cpp Simulation.cpp HResults.cpp AResults.cpp AOutput.cpp \
Output.cpp amateur.cpp
HEADERS = Global.h ReadGadget.h ExtractSnap.h TidalField.h	       \
AnalyzeSnap.h AHFrun.h MergerTreeRun.h SubhaloMatch.h HaloTrackerRun.h \
ScalarField.h GreensFunction.h Vector3.h Matrix33.h Eigensystem3.h     \
Array3.h Halo.h Distance.h Cosmology.h Simulation.h HDFIO.h NFWFit.h   \
MySTL.h OctTree.h HResults.h AResults.h AOutput.h Output.h

OBJS = Global.o ReadGadget.o ExtractSnap.o TidalField.o AnalyzeSnap.o  \
AHFrun.o MergerTreeRun.o SubhaloMatch.o HaloTrackerRun.o ScalarField.o \
GreensFunction.o Vector3.o Matrix33.o Eigensystem3.o Array3.o Halo.o   \
Distance.o Cosmology.o Simulation.o HResults.o AResults.o AOutput.o    \
Output.o amateur.o

all: amateur Makefile #needed to recompile if Makefile/Doxyfile is changed

doc: Doxyfile tags Makefile $(FILES) $(HEADERS)
	doxygen Doxyfile

tags:
	etags *.cpp *.h

amateur: $(OBJS) Makefile
	$(MYC) $(OBJS) -o amateur $(LPATH) $(LFLAGS)
	echo 'Success!'

amateur.o: amateur.cpp AHFrun.h MergerTreeRun.h SubhaloMatch.h HaloTrackerRun.h AnalyzeSnap.h ExtractSnap.h TidalField.h Global.h Makefile
	$(MYC) -c amateur.cpp

Global.o: Global.cpp Global.h Makefile
	$(MYC) -c Global.cpp


ExtractSnap.o: ExtractSnap.cpp ReadGadget.c ReadGadget.h HDFIO.h Global.h Makefile
	$(MYC) -c ExtractSnap.cpp

AHFrun.o: AHFrun.cpp Global.h Makefile
	$(MYC) -c AHFrun.cpp

MergerTreeRun.o: MergerTreeRun.cpp Global.h Makefile
	$(MYC) -c MergerTreeRun.cpp

SubhaloMatch.o: SubhaloMatch.cpp Global.h Makefile
	$(MYC) -c SubhaloMatch.cpp

HaloTrackerRun.o: HaloTrackerRun.cpp Global.h Makefile
	$(MYC) -c HaloTrackerRun.cpp

AnalyzeSnap.o: AnalyzeSnap.cpp Vector3.h Matrix33.h Eigensystem3.h Halo.o ReadGadget.h HDFIO.h Halo.h Simulation.h HResults.h AResults.h Output.h Global.h Makefile
	$(MYC) -c AnalyzeSnap.cpp

TidalField.o: TidalField.cpp TidalField.h Vector3.h Matrix33.h Eigensystem3.h Array3.h GreensFunction.h ScalarField.h PoissonSolver.h ReadGadget.h HDFIO.h Global.h Makefile
	$(MYC) -c TidalField.cpp

ReadGadget.o: ReadGadget.c ReadGadget.h Makefile
	$(MYC) -c ReadGadget.c

ScalarField.o: ScalarField.cpp ScalarField.h Makefile
	$(MYC) -c ScalarField.cpp

Vector3.o: Vector3.cpp Vector3.h Global.h Makefile
	$(MYC) -c Vector3.cpp

Matrix33.o: Matrix33.cpp Matrix33.h Vector3.h Global.h Makefile
	$(MYC) -c Matrix33.cpp

Array3.o: Array3.cpp Array3.h Makefile
	$(MYC) -c Array3.cpp

Eigensystem3.o: Eigensystem3.cpp Eigensystem3.h Vector3.h Global.h Makefile
	$(MYC) -c Eigensystem3.cpp

Halo.o: Halo.cpp Halo.h Vector3.h Matrix33.h Simulation.h Distance.h Global.h Makefile
	$(MYC) -c Halo.cpp

Distance.o: Distance.cpp Distance.h Global.h Makefile
	$(MYC) -c Distance.cpp

Simulation.o: Simulation.cpp Simulation.h Cosmology.h Global.h Makefile
	$(MYC) -c Simulation.cpp

Cosmology.o: Cosmology.cpp Cosmology.h Global.h Makefile
	$(MYC) -c Cosmology.cpp

HResults.o: HResults.cpp HResults.h Vector3.h Eigensystem3.h Global.h Makefile
	$(MYC) -c HResults.cpp

AResults.o: AResults.cpp AResults.h Vector3.h Eigensystem3.h Global.h Makefile
	$(MYC) -c AResults.cpp

AOutput.o: AOutput.cpp AOutput.h Simulation.h Output.h Makefile
	$(MYC) -c AOutput.cpp

Output.o: Output.cpp Output.h AResults.h HResults.h Simulation.h Global.h Makefile
	$(MYC) -c Output.cpp

SNAP = 3

run: amateur
	./amateur $(SNAP)

# x, a, m, s, h, z, t
runx: amateur
	./amateur $(SNAP) -x
runa: amateur
	./amateur $(SNAP) -a
runo: amateur
	./amateur $(SNAP) -o
runm: amateur
	./amateur $(SNAP) -m
runs: amateur
	./amateur $(SNAP) -s
runh: amateur
	./amateur $(SNAP) -h
runz: amateur
	./amateur $(SNAP) -z
runt: amateur
	./amateur $(SNAP) -t

clean: 
	rm -rf *.o

tidy: clean
	rm -rf amateur
	rm -rf *~
	rm -rf *.gch
	rm -rf *.bak

batch: amateur
	./amateur 15 -z
	./amateur 16 -z
	./amateur 18 -z