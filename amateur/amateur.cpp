/****************************************************************************
 *   amateur v.1.4                                                          *
 *   analyze snapshots from Gadget 2 simulations                            *
 *   Copyright (C) 2008-2010 by Pascal Steger, and Dr. Oliver Hahn          *
 *   psteger@phys.ethz.ch                                                   *
 ****************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <string>

#include <sys/resource.h>

#include "Global.h"
#include "ExtractSnap.h"
#include "AnalyzeSnap.h"
#include "AHFrun.h"
#include "AOutput.h"
#include "MergerTreeRun.h"
#include "HaloTrackerRun.h"
#include "SubhaloMatch.h"
#include "TidalField.h"
// for testing purposes:
#include "PoissonSolver.h" 
#include "ScalarField.h"
#include "Array3.h"

/**
 * \file amateur.cpp
 * \brief Main program, calls all subroutines.
 */

/**
 * \brief Help for use of commandline arguments.
 * example program calls:
 * amateur 10 -t looks at tarkin_13309_CSF_45x_snap_211 and computes the tidal field
 * amateur 1	runs all parts on tarkin_5726_10x_snap_302
 * 
 * program runs on nice value 15 by default
 * attention: option    -a requires a lot of CPU power, use a fast comoputer via ssh!
 *                      -z requires a lot of memory
 *                      -s and -x show excessive file access.
 *                         For fast execution, use local computer (no ssh).
 */
void HelpUser( void ) {
  writeln( "Usage: amateur SNAPSHOTID{1..?} [xamshzt]");
  writeln( "where SNAPSHOT is an integer between 1 and 14");
  writeln( "-a runs AMIGA's subhalo finder AHFstep");
  writeln( "-m runs AMIGA's subhalo <=> halo identifier MergerTree");
  writeln( "-o extracts relevant results from AHF output");
  writeln( "-s matches subhalos to their host halo");
  writeln( "-h runs AMIGA's HaloTracker (experimental)");
  writeln( "-t computes the tidal field of the given snapshot at 10 smoothing scales");
  writeln( "-x extracts information to a HDF5 file");
  writeln( "-z analyses the snapshot file, using AMIGA's subhalos" );
  writeln( "-T test runs for debugging purposes" );
  writeln( "default: runs successively x, a, o, m, s, (h), z, t" );
  for(int i=0; i<18; ++i)
    std::cout << "  " << i+1 << " " << Snaps(i) << std::endl;
  std::cout << std::endl;
  //printMinMaxValues();
}

/**
 * \brief main routine, calls different tools
 * The main routine of amateur presents the user a tool to select 
 * from the different possibilities to run the program, 
 * and on which snapshot file a command should be performed.
 * If no arguments are given, amateur prints a help message 
 * showing the command syntax with all possible snapshots
 * @param argc: int; number of arguments
 * @param argv: char* []; character array of arguments
 */
int main ( int argc, char *argv[ ] ) {

  //int setpriority(int which, id_t who (0 = this process), int nice);"
  setpriority(PRIO_PROCESS, 0, 0);//15);
  writeln("");
  writeln( "AMATEUR by Pascal Steger,    psteger@phys.ethz.ch");
  writeln( "=================================================");
  writeln( "");
  if ( argc <= 1 || argc > 3 ) {
    HelpUser();
    return 0;
  }
  //Array3 A(512,512,512);
  //ScalarField S(512,512,512,0.3,3);
  //PoissonSolver<3> ps(&S,&S);

  int isnap = atoi( argv[ 1 ] ) - 1;
  std::cout << "	Looking at snapshot " << Snaps( isnap ) << "." << std::endl;

  for ( int i = 0; i < 1; ++i ) {
    // std::cout << "commandline argument nr. " << i << std::endl;
    char c = getopt( argc, argv, "aomshtTxz" );
    switch ( c ) {
    case 'a':
      writeln("	Running AMIGA's AHFstep on snapshot file...");
      AHFrun( isnap );
      break;
    
    case 'o':
      writeln(" Extracting interesting output columns from AHF output...");
      AOutput( isnap );
      break;

    case 'm':
      writeln("	Searching halo/subhalo-relations with OwnMergerTree...");
      MergerTreeRun( isnap );
      break;

    case 's':
      writeln("	Matching subhalos to their host..." );
      SubhaloMatch( isnap );
      break;

    case 'h':
      writeln(" Running HaloTracker...");
      //HaloTrackerRun( isnap );
      break;

    case 't':
      writeln( " Computing tidal field..." );
      //std::cout << "Enter N for N^3 sample points: ";
      //int N;
      //std::cin >> N;
      //TidalField(isnap, N, 500.0);
      TidalField( isnap, 512, 500.0 );
      break;

    case 'x':
      writeln( " Extracting GADGET2 binary to HDF5..." );
      ExtractSnap( isnap );
      break;

    case 'z':
      writeln( " Analyzing snapshot file with AHF Halos..." );
      AnalyzeSnap( isnap );
      break;

    case 'T':
      writeln( " Running test routines..." );
//      writeln( " starting with setting values in A..." );
//      for(int i=0; i<512; ++i){
//	A.setValue(i,i,i,0.1*i);
//	std::cout << A.getValue(i,i,i);
//      }
//      cout << endl;
//      for(int i=0; i<512; ++i){
//	A.setValue(i,i,i,0.1*i);
//	std::cout << A.getValue(i,i,i);
//      }
//      cout << endl;

//      writeln( " now for ScalarField...");
//      for(int i=0; i<512; ++i){
//	S.setValue(i,i,i,i/1.0);
//	cout << S.getValue(i,i,i) << " ";
//      }
//      cout << endl;

//      writeln( " test routine from ScalarField..." );
//      S.test();

//      writeln( " test routine from PoissonSolver..." );
//      ps.test();
//      for(int i=0; i<512; ++i){
//	cout << S.getValue(i,i,i) << " ";
//      }
	cout << endl;
      break;

    default:
      writeln( " Default: Extracting GADGET2 binary to HDF5..." );
      ExtractSnap( isnap );
      writeln( " Running AMIGA's AHFstep on snapshot file..." );
      AHFrun( isnap );
      writeln( " Extracting interesting output columns from AHF output...");
      AOutput( isnap );
      writeln( " Searching halo/subhalo-relations with OwnMergerTree..." );
      MergerTreeRun( isnap );
      writeln( " Matching subhalos to their host..." );
      SubhaloMatch( isnap );
      writeln( " Running HaloTracker...");
      // HaloTrackerRun(isnap);
      writeln( " Analyzing snapshot file with Halo..." );
      AnalyzeSnap( isnap );
      writeln( " Computing tidal field..." );
      TidalField( isnap, 512, 500.0 );
    }
  }
  return EXIT_SUCCESS;
}

// The rest of this file contains documentation for compiling and
// running the code, in a format appropriate for doxygen.

/**
\mainpage Reference documentation for AMATEUR

\author Pascal S. P. Steger \n
        Department of Physics \n
        Institute of Astronomy \n
        ETH Zurich \n
        Switzerland \n 
        psteger@phys.ethz.ch \n

\author Oliver J. Hahn \n
        Department of Physics \n
        Institute of Astronomy \n
        ETH Zurich \n
        Switzerland \n 
        ohahn@phys.ethz.ch \n
\n

\section prelim Getting started

\b AMATEUR is a collection of tools to analyze snapshots from
cosmological simulations performed with GADGET-2. It uses AHFstep from
AMIGA to isolate halos and subhalos, Halo to extract
properties for custom halos out of only one particle type, and stores
all output in HDF5 or plain text files.

\section install Compilation

\b AMATEUR needs following non-standard packages to be installed on
the system:

- \b GSL - the <em>GNU scientific library</em>. This open-source
package can be obtained at http://www.gnu.org/software/gsl.

- \b HDF5 - the <em>Hierarchical Data Format</em>.  This library has
been developed by NCSA and can be obtained at
http://hdf.ncsa.uiuc.edu/HDF5 .

- \b FFTW - the <em>Fastest Fourier Transform in the West</em>. This
open-source package can be obtained at http://www.fftw.org.

- \b AMIGA - the halo finder developed by A. Knebe.

The provided makefile is compatible with GNU-make, i.e. typing \b make
should build the executable <b>amateur</b>.  If your site does not
have GNU-make, get it, or write your own Makefile.

\section howtorun Running the code

\b AMATEUR is invoked using the following syntax

<b> ./amateur ID -{xamszht}</b> or <b> ./amateur ID </b>

where ID describes a snapshot given in Global and xamszt invoke the
different commands

- \b x - Extract snapshot data into a HDF5 file
- \b a - Analyze the given snapshot using AHFstep from AMIGA
- \b o - wrap AHF output files for use in MATLAB
- \b m - Run MergerTree
- \b s - Match subhalos with their hosthalo.
- \b z - Analyze the given snapshot using Halo
- \b t - Compute the tidal field of the mass distribution
- \b h - Run HaloTracker from AMIGA. This is a planned feature and not 
implemented yet.

If only the snapshot ID is given as input parameter to \b amateur, it
will run through the commands in the order listed above.
*/








/**
\page AMATEUR-Makefile  Makefile of AMATEUR

All compilation options are covered in the Makefile.

- \b all (or none) -  compiles all different parts
- \b part.o - compiles the corresponding part
- \b run{xaomszht} - compiles amateur and runs it for a small snapshot
- \b batch - compiles amateur and runs through all high-resolution snapshots
- \b doc - generates the documentation
- \b tags - generates tags for use in Emacs
- \b clean - removes all *.o in preparation for a future distribution
- \b tidy - remove all except sourcecode

\section modifications Modifications

Whenever the Makefile is modified, all source code files are
recompiled the next time <b>make</b> is invoked.

\section doc Documentation

The documentation for AMATEUR can be (re)generated by issuing

<b> make doc </b>

doxygen output is then available in ./doc.

*/
