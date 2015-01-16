/***************************************************************************
 *   Copyright (C) 2008 by Pascal Stephan Philipp Steger                   *
 *   psteger@phys.ethz.ch                                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, ma  02111-1307, USA.             *
 ***************************************************************************/

/**
 * \file ExtractSnap.cpp
 * \brief Reads out a Gadget2 snapshot file and stores the values in HDF5
 */

#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "HDFIO.h"
#include "ReadGadget.h"
#include "Global.h"

/**
 * analysis of one single snapshot file
 * @param isnap ID of snapshot formatted as a string of the form "XXX"
 * @return error code
 */
void 
ExtractSnap( unsigned isnap ) {
  std::string SnapshotString = Snaps( isnap );

  writelndebug(DEBUG, " * Getting particle numbers..." );
  unsigned npart[ 6 ], ntot;
  unsigned id_min[ 6 ], id_max[ 6 ];
  double Redshift, Boxsize;
  double masses[ 6 ];

  std::string fn_snap = Folder("snaps") + SnapshotString;
  get_particle_number ( fn_snap.c_str(), id_min, id_max, masses, npart, ntot, Redshift, Boxsize );
  donedebug(DEBUG);

  std::cout << " * offsets: ";
  unsigned offs[ 6 ];
  for ( int x = 0; x < 6; ++x ) {
    offs[ x ] = std::accumulate( &npart[ 0 ], &npart[ x ], 0 );
    std::cout << offs[ x ] << " ";
  }
  std::cout << std::endl;

  int n;
  unsigned ngas = npart[0];


  writedebug(DEBUG, "Opening snapshot file...");
  FILE *SnapshotFile = 0;
  SnapshotFile = fopen ( fn_snap.c_str(), "r" );
  donedebug(DEBUG);

  writedebug( DEBUG, "Opening HDF5 file...");
  std::string now = "mkdir -p " + Folder("output") + SnapshotString;
  RunCommand( now );
  std::string fn_out = Folder("output") + SnapshotString + "/snapshot.hdf5";
  HDFCreateFile( fn_out );
  donedebug(DEBUG);

  writedebug( DEBUG, "Writinig header information...");
  HDFCreateGroup( fn_out, "Header" );
  HDFWriteGroupAttribute( fn_out, "Header", "Boxsize", Boxsize );
  HDFWriteGroupAttribute( fn_out, "Header", "NumParticles", ntot );
  HDFWriteGroupAttribute( fn_out, "Header", "Redshift", Redshift );
  HDFWriteGroupAttribute( fn_out, "Header", "NumGas", npart[0] );
  HDFWriteGroupAttribute( fn_out, "Header", "NumDM_1", npart[1] );
  HDFWriteGroupAttribute( fn_out, "Header", "NumDM_2", npart[2] );
  HDFWriteGroupAttribute( fn_out, "Header", "NumDM_3", npart[3] );
  HDFWriteGroupAttribute( fn_out, "Header", "NumStars", npart[4] );
  HDFWriteGroupAttribute( fn_out, "Header", "NumBoundary", npart[5] );
  donedebug(DEBUG);

  // id
  std::cout << "I/O ID...";
  unsigned * id = new unsigned [ntot];
  n = read_gadget_1( ( unsigned * ) id,   "ID  ", SnapshotFile );
  std::vector<unsigned> rid( id, id + ntot );
  free(id);
  HDFWriteDataset( fn_out,       "Data_ID",   rid );
  {
    rid.clear();
    std::vector<unsigned> v(rid);
    rid.swap( v );
    std::cout << rid.capacity() << std::endl;
  }

  // pos
  std::cout << "I/O POS ...";
  float * pos = new float [3*ntot];
  n = read_gadget_3( ( float * )    pos,  "POS ", SnapshotFile );
  std::vector<float> rpos(pos, pos + 3 * ntot );
  free(pos);
  HDFWriteDatasetVector( fn_out, "Data_Pos",  rpos );
  {
    rpos.clear();
    std::vector<float> v(rpos);
    rpos.swap( v );
    std::cout << rpos.capacity() << std::endl;
  }
  // vel
  std::cout << "I/O VEL ...";
  float * vel = new float [3*ntot];
  n = read_gadget_3( ( float * )    vel,  "VEL ", SnapshotFile );
  std::vector<float> rvel( vel, vel + 3 * ntot );
  free(vel);
  HDFWriteDatasetVector( fn_out, "Data_Vel",  rvel );
  {
    rvel.clear();
    std::vector<float> v(rvel);
    rvel.swap( v );
    std::cout << rvel.capacity() << std::endl;
  }
  // mass
  std::cout << "I/O MASS...";
  float * mass = new float [ntot];
  n = read_gadget_1( ( float * )    mass, "MASS", SnapshotFile );
  std::vector<float> rmass( mass, mass + ntot );
  free(mass);
  HDFWriteDataset( fn_out,       "Data_Mass", rmass );
  {
    rmass.clear();
    std::vector<float> v(rmass);
    rmass.swap( v );
    std::cout << rmass.capacity() << std::endl;
  }
  // temp
  std::cout << "I/O TEMP...";
  float * temp = new float [ntot];
  n = read_gadget_1( ( float * )    temp, "TEMP", SnapshotFile );
  std::vector<float> rtemp(temp, temp + ngas);
  free(temp);
  HDFWriteDataset( fn_out,       "Data_Temp", rtemp );
  {
    rtemp.clear();
    std::vector<float> v(rtemp);
    rtemp.swap( v );
    std::cout << rtemp.capacity() << std::endl;
  }
  // u
  std::cout << "I/O U   ...";
  float * u = new float [ntot];
  n = read_gadget_1( ( float * )    u,    "U   ", SnapshotFile );
  std::vector<float> ru(u, u + ngas);
  free(u);
  HDFWriteDataset( fn_out,       "Data_U",    ru );
  {
    ru.clear();
    std::vector<float> v(ru);
    ru.swap( v );
    std::cout << ru.capacity() << std::endl;
  }
  // rho
  std::cout << "I/O RHO ...";
  float * rho = new float [ntot];
  n = read_gadget_1( ( float * )    rho,  "RHO ", SnapshotFile );
  std::vector<float> rrho(rho, rho + ngas);
  free(rho);
  HDFWriteDataset( fn_out,       "Data_Rho",  rrho );
  {
    rrho.clear();
    std::vector<float> v(rrho);
    rrho.swap( v );
    std::cout << rrho.capacity() << std::endl;
  }

  // check memory usage with top -i, while the program does not do anything
  //for(;;);


  writedone();
  writeln("");

  fclose ( SnapshotFile );
}
