/***************************************************************************
 *   Copyright (C) 2010 by Pascal Stephan Philipp Steger                   *
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
 * \file AnalyzeSnap.cpp
 * \brief Analyzes a given snapshot with Halo, for all particle types
 */

////////////////////////////////////////
// inclusions
////////////////////////////////////////

#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include "Global.h"
#include "Vector3.h"
#include "Matrix33.h"
#include "Eigensystem3.h"
#include "ReadGadget.h"
#include "Simulation.h"
#include "Halo.h"
#include "HDFIO.h"
#include "HResults.h"
#include "AResults.h"
#include "Output.h"

////////////////////////////////////////
// overall definitions
////////////////////////////////////////

typedef std::vector<double> VecDouble;
typedef std::vector<std::vector <double> > VecVecDouble;

const int HaloTypeTotal = 5;
const int ParticleTypeTotal = 6;

/**
 * \brief little/big endian conversion
 * conversion between little/big endian
 * \param v
 * \return templated T value, in other byteorder
 */
template < typename T > T byteorder(T v) {
  T val;
  (reinterpret_cast < unsigned char *>(&val))[3] =
    (reinterpret_cast < unsigned char *>(&v))[0];
  (reinterpret_cast < unsigned char *>(&val))[2] =
    (reinterpret_cast < unsigned char *>(&v))[1];
  (reinterpret_cast < unsigned char *>(&val))[1] =
    (reinterpret_cast < unsigned char *>(&v))[2];
  (reinterpret_cast < unsigned char *>(&val))[0] =
    (reinterpret_cast < unsigned char *>(&v))[3];
  return val;
}

/**
 * \brief Return offsets of particle types inside snapshot.
 * Get the particle ID offsets of the different particle types
 * \param ParticleIDOffsets: unsigned*; where the offsets are stored
 * \param npart: unsigned*; counts of particles types
 * \return void
 */
void
GetOffsets(unsigned* ParticleIDOffsets, unsigned* npart)
{
  std::cout << " in GetOffsets..." << std::endl;
  std::cout << "offsets: ";
  for (int x = 0; x < ParticleTypeTotal; ++x) {
    ParticleIDOffsets[x] = std::accumulate(&npart[0], &npart[x], 0);
    std::cout << ParticleIDOffsets[x] << " ";
  }
  std::cout << std::endl;
}

/**
 * \brief Read snapshot data
 * Fill in all information from snapshot data file (id, pos, vel, mass)
 * \param isnap: int; snapshot number
 * \param TotalNOParticles: int
 * \param id: is int &
 * \param pos: is double & (3d)
 * \param vel: is double & (3d)
 * \param mass: is double & (1d)
 * \return void
 */
void 
readDataBlocks(int isnap, int TotalNOParticles,
	       int* & id, real* & pos, real* & vel, real* & mass){
  std::cout << " in ReadDataBlocks..." << std::endl;
  FILE *SnapshotFile = 0;
  std::string SnapshotFilename = Folder("snaps")+Snaps(isnap);
  SnapshotFile = fopen(SnapshotFilename.c_str(), "r");
  int n = 0;
  n = read_gadget_1((int *) id, "ID  ", SnapshotFile);
  
  float *ftemp = new float[3 * TotalNOParticles];
  n = read_gadget_3(ftemp, "POS ", SnapshotFile);
  for (int ii = 0; ii < 3 * n; ++ii)
    pos[ii] = ftemp[ii];
  
  n = read_gadget_3(ftemp, "VEL ", SnapshotFile);
  for (int ii = 0; ii < 3 * n; ++ii)
    vel[ii] = ftemp[ii];
  
  n = read_gadget_1(ftemp, "MASS", SnapshotFile);
  for (int ii = 0; ii < n; ++ii)
    mass[ii] = ftemp[ii];
  
  delete[] ftemp;
  fclose(SnapshotFile);
}

/**
 * \brief Open the ahf_out_particles file
 * take ifstream and open particles file with it
 * \param isnap: unsigned of which snapshot, see SystemCharacteristics
 * \param ParticlesFile: ifstream by ref
 * \return void
 */
void 
openParticlesFile(unsigned isnap, std::ifstream & ParticlesFile){
  std::cout << " in openParticlesFile..." << std::endl;
  std::string SnapshotParticlesFilename =
    Folder("output") + Snaps(isnap) + "/ahf_out_particles";
  ParticlesFile.open(SnapshotParticlesFilename.c_str(),
		     std::ios_base::in | std::ios_base::binary);
}


/**
 * \brief Read a halo full of snapshot particles from ahf_out_particles
 * Read in a slab of HaloPopulation particles
 * \param ParticlesFile: ifstream of ahf_out_particles
 * \param HaloPopulation: unsigned int of number of particles in this halo
 * \param pid: holds the output, all identifiers
 * \param ptype: holds the output, all types, in accordance with ParticleIndexID
 */
void
ReadSnapIdType(std::ifstream & ParticlesFile, 
	       unsigned HaloPopulation,
	       unsigned* & pid, int* & ptype){
  std::cout << " in ReadSnapIdType..." << std::endl;
  unsigned dummy_uns;
  int dummy_int;
  for (unsigned i = 0; i < HaloPopulation; ++i) {
    ParticlesFile >> dummy_uns;
    ParticlesFile >> dummy_int;
    pid[i] = dummy_uns;
    ptype[i] = dummy_int;
  }
}


/**
 * \brief preparing halo
 * extract new lists from pos, vel, mass based on id (AHF_particles)
 * \param ParticlesFile: std::ifstream& of ahf_out_particles
 * \param HaloPopulation: unsigned; number of particles total
 * \param DoubleIndexMapi: map with reverse lookup possibility
 * \param pos, vel, mass: const real *&, original data
 * \param hpos, hvel, hmass: real *&, data for halo
 * \param htype: particle type for halo
 * \return void
 */
void
GenerateHalo(std::ifstream & ParticlesFile,
	     unsigned HaloPopulation,
	     const std::map<int,unsigned>& DoubleIndexMapi,
	     real* pos,  real* vel,  real* mass,
	     real* hpos, real* hvel, real* hmass, int* htype){
  std::cout << " in GenerateHalo..." << std::endl;
  std::cout << " * * Reading id, type from _particles..." << std::endl;
  unsigned *pid   = new unsigned[HaloPopulation];
  int      *ptype = new int[HaloPopulation];
  ReadSnapIdType(ParticlesFile, HaloPopulation, pid, ptype);
  
  std::cout << " * * copying pos,vel,mass,type..." << std::endl;
  for (unsigned k = 0; k < HaloPopulation; ++k) {
    unsigned offset  = DoubleIndexMapi.find(pid[k])->second;
    
    hmass[k] = mass[offset];
    htype[k] = ptype[k];
    for (int i = 0; i < 3; ++i) {
      hpos[3*k+i] = pos[3*offset+i];
      hvel[3*k+i] = vel[3*offset+i];
    }
  }
  delete[] pid;
  delete[] ptype;
}

/**
 * \brief generate doubly indexed map
 * DIM for particle ID to occurence position
 * \param TotalNOParticles; unsigned
 * \param id: int* &; identification int for all particles in simulation
 * \param DIM: std::map<int,unsigned>&; output
 * \return void
 */
void
GenerateDIM(unsigned TotalNOParticles,
	    int* & id,
	    std::map<int,unsigned>& DIM){
  std::cout << " in GenerateDoubleIndexMap..." << std::endl;
  for (unsigned j = 0; j < TotalNOParticles; ++j)
    DIM.insert(std::make_pair(id[j], j));
}

/**
 * \brief Analysis of one single snapshot file
 * Analyzes a single snapshot file.
 * \param isnap unsigned snapshot ID
 * \return void
 */
void 
AnalyzeSnap(unsigned isnap)
{
  ////////////////////////////////////////////////////////////////////
  // preparation
  ////////////////////////////////////////////////////////////////////
  
  unsigned TotalNOParticles, npart[ParticleTypeTotal];
  unsigned id_min[ParticleTypeTotal], id_max[ParticleTypeTotal];
  double Z, BL;
  double masses[ParticleTypeTotal];
  get_particle_number((Folder("snaps")+Snaps(isnap)).c_str(), 
		      id_min, id_max, masses,
		      npart, TotalNOParticles, Z, BL);
  
  std::ifstream ParticlesFile;
  openParticlesFile(isnap, ParticlesFile);
  
  unsigned nhalo = 0;
  ParticlesFile >> nhalo;
  std::cout << nhalo << " halos totally." << std::endl;
  
  int  *id   = new int [TotalNOParticles];
  real *pos  = new real[3*TotalNOParticles];
  real *vel  = new real[3*TotalNOParticles];
  real *mass = new real[TotalNOParticles];	//, * temp, * u, * rho;
  readDataBlocks(isnap,
		 TotalNOParticles, 
		 id, pos, vel, mass);
  
  std::map<int,unsigned> DoubleIndexMap;
  GenerateDIM(TotalNOParticles, id, DoubleIndexMap);
  
  real msimtot=0.0;
  for(unsigned i=0; i<TotalNOParticles; ++i){
    msimtot += mass[i];
  }
  Simulation sim;
  sim.setIsnap(isnap);
  sim.setZsnap(Z);
  sim.setTarkinOmega();
  sim.setMass(msimtot);
  sim.setBoxlength(BL);
  sim.setNhalo(nhalo);
  Output out(sim);
  /* big loop over all halos */
  for (unsigned hit = 0;
       hit < nhalo;
       ++hit) {
    std::cout << "Halo " << hit << std::endl;
    std::cout << "+++++++++++++" << std::endl;
    
    unsigned npart;
    ParticlesFile >> npart;
    std::cout << "Halo population = " << npart << std::endl;

    int  *htype = new int [npart];
    real *hpos  = new real[3*npart];
    real *hvel  = new real[3*npart];
    real *hmass = new real[npart];
    GenerateHalo(ParticlesFile, npart, DoubleIndexMap,
		 &pos[0], &vel[0], &mass[0],
		 &hpos[0],&hvel[0],&hmass[0],&htype[0]);
    
    std::cout << " * * Creating halo from " << npart;
    std::cout << " particles" << std::endl;
    
    Halo halo;
    halo = Halo(sim, hpos, hvel, mass, htype, npart);
    
#if 0 // if true: exclude unbound particles with procedure from Halo
    writelndebug(DEBUG, " * * Cleaning unbound particles...");
    std::vector < bool > exc(npart);
    halo.CleanUnboundExc(Boxsize, exc);
    std::vector < bool > excs(npart);
    halo.CleanUnboundSimple(Boxsize, excs);
    std::cout << std::endl;
#else   // change all following calculations such that they have an exclusion list available
#warning No additional exclusion of unbound particles!
#endif
    
    /**************************************************
     * properties of halo consisting of all particles *
     **************************************************/
    real Mtot = halo.getMtot();
    real Rtot = halo.getRtot();
    Vector3 xcm = halo.getXCM();
    Vector3 vcm = halo.getVCM();
    Vector3 xmb = halo.getXMB();
    Vector3 vmb = halo.getVMB();
    real Rvir = halo.getRvir();
    real Mvir = halo.getMvir();
    std::cout << " * * Total mass:           " << Mtot << std::endl;
    std::cout << " * * maximal radius:       " << Rtot << std::endl;
    std::cout << " * * Center of mass:       " << xcm << std::endl;
    std::cout << " * * Center velocity:      " << vcm << std::endl;
    std::cout << " * * pos_mostbound:        " << xmb << std::endl;
    std::cout << " * * vel_mostbound:        " << vmb << std::endl;
    std::cout << " * * virial properties:    ";
    std::cout << "Rvir = " << Rvir << ", Mvir = " << Mvir << std::endl;
    std::cout << std::endl;
    
    out.fillHRglobal(npart,Mtot,Rtot,Mvir,Rvir,xcm,vcm,xmb,vmb);
    
    /**
     * properties depending on halo type
     * all:   1 1 1 1  15 dec
     * dm:    0 0 1 1   3 dec
     * gas:   0 1 0 0   4 dec
     * stars: 1 0 0 0   8 dec
     * bary:  1 1 0 0  12 dec
     */
    
    int ht[5] = {15, 3, 4, 8, 12};
    std::string hts[5] = {"all", "dm", "gas", "stars", "visible"};
    
    for(int type = 0; type<5; ++type){
      std::cout << " halo type: " << hts[type] << std::endl;
      
      std::cout << " * * Number of particles:          ";
      int npart = halo.calcNpart(ht[type]);
      std::cout << npart << std::endl;
      
      std::cout << " * * Mass of particles:            ";
      real mpart = halo.calcMpart(ht[type]);
      std::cout << mpart << std::endl;
      
      std::cout << " * * Angular momentum:             ";
      Vector3 J = halo.calcJ(ht[type]);
      std::cout << J << std::endl;
      
      std::cout << " * * Lambda':                      ";
      real lambda = halo.calcLambda(ht[type]);
      std::cout << lambda << std::endl;
      
      std::cout << " * * sigma:                        ";
      real sigma = halo.calcSigma(ht[type]);
      std::cout << sigma << std::endl;
      
      std::cout << " * * Axes...                       ";
      Eigensystem3 E  = halo.calcShape(ht[type]);
      std::cout << E.getEval().at(0) << ",";
      std::cout << E.getEval().at(1) << ",";
      std::cout << E.getEval().at(2) << std::endl;
      
      out.fillHRtype(type, npart, mpart, lambda, sigma, J, E);
    }
    
    std::cout << " * * Deleting arrays from halo..." << std::endl;
    delete[] hpos;
    delete[] hvel;
    delete[] hmass;
    delete[] htype;
    
    out.appendHR();
    
    //      int goon;
    //      std::cin >> goon;
  } // next Halo
  
  /* cleaning */
  delete[] id;
  delete[] pos;
  delete[] vel;
  delete[] mass;
  ParticlesFile.close();
  
  out.writeHR();
  out.writeHRstripped();
  
  writeln("Success.");
  return;
}
