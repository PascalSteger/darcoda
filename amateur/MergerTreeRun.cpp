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

/*!
 * \file MergerTreeRun.cpp
 * \brief Rewrite of MergerTree of AMIGA, implements -m of amateur
 */

#include <fstream>
#include <iostream>
#include <string>
#include "Global.h"

void RewriteParticleNumbers(unsigned isnap)
{
  writedebug(DEBUG, " * Extracting particle numbers... for #(halos) = ");
  std::ifstream fin;
  std::ofstream fout;
  fin.open((Folder("output") + Snaps(isnap) + "/ahf_out_particles").
	   c_str(), std::ios_base::in);
  fout.open((Folder("output") + Snaps(isnap) + "/ahf_out_particles2").
	    c_str(), std::ios_base::out);
  int Idummy = 0;
  int Ihalos = 0;
  int Icounter = 0;
  //  float Fdummy = 0.0f;

  fin >> Ihalos; //number of halos
  std::cout << Ihalos << std::flush;
  for (; Ihalos > 0; --Ihalos) {
    //    fin >> Idummy; //unknown
    fin >> Icounter; //number of particles in current halo
    fout << Icounter << std::endl;
    //fin >> Fdummy; //mass of current halo
    for (int i = 0; i < Icounter; ++i) {
      fin >> Idummy;// id of current particle
      fout << Idummy << std::endl;
      fin >> Idummy; // type of current particle
      //fin >> Fdummy; // mass of current particle
    }
  }

  fin.close();
  fout.close();
  donedebug(DEBUG);
}

void CreateMergeParameterFile(unsigned isnap)
{
  std::
    cout << " * Creating MergerTree parameter file..." <<
    Folder("amigabin");
  std::ofstream fout;
  fout.open((Folder("amigabin") + "parameters_MergerTree.txt").c_str(),
	    std::ios_base::out);
  // taken out _particles'''2'''
  fout << "2" << std::
    endl << Folder("output") << Snaps(isnap) << "/ahf_out_particles2"
       << std::endl;
  fout << Folder("output") << Snaps(isnap) << "/ahf_out_particles2" <<
    std::endl;
  fout << Folder("output") << Snaps(isnap) << "/ahf_out_mtree" << std::
    endl << std::endl;
  fout.close();
  writedone();
}

void RunMergerTree(void)
{
  writelndebug(DEBUG, " * Running MergerTree...");
  std::string now = "";
  now =
    Folder("amigabin") + "./MergerTree < " + Folder("amigabin") +
    "parameters_MergerTree.txt";
  RunCommand(now);
  //std::cin.get(); //wait for input
}

#include <map>
void RunOwnMergerTree(int isnap)
{
  writelndebug(DEBUG, " * Running own MergerTree..." );
  writedebug(DEBUG, " * * Extracting particle numbers... for #(halos) = " );
  std::ifstream fin;
  std::ofstream fout;
  fin.open((Folder("output") + Snaps(isnap) + "/ahf_out_particles").
	   c_str(), std::ios_base::in);
  fout.open((Folder("output") + Snaps(isnap) + "/ahf_out_particles2").
	    c_str(), std::ios_base::out);
  std::multimap < long, int >Particle_n_Halo;
  int Idummy = 0;
  long Ldummy = 0;
  int Ihalos = 0;
  int Icounter = 0;
  //  float Fdummy = 0.0f;

  fin >> Ihalos; //total number of halos
  std::cout << "Number of halos:" << Ihalos << std::endl << std::flush;
  int pop[Ihalos];
  for (int i = 0; i < Ihalos; ++i)
    pop[i] = 0;
  for (int i = 0; i < Ihalos; ++i) {
    //fin >> Idummy;		// which halo number
    fin >> Icounter;	// how many particles are inside this halo?
    pop[i] = Icounter;
    fout << Icounter << std::endl;
    std::cout <<Icounter<< " particles in halo " <<i<< std::endl;
    //fin >> Fdummy;//mass
    for (int j = 0; j < Icounter; ++j) {
      fin >> Ldummy;	// particle number
      fout << Ldummy << std::endl;
      Particle_n_Halo.insert(std::make_pair(Ldummy, i));
      fin >> Idummy;//type
      //fin >> Fdummy;//mass
    }
  }
  donedebug(DEBUG);

  writelndebug(DEBUG, " * Searching for shared particles..." );
  int *SH = new int[Ihalos * Ihalos];
  for (int m = 0; m < Ihalos; ++m)
    for (int n = 0; n < Ihalos; ++n)
      SH[m * Ihalos + n] = 0;
  writelndebug(DEBUG, " * * Starting map iterating..." );
  for (std::multimap < long, int >::iterator iter = Particle_n_Halo.begin();
       iter != Particle_n_Halo.end(); ++iter) {
    std::multimap < long, int >::iterator iter2 = iter, iter3;
    iter2 = Particle_n_Halo.find(iter->first);
    iter3 = Particle_n_Halo.upper_bound(iter->first);
    // std::cout << " * * * Checking for all halos the particle " << iter->
    // first << " belong to...";
    while (iter2 != iter3 ){ // Particle_n_Halo.end()) {
      int m = iter->second;
      int n = iter2->second;
      // std::cout << n << " ";
      SH[m * Ihalos + n]++;
      iter2++;
      // iter2 = Particle_n_Halo.find(iter->first);
    }
    // std::cout << std::endl;
  }
  std::cout << "done!" << std::endl;

  std::cout << " * Output..." << std::flush;
  std::string ff = Folder("output") + Snaps(isnap) + "/ahf_out_mtree";
  std::ofstream foutm(ff.c_str());
  for (int m = 0; m < Ihalos; ++m) {
    for (int n = 0; n < Ihalos; ++n) {
      if(SH[m * Ihalos + n] > 0)
	foutm << m << " " << pop[m] << " " << SH[m * Ihalos +
						 n] << " " << n << " " <<
	  pop[n] << std::endl;
    }
  }
  donedebug(DEBUG);

  writedebug(DEBUG, " * Closing files...");
  fin.close();
  fout.close();
  foutm.close();
  donedebug(DEBUG);

  delete[]SH;
  //std::cin.get(); //wait for input
}

void MoveMergeOut(unsigned isnap)
{
  writedebug(DEBUG, " * Moving ahf_out_ files...");
  std::string now = "mv ahf_out_* " + Folder("output") + Snaps(isnap);
  RunCommand(now);
  donedebug(DEBUG);
}

void MergerTreeRun(unsigned isnap)
{
  //RewriteParticleNumbers( isnap );
  //CreateMergeParameterFile( isnap ); //for AMIGA's MergerTree (buggy?)
  RunOwnMergerTree(isnap);
  //MoveMergeOut (isnap);
}
