/**************************************************************************
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
 * \file HaloTrackerRun.cpp
 * \brief Implementation of -h, EXPERIMENTAL!
 */

#include <fstream>
#include <iostream>
#include <string>
#include "Global.h"

void 
WrapperParticleNumbers(unsigned isnap){
  writedebug(DEBUG," * Extracting particle numbers...");
  std::ifstream fin;
  std::ofstream fout;
  fin.open ((Folder("output")+Snaps(isnap)+"/ahf_out_"+Snaps(isnap)+"_particles").c_str(), std::ios_base::in);
  fout.open((Folder("output")+Snaps(isnap)+"/ahf_out_"+Snaps(isnap)+"_particles2").c_str(), std::ios_base::out );
  int Idummy = 0;
  int Ihalos = 0;
  int Icounter = 0;
  float Fdummy = 0.0f;

  fin >> Ihalos;
  for(;Ihalos>0;--Ihalos){
    fin >> Idummy;
    fin >> Icounter;
    fout << Icounter << std::endl;
    fin >> Fdummy;
    for(int i=0; i < Icounter; ++i){
      fin >> Idummy;
      fout << Idummy << std::endl;
      fin >> Idummy;
      fin >> Fdummy;
    }
  }

  fin.close();
  fout.close();
  donedebug(DEBUG);
}

void 
CreateHaloTrackerParameterFile( unsigned isnap ){
  std::cout << " * Creating HaloTracker parameter file..." << Folder("amigabin");
  std::ofstream fout;
  fout.open ( (Folder("amigabin")+"parameters_HaloTracker.txt").c_str(), std::ios_base::out );
  fout << "0" << std::endl //HaloID
       << "ahf_out_" << std::endl //ahf prefix
       << "2" << std::endl// no. snapshots
       << Folder("snaps") << Snaps(isnap) << "/snap_" << Snaps(isnap) << std::endl
       << Folder("snaps") << Snaps(isnap) << "/snap_" << Snaps(isnap) << std::endl
       << std::endl;
  fout.close();
  writedone();
}

void 
RunHaloTracker(void){
  writelndebug(DEBUG," * Running HaloTracker...");
  std::string now = "";
  now = Folder("amigabin") + "./HaloTracker < " + Folder("amigabin") + "parameters_HaloTracker.txt";
  RunCommand(now);
  //std::cin.get(); //wait for input
}

void 
MoveHaloTrackerOutput(unsigned isnap){
  writedebug(DEBUG, " * Moving ahf_out_ files...");
  std::string now = "mv ahf_out_* " + Folder("output") + Snaps(isnap);
  RunCommand( now );
  donedebug(DEBUG);
}

void 
HaloTrackerRun( unsigned isnap ){
  // already done in MergeTree
  WrapperParticleNumbers(isnap);
  CreateHaloTrackerParameterFile(isnap);
  RunHaloTracker();
  MoveHaloTrackerOutput(isnap);
}

