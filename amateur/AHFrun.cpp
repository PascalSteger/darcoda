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
 * \file AHFrun.cpp
 * \brief Call AHFstep from AMIGA
 */

#include <fstream>
#include <iostream>
#include <string>

#include "Global.h"

/**
 * \brief Creates the parameter file for AHFstep
 * Creates the parameter file for AHFstep in output directory
 * \param isnap snapshot ID
 * \return void
 */
void 
CreateAHFParameterFile( unsigned isnap ){
  writedebug(DEBUG," * Creating AHFstep parameter file...");
  std::ofstream fout;
  fout.open ( (Folder("parameters") + "parameters_AHF_" + Snaps(isnap)).c_str(), std::ios_base::out );
  fout << Folder("snaps") << Snaps(isnap) 
       << " "
       << "60" //file type, here GADGET(2), single file
       << " "
       << "1"  //number of CPUs used to read in data
       << std::endl
       << "ahf_out_" << std::endl
       << "8"  // base grid (1dimensional)
       << std::endl
       << "4"  // refinement criterion on base grid
       << std::endl
       << "4"  // refinement criterion on refined grid
       << std::endl
       << "0" <<std::endl
       << "0" <<std::endl
       << "0" <<std::endl
       << "0" <<std::endl;
  fout.close();
  donedebug(DEBUG);
}

/**
 * \brief calls AHFstep
 * runs AHFstep parameterfile
 * @param isnap snapshot ID
 * @return void
 */
void 
RunAHFstep( unsigned isnap ){
  writelndebug(DEBUG, " * Running AHFstep..." );
  std::string now = "AHFstep " + Folder("parameters") + "parameters_AHF_" + Snaps(isnap);
  std::cout << now << std::endl;
  RunCommand( now );
  donedebug(DEBUG);
}

/**
 * \brief move output files of AHFstep
 * Moves output files from current directory to output/[snap-string]
 * @param isnap snapshot ID
 * @return void
 */
void 
MoveAHFout(unsigned isnap){
  writedebug(DEBUG, " * Moving output files...");
  std::string now;
  std::string ahfout = Folder( "output" ) + Snaps(isnap) + "/";
  now = "mkdir -p " + ahfout;
  RunCommand(now);
  now = "mv ahf_out_* " + ahfout;
  RunCommand(now);
  donedebug(DEBUG);
}

/**
 * \brief Renames AHFstep output files
 * Renames output files of AHFstep using regexpressions with bash
 * @param isnap snapshot ID
 * @return void
 */
void 
RenameAHFout(unsigned isnap){
  writedebug(DEBUG, " * Renaming output files...");
  std::string now = "";
  std::string ahfout = Folder( "output" ) + Snaps(isnap) + "/";

  now = "mv " + ahfout + "ahf*particles " + ahfout + "ahf_out_particles";
  RunCommand(now);
  now = "mv " + ahfout + "ahf*halos " + ahfout  + "ahf_out_halos";
  RunCommand(now);
  //  now = "mv " + ahfout + "ahf*mass " + ahfout + "ahf_out_mass";
  //  RunCommand(now);
  now = "mv " + ahfout + "ahf*profiles " + ahfout + "ahf_out_profiles";
  RunCommand(now);
  now = "mv " + ahfout + "ahf*substructure " + ahfout + "ahf_out_substructure";
  RunCommand(now);
  now = "mv " + ahfout + "ahf*centres " + ahfout + "ahf_out_centres";
  RunCommand(now);
  now = "tail -n+2 " + ahfout + "ahf_out_halos > " + ahfout  + "ahf_out_halos2";
  RunCommand(now);

  donedebug(DEBUG);
}

/**
 * \brief Controls execution of AHFstep
 */
/**
 * Creates the parameter file,
 * invokes AHFstep,
 * moves output files
 * and renames them
 * @param isnap snapshot ID
 * @return void
 */
void 
AHFrun( unsigned isnap ){

  CreateAHFParameterFile(isnap);

  RunAHFstep(isnap);

  MoveAHFout(isnap);
  RenameAHFout(isnap);
}
