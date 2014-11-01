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
 * \file Output.h
 * \brief filling output classes with values and writing them to files
 */

#ifndef OUTPUT_H
#define OUTPUT_H

#include <vector>
#include "Global.h"
#include "Simulation.h"
#include "AResults.h"
#include "HResults.h"

class Output{
 private:
  Simulation sim_;
  
  std::vector<AResults> ar_;
  std::vector<HResults> hr_;
  HResults tmp_;
  
 public:
  Output();
  Output( Simulation sim );
  ~Output();

  void getAR( void );
  void writeAR( void );

  void fillHRglobal( int Ntot, real Mtot, real Rtot, real Mvir, real Rvir,
		     Vector3 xcm, Vector3 vcm, Vector3 xmb, Vector3 vmb );
  void fillHRtype( int type, int npart, real mpart, real lambda, real sigma,
		   Vector3 J, Eigensystem3 E );
  void appendHR( void );
  void writeHR( void );
  void writeHRstripped( void );
};

#endif
