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
 * \file TidalField.h
 * \brief Class declarations for invocation of amateur .. -t
 */

#ifndef TIDALFIELD_H
#define TIDALFIELD_H

#define EPS_ACC 1.0e-6

#include <math.h>
#include <iostream>
#include "Vector3.h"
#include "Matrix33.h"
#include "Eigensystem3.h"
#include "Array3.h"
#include "GreensFunction.h"
#include "ScalarField.h"

void TidalField( unsigned isnap, int N, double Radius );
void TidalField( std::string snaphdf5, unsigned N, double Radius, 
		 std::string extra, std::string tensorhdf5 );
void StoreTensorField( std::string filename, unsigned nPoints,
		       const std::vector<float>& coordinates,
		       const std::vector< Matrix33 > &aTensorField,
		       float SmoothingMassScale );

#endif
