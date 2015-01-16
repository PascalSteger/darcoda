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
 * \file ScalarField.h
 * \brief Phi
 */

#ifndef SCALARFIELD_H
#define SCALARFIELD_H

#include <vector>
#include "Matrix33.h"
#include "Array3.h"
//#include "Stencil.h"

class ScalarField {
 private:
  Array3 phi_;
  //  std::vector<ScalarField> dy;
  float dx_;
  unsigned nx_, ny_, nz_;
  int nStencil_, nStencil2_;
  //Stencil pStencil1_, pStencil2_;
  float *pStencil1_, *pStencil2_;

  void diff_primitive( ScalarField*, bool bStagger = false );

 public:
  ScalarField( void );
  ScalarField( unsigned nx, unsigned ny, unsigned nz, 
	       float dx, int nstencil );
  ~ScalarField( void );
  ScalarField* operator=(ScalarField* s);
  void recreate( unsigned nx, unsigned ny, unsigned nz, 
		 float dx, int nstencil );

  Array3* getPhi( void );
  void setPhi(Array3* phi);
  int getNx( void );
  int getNy( void );
  int getNz( void );
  float getDx( void );
  int getNStencil( void );
  int getNStencil2(void);
  void test( void );

  //float operator()(int ix, int iy, int iz);
  float getValue(int ix, int iy, int iz);
  void setValue(int ix, int iy, int iz, float val);

  void ChooseStencil( void );

  void put_CIC( unsigned npart, float *coord, float *wpar );
  float get_CIC( float coordx, float coordy, float coordz );

  //  void grad( ScalarField *dy );
  //  void diff( VectorField* dy );

  void lookup( float *Coord, unsigned nCoord, 
	       std::vector< float >& values );
  void lookup( float *Coord, unsigned nCoord, float ll, unsigned N, 
	       std::vector< float >& values );

  float diff1( int ix, int iy, int iz, int i );
  float diff2( int ix, int iy, int iz, int i, int j );
  void diff2atCIC( float *Coord, unsigned nPoints, 
		   std::vector< Matrix33 > &t_ij );

  template < typename T >
    T get_NGP( T *coord ) {
    int ix = ( int ) coord[ 0 ];
    int iy = ( int ) coord[ 1 ];
    int iz = ( int ) coord[ 2 ];
    return this->getValue( ix, iy, iz );
  }

};

#endif
