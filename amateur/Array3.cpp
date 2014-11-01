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
 * \file Array3.cpp
 * \brief threedimensional array
 */

#include <algorithm> // for max
#include <iostream>
#include "Array3.h"

Array3::Array3( void ){
  data_ = new float();
  nx_ = 0;
  ny_ = 0;
  nz_ = 0;
  RealContainer_ = false;
}

Array3::Array3( float *dataptr, int nx, int ny, int nz ){
  cout << " in Array3 constructor with dataptr..." << endl;
  data_ = dataptr;
  nx_ = nx;
  ny_ = ny;
  nz_ = nz;
  RealContainer_ = true;
}

Array3::Array3( int nx, int ny, int nz ){
  std::cout << " in Array3 constructor..." << std::endl;
  nx_ = nx;
  ny_ = ny;
  nz_ = nz;
  data_ = new float[ nx_*ny_*nz_ ];
  RealContainer_ = true;
  for ( long long i = 0; i < nx_*ny_*nz_; ++i ){
    data_[ i ] = 0.0f;
  }
}

Array3::~Array3(){
  if ( RealContainer_ ){
    delete[] data_;
  }
}

Array3*
Array3::operator=(Array3* a){
  nx_ = a->getNx();
  ny_ = a->getNy();
  nz_ = a->getNz();
  data_ = a->getData();
  RealContainer_ = a->getRealContainer();
  return this;
}

//float 
//Array3::operator() ( int ix, int iy, int iz ) {
//  long long i =  ( ix * ny_ + iy ) * nz_ + iz;
//  return data_[i];
//}

float
Array3::getValue(int ix, int iy, int iz){
  return data_[(ix*ny_+iy)*nz_+iz];
}

void
Array3::setValue(int ix, int iy, int iz, float val){
  long long i =  ( ix * ny_ + iy ) * nz_ + iz;
  data_[i] = val;
}

void
Array3::addValue(int ix, int iy, int iz, float val){
  long long i =  ( ix * ny_ + iy ) * nz_ + iz;
  data_[i] += val;
} 

int
Array3::get_size( void ) {
  return std::max( std::max( nx_, ny_ ), nz_ );
}

void 
Array3::recreate( int nx, int ny, int nz ) {
  if ( RealContainer_ ){
    delete[] data_;
  } else {
    RealContainer_ = true;
    nx_ = nx;
    ny_ = ny;
    nz_ = nz;
    data_ = new float[ nx_ * ny_ * nz_ ];
	
    for ( long long i = 0; i < nx*ny*nz; ++i ){
      data_[ i ] = 0.0f;
    }
  }
}

