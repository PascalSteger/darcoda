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
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/**
 * \file Vector3.h
 * \brief All publicly known functions of Vector3.cpp
 * A three-dimensional vector (a,b,c)
 */

#ifndef VECTOR3_H
#define VECTOR3_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <list>
#include <vector>
#include <iterator>
#include <cassert>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "Global.h"

#ifdef LONGIDS
typedef unsigned long Identifier;
typedef unsigned long Size;
#define HDF_UINT H5T_NATIVE_UINT64
#define MPI_UINT MPI_UNSIGNED_LONG

#else
typedef unsigned Identifier;
typedef unsigned Size;
#define HDF_UINT H5T_NATIVE_UINT
#define MPI_UINT MPI_UNSIGNED
#endif

#define HDF_FLOAT H5T_NATIVE_FLOAT
typedef float Scalar;
typedef float Weight;

#define FIX(x) (Identifier)(x+0.5)
// implements a vector of fixed length 3
class Vector3 {
 protected:
  real* data_;
 public:
  Vector3( void );
  Vector3( const Vector3& v );
  Vector3( real x);
  Vector3( real x, real y, real z );

  template < typename real_t > Vector3( const real_t* x0 );

  ~Vector3();

  Vector3& operator=( const Vector3& b );
  Vector3& operator=( const real& x );
  Vector3& operator-( void );
  Vector3 operator+( const Vector3& x ) const;
  Vector3 operator-( const Vector3& x ) const;
  Vector3 operator*( const real& x );
  Vector3 operator/( const real& x );
  Vector3& operator+=( const Vector3& x );
  Vector3& operator-=( const Vector3& x );
  Vector3& operator*=( const real& x );
  Vector3& operator/=( const real& x );
  real& operator() ( const unsigned i );
  //dot product
  real operator*( const Vector3& vb ) const;
  // cross product
  Vector3 operator^( const Vector3& vb ) const;
  real norm( void ) const;
  real norm2( void ) const;
  Vector3& normalize( void );

  std::string toString( void );
};

std::ostream& 
operator<<( std::ostream& s, Vector3 v );

#endif
