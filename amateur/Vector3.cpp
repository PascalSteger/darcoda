/***************************************************************************
 *   Copyright (C) 2008 by Pascal Stephan Philipp Steger                   *
 *   psteger@student.ethz.ch                                               *
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
 * \file Vector3.cpp
 * \brief Vector3 implementations
 */

#include <math.h>
#include <iterator>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <list>
#include <vector>
#include <iterator>
#include <cassert>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "Global.h"
#include "Vector3.h"

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

Vector3::Vector3( void ) {
  data_ = new real[ 3 ];
  for ( unsigned i = 0; i < 3; i++ )
    data_[ i ] = 0.0;
}

Vector3::Vector3( const Vector3& v ) {
  data_ = new real[ 3 ];
  for ( unsigned i = 0; i < 3; i++ )
    data_[ i ] = v.data_[ i ];
}

Vector3::Vector3( real x){
  data_ = new real[ 3 ];
  for ( unsigned i = 0; i < 3; i++ )
    data_[ i ] = x;
}

Vector3::Vector3( real x, real y, real z ) {
  data_ = new real[ 3 ];
  data_[ 0 ] = x;
  data_[ 1 ] = y;
  data_[ 2 ] = z;
}

template < typename real_t > 
Vector3::Vector3( const real_t* x0 ) {
  data_ = new real[ 3 ];
  for ( unsigned i = 0; i < 3; i++ )
    data_[ i ] = x0[ i ];
}

Vector3::~Vector3() {
  delete[] data_;
}

//inline
Vector3& 
Vector3::operator=( const Vector3& b ) {
  for ( unsigned i = 0; i < 3; i++ )
    data_[ i ] = b.data_[ i ];
  return *this;
}

//inline 
Vector3& 
Vector3::operator=( const real& x ) {
  for ( unsigned i = 0; i < 3; i++ )
    data_[ i ] = x;
  return *this;
}

//inline 
Vector3& 
Vector3::operator-( void ) {
  for ( unsigned i = 0; i < 3; i++ )
    data_[ i ] = -data_[ i ];
  return *this;
}

//inline
Vector3 
Vector3::operator+( const Vector3& x ) const {
  Vector3 temp;
  for ( unsigned i = 0; i < 3; i++ )
    temp.data_[ i ] = data_[ i ] + x.data_[ i ];
  return temp;
}

//inline 
Vector3 
Vector3::operator-( const Vector3& x ) const {
  Vector3 temp;
  for ( unsigned i = 0; i < 3; i++ )
    temp.data_[ i ] = data_[ i ] - x.data_[ i ];
  return temp;
}

//inline 
Vector3 
Vector3::operator*( const real& x ) {
  Vector3 temp;
  for ( unsigned i = 0; i < 3; i++ )
    temp.data_[ i ] = data_[ i ] * x;
  return temp;
}

//inline 
Vector3 
Vector3::operator/( const real& x ) {
  Vector3 temp;
  for ( unsigned i = 0; i < 3; i++ )
    temp.data_[ i ] = data_[ i ] / x;
  return temp;
}

//inline 
Vector3& 
Vector3::operator+=( const Vector3& x ) {
  for ( unsigned i = 0; i < 3; i++ )
    data_[ i ] += x.data_[ i ];
  return *this;
}

//inline 
Vector3& 
Vector3::operator-=( const Vector3& x ) {
  for ( unsigned i = 0; i < 3; i++ )
    data_[ i ] -= x.data_[ i ];
  return *this;
}

//inline 
Vector3& 
Vector3::operator*=( const real& x ) {
  for ( unsigned i = 0; i < 3; i++ )
    data_[ i ] *= x;
  return *this;
}

//inline 
Vector3& 
Vector3::operator/=( const real& x ) {
  for ( unsigned i = 0; i < 3; i++ )
    data_[ i ] /= x;
  return *this;
}

//inline 
real& 
Vector3::operator() ( const unsigned i ) {
  return data_[ i ];
}

// dot product
//inline 
real 
Vector3::operator*( const Vector3& vb ) const {
  real res = 0.0;
  for ( unsigned i = 0; i < 3; i++ )
    res += data_[ i ] * vb.data_[ i ];
  return res;
}

// cross product
//inline 
Vector3 
Vector3::operator^( const Vector3& vb ) const {
  Vector3 temp;

  temp.data_[ 0 ] = data_[ 1 ] * vb.data_[ 2 ] - data_[ 2 ] * vb.data_[ 1 ];
  temp.data_[ 1 ] = data_[ 2 ] * vb.data_[ 0 ] - data_[ 0 ] * vb.data_[ 2 ];
  temp.data_[ 2 ] = data_[ 0 ] * vb.data_[ 1 ] - data_[ 1 ] * vb.data_[ 0 ];

  return temp;
}

//inline 
real 
Vector3::norm( void ) const {
  return sqrt( ( *this ) * ( *this ) );
}

//inline
real 
Vector3::norm2( void ) const {
  return ( *this ) * ( *this );
}

//inline 
Vector3& 
Vector3::normalize( void ) {
  if ( norm2() < 1e-10 )
    return *this;
  *this /= norm();
  return *this;
}

// 	inline Vector& moduloBox( const Vector& vm ) {
// 		for( Identifier i=0; i<N; ++i ) {
// 			if(data_[i] < -vm[i]/2) data_[i] += vm[i];
// 			if(data_[i] > vm[i]/2) data_[i] -= vm[i];
// 		}
// 	}

std::string
Vector3::toString(){
  std::ostringstream out;
  for(int i=0; i<3; ++i){
    out << data_[i] << " ";
  }
  return out.str();
}

// output
std::ostream& 
operator<<( std::ostream& s, Vector3 v ) {
  s << v( 0 );
  for ( unsigned i = 1; i < 3 && s << " "; ++i )
    s << v( i );
  return s;
}
