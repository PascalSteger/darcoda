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
 * \file Matrix33.cpp
 * \brief Matrix33 implementations
 */
#include <math.h>
#include <iterator>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <list>
#include <vector>
#include <iterator>
#include <cassert>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "Vector3.h"
#include "Matrix33.h"

Matrix33::Matrix33() {
  data_ = new double[ 3 * 3 ];
}

Matrix33::Matrix33( const Matrix33& m ) {
  data_ = new double[ 3 * 3 ];
  for ( int i = 0; i < 3*3; i++ )
    data_[ i ] = m.data_[ i ];
}

Matrix33::~Matrix33() {
  delete[] data_;
};

//inline 
double& 
Matrix33::operator() ( int row, int col ) {
  return data_[ col + row * 3 ];
}

//inline 
double 
Matrix33::operator() ( int row, int col ) const {
  return data_[ col + row * 3 ];
}

//inline 
Vector3 
Matrix33::operator() ( int col ) const {
  Vector3 v;
  for ( unsigned i = 0; i < 3; ++i )
    v( i ) = ( *this ) ( i, col );
  return v;
}

//inline 
Matrix33& 
Matrix33::operator= ( Matrix33 const &M ) {
  int i;
  for ( i = 0;i < 3*3;i++ )
    data_[ i ] = M.data_[ i ];
  return *this;
}

//inline 
Matrix33& 
Matrix33::operator= ( const double& x ) {
  int i;
  for ( i = 0;i < 3*3;i++ )
    data_[ i ] = x;
  return *this;
}

//inline 
Vector3 
Matrix33::operator[] ( int col ) const {
  Vector3 v;
  for ( int i = 0;i < 3;++i )
    v( i ) = data_[ col + i * 3 ];
  return v;
}

//inline 
Matrix33& 
Matrix33::operator+=( const Matrix33& m ) {
  for ( int i = 0; i < 3*3; i++ )
    data_[ i ] += m.data_[ i ];
  return *this;
}

//inline 
Matrix33& 
Matrix33::operator-=( const Matrix33& m ) {
  for ( int i = 0; i < 3*3; i++ )
    data_[ i ] -= m.data_[ i ];
  return *this;
}

//inline 
Matrix33 
Matrix33::operator*( const double& x ) const {
  Matrix33 temp;
  for ( int i = 0; i < 3*3; i++ )
    temp.data_[ i ] = data_[ i ] * x;
  return temp;
}

//inline 
Matrix33 
Matrix33::operator/( const double& x ) const {
  Matrix33 temp;
  double ix = 1.0 / x;
  for ( int i = 0; i < 3*3; i++ )
    temp.data_[ i ] = data_[ i ] * ix;
  return temp;
}

//inline 
void 
Matrix33::set_diag( Vector3 diag ) {
  *this = 0.0;
  for ( unsigned i = 0; i < 3; ++i ) {
    ( *this ) ( i, i ) = diag( i );
  }
}

void 
Matrix33::Eigen( std::vector<real>& Eval, std::vector<Vector3>& Evec ) const {
  //  std::cout << " in Eigen()..." << std::endl;
  gsl_matrix_view m = gsl_matrix_view_array( data_, 3, 3 );
  gsl_vector *eval = gsl_vector_alloc ( 3 );
  gsl_matrix *evec = gsl_matrix_alloc ( 3, 3 );
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc ( 3 );

  gsl_eigen_symmv ( &m.matrix, eval, evec, w );
  gsl_eigen_symmv_free ( w );
  //gsl_eigen_symmv_sort ( eval, evec, GSL_EIGEN_SORT_ABS_ASC );
  gsl_eigen_symmv_sort ( eval, evec, GSL_EIGEN_SORT_VAL_ASC );

  for ( unsigned i = 0; i < 3; ++i ) {
    //double eval_i = gsl_vector_get (eval, i);
    gsl_vector_view evec_i = gsl_matrix_column (evec, i);
    Eval.at(i) = gsl_vector_get ( eval, i );
		
    for ( unsigned j = 0; j < 3; ++j ){
      (Evec.at(i))(j) = gsl_vector_get( &evec_i.vector, j );
    }
  }
}

double 
Matrix33::det( void ) const {
  std::vector<real> Eval;
  std::vector<Vector3> Evec;
  this->Eigen( Eval, Evec );
  double d = 1.0;
  for ( unsigned i = 0; i < 3; ++i )
    d *= Eval.at( i );
  return d;
}

Vector3 
Matrix33::get_invariants(std::vector<Vector3> Evec) const {
  std::vector<real> Eval;
  Vector3 inv;
  this->Eigen(Eval,Evec);
  inv(0)=Eval[0]+Eval[1]+Eval[2];
  inv(1)=Eval[0]*Eval[1]+Eval[1]*Eval[2]+Eval[0]*Eval[2];
  inv(2)=Eval[0]*Eval[1]*Eval[2];
  return inv;
}

void 
ROTATE(Matrix33 a,int i,int j,int k,int l,double s, double tau){
  double g=(a)(i,j);
  double h=(a)(k,l);
  (a)(i,j)=g-s*(h+g*tau);
  (a)(k,l)=h+s*(g-h*tau);
}

void 
ROTATE(std::vector<Vector3> a,int i,int j,int k,int l, double s, double tau){
  double g=a.at(i)(j);
  double h=a.at(k)(l);
  a.at(i)(j)=g-s*(h+g*tau);
  a.at(k)(l)=h+s*(g-h*tau);
}

/*
 * calculates Eigenvalues using Jacobi algorithm
 * use only if Eigen() is not applicable
 * cf. http://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm
 */
void 
Matrix33::Jacobi( std::vector<real>& Eval, std::vector<Vector3>& Evec, int &nrot ) {
  int j, iq, ip, i;
  double tresh, theta, tau, t, sm, s, h, g, c;
  double *b = new double[ 3 ],
    *z = new double[ 3 ];

  //... set output matrix to identity matrix ...
  for ( ip = 0; ip < 3; ip++ ) {
    for ( iq = 0; iq < 3; iq++ )
      Evec.at(ip)( iq ) = 0.0;
    Evec.at(ip)( ip ) = 1.0;
  }

  for ( ip = 0; ip < 3; ip++ ) {
    b[ ip ] = Eval[ip] = ( *this ) ( ip, ip );
    z[ ip ] = 0.0;
  }
  nrot = 0;
  for ( i = 0; i < 50;i++ ) {
    sm = 0.0;
    for ( ip = 0; ip < 2;ip++ ) {
      for ( iq = ip + 1;iq < 3; ++iq )
	sm += fabs( ( *this ) ( ip, iq ) );
    }
    if ( sm == 0.0 ) {
      delete[] b;
      delete[] z;
      return ;
    }
    if ( i < 4 )
      tresh = 0.2 * sm / ( 3 * 3 );
    else
      tresh = 0.0;
    for ( ip = 0; ip < 2; ip++ ) {
      for ( iq = ip + 1; iq < 3; iq++ ) {
	g = 100.0 * fabs( ( *this ) ( ip, iq ) );
	if ( i > 4 && ( double ) ( fabs( Eval[ip] ) + g ) == ( double ) fabs( Eval[ip] )
	     && ( double ) ( fabs( Eval[iq] ) + g ) == ( double ) fabs( Eval[iq] ) )
	  ( *this ) ( ip, iq ) = 0.0;
	else if ( fabs( ( *this ) ( ip, iq ) ) > tresh ) {
	  h = Eval[iq] - Eval[ip];
	  if ( ( double ) ( fabs( h ) + g ) == ( double ) fabs( h ) )
	    t = ( ( *this ) ( ip, iq ) ) / h;
	  else {
	    theta = 0.5 * h / ( ( *this ) ( ip, iq ) );
	    t = 1.0 / ( fabs( theta ) + sqrt( 1.0 + theta * theta ) );
	    if ( theta < 0.0 )
	      t = -t;
	  }
	  c = 1.0 / sqrt( 1 + t * t );
	  s = t * c;
	  tau = s / ( 1.0 + c );
	  h = t * ( ( *this ) ( ip, iq ) );
	  z[ ip ] -= h;
	  z[ iq ] += h;
	  Eval[ip] -= h;
	  Eval[iq] += h;
	  ( *this ) ( ip, iq ) = 0.0;
	  for ( j = 0; j < ip - 1; ++j ) {
	    ROTATE( *this, j, ip, j, iq, s, tau );
	  }
	  for ( j = ip + 1;j <= iq - 1; ++j ) {
	    ROTATE( *this, ip, j, j, iq, s, tau );
	  }
	  for ( j = iq + 1;j < 3; ++j ) {
	    ROTATE( *this, ip, j, iq, j, s, tau );
	  }
	  for ( j = 0;j < 3; ++j ) {
	    ROTATE( Evec, j, ip, j, iq, s, tau );
	  }
	  ++nrot;
	}
      }
    }
    for ( ip = 0; ip < 3; ++ip ) {
      b[ ip ] += z[ ip ];
      Eval[ip] = b[ ip ];
      z[ ip ] = 0.0;
    }
  }
  std::cerr << "Error : Too many iterations in routine Jacobi";
  delete[] b;
  delete[] z;
}

std::vector<Vector3>
Matrix33::toVecVector3() const{
  std::vector<Vector3> tvv;
  tvv.push_back(this->operator()(0));
  tvv.push_back(this->operator()(1));
  tvv.push_back(this->operator()(2));

  return tvv;
}
