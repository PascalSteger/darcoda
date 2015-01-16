/***************************************************************************
 *   Copyright (C) 2010 by Pascal Stephan Philipp Steger                   *
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
 * \file Matrix33.h
 * \brief All publicly known functions of Matrix33.cpp
 * presents a 3x3 Matrix
 */

#ifndef MATRIX33_H
#define MATRIX33_H

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
#include "Vector3.h"

class Matrix33 {
 protected:
  double *data_;
 public:
  Matrix33();
  Matrix33( const Matrix33& m );
  ~Matrix33();

  double& operator() ( int row, int col );
  double operator() ( int row, int col ) const;
  Vector3 operator() ( int col ) const;
  Matrix33& operator= ( Matrix33 const &M );
  Matrix33& operator= ( const double& x );
  Vector3 operator[] ( int col ) const;
  Matrix33& operator+=( const Matrix33& m );
  Matrix33& operator-=( const Matrix33& m );
  Matrix33 operator*( const double& x ) const;
  Matrix33 operator/( const double& x ) const;
  void set_diag( Vector3 diag );
  void Jacobi( std::vector<real>& Eval, std::vector<Vector3>& Evec, int& nrot );
  void Eigen( std::vector<real>& Eval, std::vector<Vector3>& Evec ) const;
  double det( void ) const;
  Vector3 get_invariants( std::vector<Vector3> Evec ) const;

  std::vector<Vector3> toVecVector3() const;
};

#endif
