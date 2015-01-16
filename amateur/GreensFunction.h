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
 * \file GreensFunction.h
 * \brief implementation of Green's function G
 */
#ifndef GREENSFUNCTION_H
#define GREENSFUNCTION_H

#include <cmath>


template < int nOrder >
class GreensFunction {
  float dx_;
 public:
  GreensFunction( float dx ): dx_(dx){};
  float apply( int k );
};

template <>
inline float GreensFunction<3>::apply( int k ) {
  float ss = sin( 0.5 * dx_ * ( float ) k );
  return 4.0*ss*ss;
}

template <>
inline float GreensFunction<5>::apply( int k ) {
  const static float a0 = 0.6666666666666666667;
  float ss = sin( 0.5 * dx_ * ( float ) k );
  return -a0*ss*ss*( cos( dx_ * ( float ) k ) - 7.0 );
}

template <>
inline float GreensFunction<7>::apply( int k ) {
  const static float a0 = -0.04444444444444444444;
  float ss = sin( 0.5 * dx_ * k );
  return a0*ss*ss*(23.0*cos( dx_*(float)k)
		   -2.0*cos(2.0*dx_*(float)k)-111.0);
}

template <>
inline float GreensFunction<9>::apply( int k ) {
  const static float a0 = -1.58730159e-3;
  float ss = sin( 0.5 * dx_ * ( float ) k );
  return a0*ss*ss*( 779.0 * cos( dx_ * ( float ) k ) 
		    - 110.0 * cos( 2.0 * dx_ * ( float ) k )
		    + 9.0 * cos( 3.0 * dx_ * ( float ) k ) - 3198.0 );
}

template <>
inline float GreensFunction<11>::apply( int k ) {
  const static float a0 = -3.17460317e-4;
  float ss = sin( 0.5 * dx_ * ( float ) k );
  return a0*ss*ss*( 4343.0 * cos( dx_ * ( float ) k ) 
		    - 774.0 * cos( 2.0 * dx_ * ( float ) k )
		    + 109.0 * cos( 3.0 * dx_ * ( float ) k ) 
		    - 8.0 * cos( 4.0 * dx_ * ( float ) k )
		    - 16270.0 );
}

#endif
