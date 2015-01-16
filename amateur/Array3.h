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
 * \file Array3.h
 * \brief a threedimensional array (x,y,z)
 */

#ifndef ARRAY3_H
#define ARRAY3_H

#include "Global.h"

class Array3 {
 protected:
  float *data_;
  int nx_, ny_, nz_;
  bool RealContainer_;
 public:
  
  Array3( void );
  Array3( float *dataptr, int nx, int ny, int nz );
  Array3( int nx, int ny, int nz );
  ~Array3();
    
  //float& operator() ( int ix, int iy, int iz );
  //float operator() ( int ix, int iy, int iz );
  Array3* operator=(Array3* a);
  int getNx() {return nx_;};
  int getNy() {return ny_;};
  int getNz() {return nz_;};
  float* getData() {return data_;};
  bool getRealContainer() {return RealContainer_;};
  float getValue(int ix, int iy, int iz);
  void setValue(int ix, int iy, int iz, float val);
  void addValue(int ix, int iy, int iz, float val);
  int get_size( void );
  void recreate( int nx, int ny, int nz );
};

#endif
