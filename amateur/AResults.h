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
 * \file AResults.h
 * \brief storage class for AHF output
 */

#ifndef ARESULTS_H
#define ARESULTS_H

#include <string>
#include "Global.h"
#include "Vector3.h"
#include "Eigensystem3.h"

class AResults{
 private:
  unsigned npart_;
  unsigned nvpart_;
  Vector3 xcm_;
  Vector3 vcm_;
  float Mvir_, Rvir_, Vmax_, Rmax_, sigV_;
  float lambda_;
  Vector3 J_;
  Eigensystem3 E_;
  float ovdens_, Redge_;
  int nbins_;
  float Ekin_, Epot_;
  float mbp_offset_, com_offset_;
  float r2_, lambdaE_;

 public:
  AResults();
  AResults(unsigned npart, unsigned nvpart, float Xc, float Yc, float Zc,
	   float VXc, float VYc, float VZc,
	   float Mvir, float Rvir,
	   float Vmax, float Rmax,
	   float sigV, float lambda, float Jx, float Jy, float Jz,
	   float a, float Eax, float Eay, float Eaz,
	   float b, float Ebx, float Eby, float Ebz,
	   float c, float Ecx, float Ecy, float Ecz,
	   float ovdens, float Redge,
	   int nbins,
	   float Ekin, float Epot,
	   float mbp_offset, float com_offset,
	   float r2, float lambdaE);
  ~AResults();

  std::string toString();

};

#endif
