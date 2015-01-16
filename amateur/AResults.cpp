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
 * \file AResults.cpp
 * \brief implementation of AHF results storage class
 */

#include <string>
#include <sstream>
#include "AResults.h"
#include "Vector3.h"
#include "Eigensystem3.h"

AResults::AResults(){

}

AResults::AResults(unsigned npart, unsigned nvpart, float Xc, float Yc, float Zc,
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
		   float r2, float lambdaE){
  npart_ = npart;  nvpart_ = nvpart;
  xcm_ = Vector3(Xc, Yc, Zc);
  vcm_ = Vector3(VXc, VYc, VZc);
  Mvir_ = Mvir;  Rvir_ = Rvir;
  Vmax_ = Vmax;  Rmax_ = Rmax;
  sigV_ = sigV;
  lambda_ = lambda;
  J_ = Vector3(Jx, Jy, Jz);
  std::vector<real> Eval;
  Eval.push_back(a); Eval.push_back(b); Eval.push_back(c);
  std::vector<Vector3> Evec;
  Evec.push_back(Vector3(Eax,Eay,Eaz));
  Evec.push_back(Vector3(Ebx,Eby,Ebz));
  Evec.push_back(Vector3(Ecx,Ecy,Ecz));
  E_ = Eigensystem3(Eval, Evec);

  ovdens_ = ovdens;
  Redge_ = Redge;
  nbins_ = nbins;
  Ekin_ = Ekin;  Epot_ = Epot;
  mbp_offset_ = mbp_offset;
  com_offset_ = com_offset;
  r2_ = r2;
  lambdaE_ = lambdaE;
}

AResults::~AResults(){

}

std::string
AResults::toString(){
  std::ostringstream os;
  os << npart_ << " ";
  os << nvpart_ << " ";
  os << Mvir_ << " ";
  os << Rvir_ << " ";
  os << xcm_.toString() << " ";
  os << vcm_.toString() << " ";
  os << lambda_ << " ";
  os << sigV_ << " ";
  os << J_.toString() << " ";
  os << E_.toString() << " ";
  return os.str();
}
