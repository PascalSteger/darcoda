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
 * \file HResults.cpp
 * \brief all routines for HResults.h
 * contains all different routines for storing and writing results from Halo
 */

#include <sstream>
#include "Global.h"
#include "HResults.h"
#include "Vector3.h"
#include "Eigensystem3.h"

HResults::HResults(){

}

HResults::~HResults(){

}

void
HResults::setNtot(int ntot){
  ntot_ = ntot;
}

void
HResults::setMtot(real Mtot){
  Mtot_ = Mtot;
}

void
HResults::setRtot(real Rtot){
  Rtot_ = Rtot;
}

void
HResults::setMvir(real Mvir){
  Mvir_ = Mvir;
}

void
HResults::setRvir(real Rvir){
  Rvir_ = Rvir;
}

void
HResults::setXCM(Vector3 xcm){
  xcm_ = xcm;
}

void
HResults::setVCM(Vector3 vcm){
  vcm_ = vcm;
}

void
HResults::setXMB(Vector3 xmb){
  xmb_ = xmb;
}

void
HResults::setVMB(Vector3 vmb){
  vmb_ = vmb;
}

void
HResults::setNpart(int npart, int type){
  npart_[type] = npart;
}

void
HResults::setMpart(real mpart, int type){
  mpart_[type] = mpart;
}

void
HResults::setJ(Vector3 J, int type){
  J_[type] = J;
}

void
HResults::setLambda(real lambda, int type){
  lambda_[type] = lambda;
}

void
HResults::setSigma(real sigma, int type){
  sigma_[type] = sigma;
}

void
HResults::setE(Eigensystem3 E, int type){
  E_[type] = E;
}

/**
 * \brief output a HResult as string
 * concatenate first all global values, then all halo type specific properties
 */
std::string
HResults::toString(){
  std::ostringstream out;
  out << ntot_ << " ";
  out << Mtot_ << " ";
  out << Rtot_ << " ";
  out << Mvir_ << " ";
  out << Rvir_ << " ";

  out << xcm_.toString() << " ";
  out << vcm_.toString() << " ";
  out << xmb_.toString() << " ";
  out << vmb_.toString() << " ";

  for(int i=0;i<HTotal;++i){
    out << " ";
    out << npart_[i] << " ";
    out << mpart_[i] << " ";
    out << lambda_[i] << " ";
    out << sigma_[i] << " ";
    out << J_[i].toString() << " ";
    out << E_[i].toString() << " ";
  }
  return out.str();
}
