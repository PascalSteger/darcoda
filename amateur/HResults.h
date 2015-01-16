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
 * \file HResults.h
 * \brief storage class for amateur::Halo output
 */

#ifndef HRESULTS_H
#define HRESULTS_H

#include "Global.h"
#include "Vector3.h"
#include "Eigensystem3.h"

const int HTotal = 5;

class HResults {
 private:
  int ntot_;

  real Mtot_;
  real Rtot_;
  real Mvir_;
  real Rvir_;
  Vector3 xcm_;
  Vector3 vcm_;
  Vector3 xmb_;
  Vector3 vmb_;

  int npart_[HTotal];
  real mpart_[HTotal];
  real lambda_[HTotal];
  real sigma_[HTotal];
  Vector3 J_[HTotal];
  Eigensystem3 E_[HTotal];
  
 public:
  HResults();
  ~HResults();

  void setNtot(int ntot);
  void setMtot(real Mtot);
  void setRtot(real Rtot);
  void setMvir(real Mvir);
  void setRvir(real Rvir);
  void setXCM(Vector3 xcm);
  void setVCM(Vector3 vcm);
  void setXMB(Vector3 xmb);
  void setVMB(Vector3 vmb);

  void setNpart(int npart, int type);
  void setMpart(real mpart, int type);
  void setLambda(real lambda, int type);
  void setSigma(real sigma, int type);
  void setJ(Vector3 J, int type);
  void setE(Eigensystem3 E, int type);

  std::string toString();
};

#endif
