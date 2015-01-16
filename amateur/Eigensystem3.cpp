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

/*
 * \file Eigensystem3.cpp
 * \brief implementations of Eigenvalue problem routine
 */

#include <cmath>
#include <vector>
#include <sstream>
#include <iostream>
#include "Global.h"
#include "Eigensystem3.h"
#include "Vector3.h"

Eigensystem3::Eigensystem3(){
  eigen_ = false;
}

Eigensystem3::Eigensystem3(Matrix33 M){
  M_ = M;
  eigen_ = true;
  for(int i=0; i<3; ++i){
    Eval_.push_back(0.0);
    Evec_.push_back(Vector3(0.0,0.0,0.0));
  }
}

Eigensystem3::Eigensystem3(const Eigensystem3& E){
  *this = E;
}

Eigensystem3::Eigensystem3(std::vector<real> Eval, std::vector<Vector3> Evec){
  eigen_ = true;
  Eval_ = Eval;
  Evec_ = Evec;
}

void
Eigensystem3::calcEV(){
  //  std::cout << " in calcEV()..." << std::endl;
  if(eigen_){
    M_.Eigen(Eval_, Evec_);
  } else {
    std::cerr << "Eigensystem3 not able to compute Eval and Evec, ";
    std::cerr << "supply a matrix first." <<std::endl;
    abort();
  }
}

std::vector<real>
Eigensystem3::getEval() const{
  //  std::cout << " in getEval()..." << std::endl;
  return Eval_;
}

real
Eigensystem3::getEval(int index) const{
  //  std::cout << " in getEval(int)..." << std::endl;
  return Eval_.at(index);
}

void
Eigensystem3::setEval(std::vector<real> Eval) {
  Eval_ = Eval;
}

std::vector<Vector3>
Eigensystem3::getEvec() const{
  //  std::cout << " in getEvec()..." << std::endl;
  return Evec_;
}

Vector3
Eigensystem3::getEvec(int index) const{
  //  std::cout << " in getEvec(int)..." << std::endl;
  return Evec_.at(index);
}

void
Eigensystem3::setEvec(std::vector<Vector3> Evec){
  Evec_ = Evec;
}

void
Eigensystem3::normalize(real factor){
  //  std::cout << " in normalize(real)..." << std::endl;
  if(eigen_){
    real la = Eval_.at(0);
    real lb = Eval_.at(1);
    real lc = Eval_.at(2);
    Eval_.at(0) = factor * std::sqrt( - la + lb + lc );
    Eval_.at(1) = factor * std::sqrt( + la - lb + lc );
    Eval_.at(2) = factor * std::sqrt( + la + lb - lc );
  } else {
    std::cerr << "Attempting to normalize ";
    std::cerr << "before Eigenvalues are determined." << std::endl;
    abort();
  }
}

std::string
Eigensystem3::toString(){
  if(!eigen_){
    std::cerr << "no Eigenvalues determined yet";
    abort();
  }
  std::ostringstream out;
  for(int i=0; i<3; ++i){
    if(std::isnan(Eval_.at(i))){
      out << 0 << " ";
    }else{
      out << Eval_.at(i) << " ";
    }
    out << Evec_.at(i).toString() << " ";
  }
  return out.str();
}
