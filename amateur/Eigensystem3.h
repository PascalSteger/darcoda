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
 * \file Eigensystem3.h
 * \ brief implements a framework for finding eigenvalues and eigenvectors * of a 3x3-matrix
 */

#ifndef EIGENSYSTEM3_H
#define EIGENSYSTEM3_H


#include <vector>
#include "Global.h"
#include "Vector3.h"
#include "Matrix33.h"

class Eigensystem3 {
 private:
  bool eigen_;
  std::vector<real> Eval_;
  std::vector<Vector3> Evec_;
  Matrix33 M_;
  

 public:
  Eigensystem3();
  Eigensystem3(Matrix33 M);
  Eigensystem3(const Eigensystem3& E);
  Eigensystem3(std::vector<real> Eval, std::vector<Vector3> Evec);

  bool getEigen();
  std::vector<real> getEval() const;
  void setEval(std::vector<real> Eval);
  real getEval(int index) const;
  std::vector<Vector3> getEvec() const;
  void setEvec(std::vector<Vector3> Evec);
  Vector3 getEvec(int index) const;

  void calcEV();
  void normalize(real factor);

  std::string toString( void );
};

#endif
