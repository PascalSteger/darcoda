/*****************************************************************************************************
 * class to handle halos
 *****************************************************************************************************/
/*
 * \file Halo.h
 * basic properties of a halo
 */

#ifndef HALO_H
#define HALO_H

#include <vector>

#include "Global.h"
#include "Vector3.h"
#include "Matrix33.h"
#include "Eigensystem3.h"
#include "Cosmology.h"
#include "Simulation.h"
#include "Distance.h"
#include "OctTree.h"

#define OCTTREE
#define FORCE_DIRECT_TREE_SPLIT 4000

class Halo {
 private:
  Simulation sim_;

  real* pos_;
  real* vel_;
  real* mass_;
  int*  type_;
  unsigned npart_;

  real mtot_;
  real rtot_;
  real mvir_;
  real rvir_;

  std::vector<Distance> dvec_;
  Vector3 xcm_;
  Vector3 vcm_;
  Vector3 xmb_;
  Vector3 vmb_;

 private:
  void calcMtot();
  void calcVir(real rhobox);
  void calcDvec();
  void centerOnCoord( Vector3 xc );
  void centerCM();
  void centerMB();
  void centerVel();

  real calcEkin( std::vector<real>& kin );
  real calcEkinNoHubble( std::vector<real>& kin );
  real calcPotDirect( real plumsoft, std::vector<real>& phi );
  real calcPotTree( real theta, real plumsoft, std::vector<real>& phi );

  bool consider(int ptype, int htype);

 public:
  Halo();
  ~Halo();
  Halo(Simulation sim, real* pos, real* vel, real* mass, int* type, unsigned npart);

  real getMtot();
  real getRtot();
  real getMvir();
  real getRvir();

  Vector3 getXCM();
  Vector3 getVCM();
  Vector3 getXMB();
  Vector3 getVMB();

  int calcNpart(int htype);
  real calcMpart(int htype);
  real calcSigma(int htype);
  real calcLambda(int htype);
  Vector3 calcJ(int htype);
  Eigensystem3 calcShape(int htype);
  
};

#endif
