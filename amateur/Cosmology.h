/*
 * \file Cosmology.h
 * global properties of a LCDM cosmology
 */

#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include "Global.h"

class Cosmology{

 private:
  real H0_;
  real G_;
  real OmegaL_;
  real OmegaM_;
  real OmegaR_;
  real OmegaK_;

  void setConcordanceG();
  void setConcordanceH0();
  void setConcordanceOmega();
  void setOmegaK();

 public:
  Cosmology();
  Cosmology(real H0);
  Cosmology(real OmegaM, real OmegaR, real OmegaL);
  Cosmology(real H0, real OmegaM, real OmegaR, real OmegaL);
  void setH0(real H0);
  real getH0();
  void setG(real G);
  real getG();
  void setTarkinOmega();
  void setOmegaL(real OmegaL);
  real getOmegaL();
  void setOmegaM(real OmegaM);
  real getOmegaM();
  void setOmegaR(real OmegaR);
  real getOmegaR();
  // we adjust OmegaK by private setOmegaK() if new values for O_M,R,L appear
  real getOmegaK();
};

#endif
