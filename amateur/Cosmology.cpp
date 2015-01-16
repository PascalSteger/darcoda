#include "Cosmology.h"
#include "Global.h"

Cosmology::Cosmology(){
  setConcordanceG();
  setConcordanceH0();
  setConcordanceOmega();
}

Cosmology::Cosmology(real H0){
  H0_=H0;
  setConcordanceG();
  setConcordanceOmega();
}

Cosmology::Cosmology(real OmegaM, real OmegaR, real OmegaL){
  setConcordanceH0();
  setConcordanceG();
  OmegaL_ = OmegaL;
  OmegaM_ = OmegaM;
  OmegaR_ = OmegaR;
  setOmegaK();
}

Cosmology::Cosmology(real H0, real OmegaM, real OmegaR, real OmegaL){
  H0_ = H0;
  setConcordanceG();
  OmegaL_ = OmegaL;
  OmegaM_ = OmegaM;
  OmegaR_ = OmegaR;
  setOmegaK();
}

void
Cosmology::setConcordanceH0(){
  H0_ = 70.8;
}

void
Cosmology::setConcordanceG(){
  G_ = 43011.7902;
}

void
Cosmology::setConcordanceOmega(){
  OmegaL_ = 0.726;
  OmegaM_ = 0.228 + 0.0456;
  OmegaR_ = 0.0004;
  setOmegaK();
}

void
Cosmology::setTarkinOmega(){
  OmegaL_ = 0.7;
  OmegaM_ = 0.3;
  OmegaR_ = 0.0;
  setOmegaK();
}

void
Cosmology::setOmegaK(){
  OmegaK_ = 1.0 - OmegaL_ - OmegaM_ - OmegaR_;
}

void
Cosmology::setH0(real H0){
  H0_ = H0;
}

real
Cosmology::getH0(){
  return H0_;
}

void
Cosmology::setG(real G){
  G_ = G;
}

real
Cosmology::getG(){
  return G_;
}

void
Cosmology::setOmegaL(real OmegaL){
  OmegaL_ = OmegaL;
  setOmegaK();
}

real
Cosmology::getOmegaL(){
  return OmegaL_;
}

void
Cosmology::setOmegaM(real OmegaM){
  OmegaM_ = OmegaM;
  setOmegaK();
}

real 
Cosmology::getOmegaM(){
  return OmegaM_;
}

void
Cosmology::setOmegaR(real OmegaR){
  OmegaR_ = OmegaR;
  setOmegaK();
}

real
Cosmology::getOmegaR(){
  return OmegaR_;
}

real
Cosmology::getOmegaK(){
  return OmegaK_;
}

