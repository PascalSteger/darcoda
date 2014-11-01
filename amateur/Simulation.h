#ifndef SIMULATION_H
#define SIMULATION_H

#include "Global.h"
#include "Cosmology.h"

class Simulation{
 private:
  Cosmology c_;
  int isnap_;
  real asnap_;
  real zsnap_;
  int nhalo_;
  real mass_;
  real boxlength_;
  real rhobox_;

  void calcRhobox();
  void setRhobox();

 public:
  Simulation();
  Simulation(Simulation & sim);
  ~Simulation();
  void setCosmology(Cosmology c);
  void setIsnap(int isnap);
  void setAsnap(real asnap);
  void setZsnap(real zsnap);
  void setNhalo(int nhalo);
  void setMass(real mass);
  void setTarkinOmega();
  void setBoxlength(real boxlength);

  int getIsnap();
  real getAsnap();
  real getZsnap();
  int getNhalo();
  real getMass();
  real getBoxlength();
  real getRhobox();

  real getH0();
  real getG();
  real getOmegaL();
  real getOmegaM();
  real getOmegaR();
  real getOmegaK();
};

#endif
