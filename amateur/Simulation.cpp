/*
 * \file Simulation.cpp
 * characteristics for the simulation
 */

#include <cmath>
#include <cassert>

#include "Global.h"
#include "Cosmology.h"
#include "Simulation.h"

Simulation::Simulation(){
  //nothing to do

}

Simulation::Simulation(Simulation & sim){
  *this = sim;
}

Simulation::~Simulation(){
  //nothing to do as well
}

void 
Simulation::setCosmology(Cosmology c){
  c_ = c;
  setRhobox(); //which is a constant
}

void
Simulation::setTarkinOmega(){
  c_.setTarkinOmega();
}

void
Simulation::setIsnap(int isnap){
  isnap_ = isnap;
}

void 
Simulation::setAsnap(real asnap){
  asnap_ = asnap;
  zsnap_ = 1.0/asnap-1.0;
}

void 
Simulation::setZsnap(real zsnap){
  zsnap_ = zsnap;
  asnap_ = 1.0/(1.0+zsnap);
}

void
Simulation::setNhalo(int nhalo){
  nhalo_ = nhalo;
}

void
Simulation::setMass(real mass){
  mass_ = mass;
  calcRhobox();
}

void 
Simulation::setBoxlength(real boxlength){
  boxlength_ = boxlength;
  calcRhobox();
}

void
Simulation::calcRhobox(){
  rhobox_ = mass_ / std::pow(boxlength_,3);
}

//rhobox gives comoving density, independent of redshift of snapshot!
void
Simulation::setRhobox(){
  //assert(mass_ > 0.0);
  //assert(boxlength_ > 0.0);
  //assert(asnap_ > 0.0);
  //rhobox_ = mass_/(std::pow(boxlength_,3));
  //rhobox_ = 1.8791E-29 * Omega_M; // in gm/cm^3
  real h = c_.getH0()/100.0;
  rhobox_ = 2.7757819e-8 * c_.getOmegaM() / (h * h);	//in 10^10 msun/(h kpc)^3, in the units of tarkin
  //#warning overdensity criterion neglected, using all particles from AHF
  //rhobox_ = 0.0;	// gives loops over all particles, independent of overdensity
}






int
Simulation::getIsnap(){
  return isnap_;
}

real
Simulation::getAsnap(){
  return asnap_;
}

real 
Simulation::getZsnap(){
  return zsnap_;
}

int
Simulation::getNhalo(){
  return nhalo_;
}

real
Simulation::getMass(){
  return mass_;
}

real
Simulation::getBoxlength(){
  return boxlength_;
}

real
Simulation::getRhobox(){
  return rhobox_;
}

real
Simulation::getH0(){
  return c_.getH0();
}

real
Simulation::getG(){
  return c_.getG();
}

real
Simulation::getOmegaL(){
  return c_.getOmegaL();
}

real
Simulation::getOmegaM(){
  return c_.getOmegaM();
}

real
Simulation::getOmegaR(){
  return c_.getOmegaR();
}

real
Simulation::getOmegaK(){
  return c_.getOmegaK();
}
