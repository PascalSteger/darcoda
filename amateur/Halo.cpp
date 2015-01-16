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


/****************************************************************************
 * \file Halo.cpp
 * \brief all methods of Halo.h
 ****************************************************************************/
#include <algorithm>
#include <vector>

#include "Global.h"
#include "Distance.h"
#include "Vector3.h"
#include "Matrix33.h"
#include "Halo.h"
#include "Simulation.h"

Halo::Halo(){
  //nothing to do
}

Halo::~Halo(){
  // nothing to do as well
}

Halo::Halo(Simulation sim, real* pos, real* vel, real* mass, int* type, unsigned npart){
  sim_ = sim;
  pos_ = pos;
  vel_ = vel;
  mass_ = mass;
  type_ = type;
  npart_ = npart;
  
  calcMtot();  // independent on centering
  centerCM();  // first centering on reference point, no bndry check anymore
  centerVel(); // used for Ekin in centerMB
  centerMB();  // most bound gives one halo in a binary as the main ref. pt
  centerVel(); // necessary for J,Lambda'
  calcDvec();  // distances from CMB
  calcVir(sim_.getRhobox()); // rvir, mvir
}

//private

/*
 * sum all particle masses
 */
void
Halo::calcMtot(){
  std::cout << " in calcMtot..." << std::endl;
  // sum over all particles inside halo
  mtot_ = 0.0;
  for(unsigned i = 0; i<npart_; ++i){
    mtot_ += mass_[i];
  }
}

/*
 * determine virial radius and mass
 * from spherical overdensity criterion
 * stopping at radius of first particle outside
 * required density
 * TODO: better approximation using ellipsoid,
 * or using all particles in overdensity, median value
 */
void
Halo::calcVir(real rhobox){
 std::cout << " in calcVir..." << std::endl;
  // using VirialProperties
  mvir_ = 0.0;
rvir_ = 0.0;

  real OM = sim_.getOmegaM();
  real OL = sim_.getOmegaL();
real a = sim_.getAsnap();

  real pi43 = 4.0 * M_PI / 3.0;
  real fz =  OM / ( a*a*a ) / (OM/( std::pow(a,3))+ OL ) - 1.0;

  real rhocrit = rhobox / OM * std::pow(a,3);
  // mean density / critical density (e.g. Bryan, Norman 1998)
  real delta = ( 18.0 * PI * PI + 82.0 * fz - 39.0 * fz * fz );

  // for density delta following Bryan, Norman 1998
  real rhodelta = delta * rhocrit;
  // for simple 200* overdensity criterion
  //real rhodelta = 200.0 * rhobox;

  real Vc;
  real Vcmax = -1.0e+99;

  real minside = 0.0;
  real R = 0.0;
  real rhomean = 1.0e+99;

  std::cout << "dvec_ has size: " << dvec_.size() << std::endl;
  std::cout << "asnap is: " << sim_.getAsnap() << std::endl;
  for( unsigned ninside=0;
       rhomean > rhodelta && ninside < npart_;
       ++ninside) {
    //    if(ninside % 1000 == 10){
    //  std::cout << ninside << "/" << npart_ << ":" << rhodelta << " " << rhomean << std::endl;
    //std::cout << minside << std::endl;
    //}
    R = sim_.getAsnap() * dvec_[ninside].getr();
    minside += mass_[dvec_[ninside].geti()];

    rhomean = minside / ( pi43 * std::pow(R,3) );

    Vc = minside / R;
    if ( Vc > Vcmax )
      Vcmax = Vc;
  }

  Vcmax = Vcmax * sim_.getG();
  mvir_ = minside;
  rvir_ = R;
  //rvir_ = sim_.getAsnap()*dvec_[npart_-1].getr();

}

/*
 * given a coordinate xc
 * subtract it from all particle positions,
 * correcting for periodic boundaries
 */
void 
Halo::centerOnCoord( Vector3 xc ){
  std::cout << " in centerOnCoord..." << std::endl;
  //real boxlength = sim_.getBoxlength();
  for ( unsigned i = 0; i < npart_; ++i ) {
    for ( unsigned j = 0; j < 3; ++j ) {
      pos_[3*i+j] = pos_[3*i+j] - xc(j);//TODO: bndry conditions
    }
  }
}

/*
 * calculate center of mass (also for halos across boundary)
 * and subtract it from all positions of particles in halo
 */
void
Halo::centerCM(){
  std::cout << " in centerCM..." << std::endl;
  xcm_ = 0.0;

  Vector3 rfirst;
  rfirst(0) = pos_[0];
  rfirst(1) = pos_[1];
  rfirst(2) = pos_[2];

  real boxlength = sim_.getBoxlength();
  //real hBL = boxlength / 2.0;
  for ( unsigned i = 0; i < npart_; ++i ) {
    for ( unsigned j = 0; j < 3; ++j ) {
      real dr = pos_[3*i+j] - rfirst(j);
      //if ( dr >= hBL ) {
      //	dr -= boxlength;
      //} else if ( dr < -hBL ) {
      //dr += boxlength;
      //}
      xcm_(j) += dr * mass_[i];
    }
  }
  xcm_ = xcm_ / mtot_ + rfirst;

  centerOnCoord( xcm_ );
}

/**
 * \brief calculate position of most bound particle
 * subtract position of MB from all particle positions in this halo
 * \param void
 * \return void
 */
void
Halo::centerMB(){
  std::cout << " in centerMB..." << std::endl;
  std::vector<real> kin;
  std::vector<real> pot;
  real totkin = 0.0;
  real totpot = 0.0;

  totkin = calcEkinNoHubble( kin );
  if ( npart_ > FORCE_DIRECT_TREE_SPLIT ) {
    totpot = calcPotTree( 0.75, 0.001, pot );
  } else {
    totpot = calcPotDirect( 0.001, pot );
  }

  real Emin = 1.0e+99;
  unsigned imb = 0;

  for ( unsigned i = 0; i < npart_; ++i ) {
    real Eres = kin[i] + pot[i];
    if ( Eres < Emin ) {
      Emin = Eres;
      imb = i;
    }
  }
  xmb_ = Vector3(pos_[3*imb+0], pos_[3*imb+1], pos_[3*imb+2]);
  centerOnCoord( xmb_ );
  xmb_ += xcm_;
  vmb_ = Vector3(vel_[3*imb],vel_[3*imb+1],vel_[3*imb+2]);
}


/**
 * \brief calculate mean velocity
 * subtract mean velocity from all velocities
 * \param void
 * \return void
 */
void
Halo::centerVel(){
  std::cout << " in centerVel..." << std::endl;
  // set vcm to 0.0
  vcm_ = 0.0;
	
  // add up velocities with weight given by mass
  for ( unsigned i = 0; i < npart_; ++i ) {
    for ( unsigned j = 0; j < 3; ++j ) {
      vcm_(j) += vel_[3*i+j] * mass_[i];
    }
  }

  // normalize velocity and subtract it from all other velocities
  assert( mtot_ > 0.0 );
  vcm_ /= mtot_;
  for ( unsigned i = 0; i < npart_; ++i ){
    for(unsigned j=0; j<3; ++j){
      vel_[3*i+j] -= vcm_( j );
    }
  }
}

/**
 * \brief generate list of radial distances
 * radial distances calculated after a centerCM() or centerMB()
 * \param void
 * \return void
 */
void 
Halo::calcDvec() {
  std::cout << " in calcDvec..." << std::endl;
  //dvec_.clear();

  for ( unsigned i = 0; i < npart_; ++i) {
    real ds2 = 0.0;
    for ( int j = 0; j < 3 ;++j )
      ds2 += pos_[3*i+j] * pos_[3*i+j];
    dvec_.push_back( Distance( std::sqrt( ds2 ), i ) );
  }

  rtot_ = dvec_.back().getr();
  std::sort( dvec_.begin(), dvec_.end() );
  //  for(unsigned i = 0; i<100 && i<npart_; ++i)
  //  std::cout << dvec_[i].getr() << ",";
  std::cout << std::endl;
}

// public
real 
Halo::getMtot(){
  return mtot_;
}

real
Halo::getRtot(){
  return rtot_;
}

real 
Halo::getMvir(){
  return mvir_;
}

real 
Halo::getRvir(){
  return rvir_;
}

Vector3
Halo::getXCM(){
  return xcm_;
}

Vector3
Halo::getVCM(){
  return vcm_;
}

Vector3
Halo::getXMB(){
  return xmb_;
}

Vector3
Halo::getVMB(){
  return vmb_;
}

/**
 * \brief do we have to look at a particle in this halo type?
 * determines whether particle is to be considered in a given halo type
 * \param int ptype: particle type, (0,1,2,3,4,5) for (gas,halo,disk,bulge,stars,bndry)
 * \param int htype: halo type, (15,3,4,8,12) for (all,dm,gas,stars,visible)
 */
bool
Halo::consider(int ptype, int htype){
  switch(htype){
  case 15:
    return true;
  case 3:
    if(ptype==1 || ptype==2 || ptype==3 || ptype == 5){
      return true;
    }else {
      return false;
    }
  case 4:
    if(ptype == 0){
      return true;
    }else{
      return false;
    }
  case 8:
    if(ptype == 4){
      return true;
    }else{
      return false;
    }
  case 12:
    if(ptype==0 || ptype==4){
      return true;
    }else{
      return false;
    }
  default:
    return false;
  } 
}



/**
 * \brief calculate kinetic NRJ
 * kinetic NRJ using the classical definition
 * \param std::vector<real>& kin: vector of output kinetic NRJ
 * \return total kinetic NRJ
 */
real 
Halo::calcEkin( std::vector<real>& kin ) {
  std::cout << " in calcEkin..." << std::endl;
  real kintot = 0.0;
  kin.clear();

  Vector3 v;
  real asnap = sim_.getAsnap();
  real sqrta = std::sqrt( asnap );
  real Hza = sim_.getH0() * asnap 
    * std::sqrt( sim_.getOmegaM()/std::pow( asnap,3 ) + sim_.getOmegaL() );

  for ( unsigned i = 0; i < npart_; ++i ) {
    v = 0.0;
    real k = 0;
    for ( unsigned j = 0; j < 3; ++j ) {
      v(j) = sqrta * vel_[3*i+j] + Hza * pos_[3*i+j];
    }
    k = 0.5 * mass_[i] * v.norm2();
    
    kin.push_back( ( real ) k );
    kintot += k;
  }

  return kintot;
}

/**
 * \brief calculate kinetic NRJ
 * Ekin with velocity corrected for Hubble expansion
 * \param std::vector<real>& kin: vector of output kinetic NRJ
 * \return total kinetic NRJ
 */
real 
Halo::calcEkinNoHubble( std::vector<real>& kin ) {
  std::cout << " in calcEkinNoHubble..." << std::endl;
  real kintot = 0.0;
  kin.clear();

  Vector3 v;
  real sqrta = std::sqrt( sim_.getAsnap() );
  
  // need for memory
  kin.reserve(npart_);
  for ( unsigned i = 0; i < npart_; ++i ) {
    v = 0.0;
    real k = 0;
    for ( unsigned j = 0; j < 3; ++j ) {
      v(j) = sqrta * vel_[3*i+j];
    }
    
    k = 0.5 * mass_[i] * v.norm2();
    
    kin.push_back( real(k) );
    kintot += k;
  }
  return kintot;
}

/**
 * \brief calculate potential by adding up inter-particle potentials
 * slow method, using all interactions
 * \param real plumsoft: smoothing scale
 * \param std::vector<real>& phi: holding all 
 */
real
Halo::calcPotDirect( real plumsoft, std::vector<real>& phi ) {
  std::cout << " in calcPotDirect..." << std::endl;
  real phitot = 0.0;
  real asnap = sim_.getAsnap();

  plumsoft *= asnap;
  for ( unsigned i = 0; i < npart_; ++i ) {
    real phii = 0.0;
    for ( unsigned j = i + 1; j < npart_; ++j ) {
      real
	dx = asnap * ( pos_[3*i+0] - pos_[3*j+0] ),
	dy = asnap * ( pos_[3*i+1] - pos_[3*j+1] ),
	dz = asnap * ( pos_[3*i+2] - pos_[3*j+2] );

      real r2 = dx * dx + dy * dy + dz * dz;
      if ( j != i )
	phii += -sim_.getG()*mass_[i]*mass_[i]/( std::sqrt( r2 )+plumsoft );
    }

    phi.push_back( phii );
    phitot += phii;
  }
  return phitot;
}

/*
 * \brief calculate potential using a OctTree
 * with opening angle theta
 * \param real theta: opening angle; determining when new leave is added
 * \param real plumsoft:
 * \param std::vector<real>& phi: potential getting output
 * \returns total potential NRJ
 */
real
Halo::calcPotTree( real theta, real plumsoft, std::vector<real>& phi ) {
  std::cout << " in calcPotTree..." << std::endl;
  real phitot = 0.0;
  phi.clear();
  real asnap = sim_.getAsnap();
  plumsoft *= asnap;

  std::vector< PseudoBody > barray;

  real xmin = 1.0e+99; //real xmax = -1.0e+99;
  real ymin = 1.0e+99; //real ymax = -1.0e+99;
  real zmin = 1.0e+99; //real zmax = -1.0e+99;

  for ( unsigned i = 0; i < npart_; ++i ) {
    xmin = ( pos_[3*i+0] < xmin ) ? pos_[3*i+0] : xmin;
    //xmax = ( pos_[3*i+0] > xmax ) ? pos_[3*i+0] : xmax;
    ymin = ( pos_[3*i+1] < ymin ) ? pos_[3*i+1] : ymin;
    //ymax = ( pos_[3*i+1] > ymax ) ? pos_[3*i+1] : ymax;
    zmin = ( pos_[3*i+2] < zmin ) ? pos_[3*i+2] : zmin;
    //zmax = ( pos_[3*i+2] > zmax ) ? pos_[3*i+2] : zmax;
  }

  real x, y, z;
  for ( unsigned i = 0, i3 = 0; i < npart_; ++i, i3 = 3 * i ) {
    x = asnap * ( pos_[3*i+0] - xmin );
    y = asnap * ( pos_[3*i+1] - ymin );
    z = asnap * ( pos_[3*i+2] - zmin );

    barray.push_back( PseudoBody( mass_[i], x, y, z ) );
  }

  //real length = std::max( std::max(xmax-xmin,ymax-ymin),zmax-zmin );

  //#ifndef OCTTREE
  //	kDTreeNode *BHroot = new kDTreeNode( barray );
  //#else
  OctTreeNode *BHroot = new OctTreeNode( barray, plumsoft ); //, length );
  //#endif

  for ( unsigned i = 0; i < npart_; ++i ) {
    real phii = sim_.getG()*mass_[i] * BHroot->traverse( barray[i].getx(),
							 barray[i].gety(),
							 barray[i].getz(),
							 theta, plumsoft );
    phitot += phii;
    phi.push_back( real( phii ) );
  }

  delete BHroot;

  return 0.5*phitot;
}


/**
 * \brief return number of particles of a given halo type
 * \param htype: int
 * \return total number of htype particles
 */
int
Halo::calcNpart(int htype){
  int count=0;
  for(unsigned i=0; i<npart_; ++i){
    if(consider(type_[i],htype)){
    ++count;
    }
  }
  return count;
}

/**
 * \brief return mass of particles of a given halo type
 * \param htype: int
 * \return total mass of htype particles
 */
real
Halo::calcMpart(int htype){
  real ms=0;
  for(unsigned i=0; i<npart_; ++i){
    if(consider(type_[i],htype)){
    ms+=mass_[i];
    }
  }
  return ms;
}


/**
 * \brief calculate velocity dispersion
 * sigma^2 = mean(v^2) - mean(v)^2
 * \param htype: int
 * \return sigma: real
 */
real
Halo::calcSigma(int htype){
  //std::cout << " in calcSigma..." << std::endl;
  real sigma = 0.0;
  // either go on with norm of velocity vector
  // use method minimizing roundoff errors
  real a = 0.0, aold = 0.0;
  real q = 0.0, qold = 0.0;
  unsigned count = 0;
  for(unsigned i = 0; i<npart_; ++i){
    if(!consider(type_[i],htype))
      continue;
    ++count;

    real vx = vel_[3*i+0];
    real vy = vel_[3*i+1];
    real vz = vel_[3*i+2];
    real v2 = std::sqrt(vx*vx + vy*vy + vz*vz);

    //a_0 = 0
    //a_i = a_{i-1}+(x_i-a_{i-1})/i
    a = aold + (v2-aold)/(i+1);
    //q_0 = 0
    //q_i = q_{i-1}+(x_i-a_{i-1})(x_i-a_i)
    q = qold + (v2 - aold)*(v2 - a);
    
    aold = a;
    qold = q;
  }

  if(count == 0){
    sigma = 0.0;
  } else {
    sigma = std::sqrt(q/count);
  }

  // or step through all euclidean components
  // or through (r,phi,theta) components
  // or use other definition of sigma
  return sigma;
}

/**
 * \brief computes the simplified spin parameter lambda prime
 * computes the halo spin parameter lambda' (Bullock et al. 2001).
 * assumption: POS AND VEL CENTERED ON REFERENCE POINT (centerCM, centerMB)
 * \param htype: int
 * \return lambda'
 */
real
Halo::calcLambda(int htype){
  //std::cout << " in calcLambda..." << std::endl;
  real lambda = -1.0;

  Vector3 J = 0.0;
  Vector3 d = 0.0;
  Vector3 v = 0.0;
  real minside = 0.0;
  real asnap = sim_.getAsnap();
  real sqrta = std::sqrt( asnap );

  unsigned count = 0;//count valid particles

  for( unsigned ninside=0; ninside < npart_
	 && asnap * dvec_[ninside].getr() < rvir_;
       //#warning Not using virial radius for calcLambda
       ++ninside ){
    if(!consider(type_[ninside],htype)){
      continue;
    }
    ++count;
    for ( unsigned j = 0; j < 3; ++j ) {
      v(j) = sqrta * vel_[3*dvec_[ninside].geti()+j];
      d(j) = asnap * pos_[3*dvec_[ninside].geti()+j];
    }
    minside += mass_[dvec_[ninside].geti()];
    J += (d^v)*mass_[dvec_[ninside].geti()]; // d^v is cross product
  }
  
  if( count == 0 ){
    std::cerr << "No particles in halo" << std::endl;
    lambda = -1.0;
  } else {
    assert( mvir_ > 0.0 );
    assert( rvir_ > 0.0 );
    lambda = real( J.norm() 
		   / ( std::sqrt( 2.0*sim_.getG()*rvir_*mvir_ ) * mvir_ ) );
  }
  return lambda;
}


/**
 * \brief calculate angular momentum
 * Angular momentum without considering Hubble flow (ASS: too short scales)
 * \param htype: int
 * \return J: Vector3
 */
Vector3 
Halo::calcJ(int htype){
  //std::cout << " in calcJ..." << std::endl;
  Vector3 J = 0.0;

  Vector3 r,v;
  unsigned counter = 0;
  real asnap = sim_.getAsnap();
  // real Hza = sim_.getH0() * asnap * sqrt( sim_.getOmegaM() / ( asnap * asnap * asnap ) + sim_.getOmegaL() );
  for ( unsigned i = 0; i < npart_; ++i ) {
    if(!consider(type_[i],htype )){
      continue;
    }
    counter++;
    for ( unsigned j = 0; j < 3; ++j ) {
      v(j) = std::sqrt(asnap)*vel_[3*i+j];//+Hza*pos_[3*i+j]; //HubbleFlow
      r(j) = asnap * pos_[3*i+j];
    }
    J += (r^v)*mass_[i];
  }
  //Jtot = sqrt( J( 0 ) * J( 0 ) + J( 1 ) * J( 1 ) + J( 2 ) * J( 2 ) );
  return J;
}




/**
 * \brief calculate halo shape
 * determine inertia tensor of all particles of a given type,
 * search eigensystem
 * \param htype: int
 */
Eigensystem3
Halo::calcShape(int htype){
  //std::cout << " in calcShape..." << std::endl;
  int count = 0;
  Matrix33 I; I = 0.0;
  for ( unsigned i = 0; i < npart_; ++i) {
    if(!consider(type_[i],htype)){
      continue;
    }
    count++;
    I(0,0)+= mass_[i]*( pos_[3*i+1]*pos_[3*i+1] + pos_[3*i+2]*pos_[3*i+2]);
    I(0,1)+= mass_[i]*(-pos_[3*i+0]*pos_[3*i+1]);
    I(0,2)+= mass_[i]*(-pos_[3*i+0]*pos_[3*i+2]);

    I(1,0)+= mass_[i]*(-pos_[3*i+1]*pos_[3*i+0]);
    I(1,1)+= mass_[i]*( pos_[3*i+0]*pos_[3*i+0] + pos_[3*i+2]*pos_[3*i+2]);
    I(1,2)+= mass_[i]*(-pos_[3*i+1]*pos_[3*i+2]);

    I(2,0)+= mass_[i]*(-pos_[3*i+2]*pos_[3*i+0]);
    I(2,1)+= mass_[i]*(-pos_[3*i+2]*pos_[3*i+1]);
    I(2,2)+= mass_[i]*( pos_[3*i+0]*pos_[3*i+0] + pos_[3*i+1]*pos_[3*i+1]);
  }
  if(count==0){
    std::cerr << "No particles for calculation of Eigensystem" << std::endl;
    I.set_diag(Vector3(-1,-1,-1)/0.0);
  }
  Eigensystem3 E = Eigensystem3(I);
  E.calcEV();
  assert(mtot_>0.0);
  E.normalize(std::sqrt(5.0/(2.0*mtot_)));

  return E;
}


