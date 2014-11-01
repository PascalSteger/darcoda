//HaloProperties (c) 2008 by Oliver Hahn, hahn@phys.ethz.ch

/**
 * \file HaloProperties.h
 * \brief Class definition with all publicly known functions from HaloProperties
 */

#ifndef HALOPROPERTIES_H
#define HALOPROPERTIES_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <ctime>

#include "DataStruct.h"
#include "OctTree.h"

#define OCTTREE
#define FORCE_DIRECT_TREE_SPLIT 4000

// typedef float real;
//typedef double gsl_real;

inline real 
fbmod( real x, real m ) {
  return x - floor( x / m + 0.5 ) * m;
}

template< typename real_t >
inline real_t 
SPHK3D( const real_t r, const real_t h ){
  const real_t xi(r/h);
  const real_t norm(8.0/M_PI/(h*h*h));
	
  if( xi < 0.5 )
    return norm*(6.0*xi*xi*(xi-1.0)+1.0);
  else if( xi <= 1.0 ){
    real_t t(1.0-xi);
    return norm*2.0*t*t*t;
  }
  return 0.0;
}

class HaloProperties {
 protected:
  real m_xcm[ 3 ], m_xmb[ 3 ], m_vcm[ 3 ], m_vmb[ 3 ];
  real m_rmbcm;
  real *m_pPos, *m_pVel, *m_pMass;
  real m_mtot;
  unsigned m_nPart;

  real m_alpha;
  unsigned m_imostbound;

  struct Cosm {
    real OmegaM0, OmegaL0, H0, G, asnap;
  };

  Cosm m_Cosm;

 public:

  class Distance {
    real m_r;
    unsigned m_i;
  public:
    Distance( real r, unsigned i ) : m_r( r ), m_i( i ) { }
    Distance( const Distance& d ) : m_r( d.m_r ), m_i( d.m_i ) { }
    bool operator<( const Distance& d2 ) const {
      return m_r < d2.m_r;
    }
    real getr( void ) const {
      return m_r;
    }
    unsigned geti( void ) const {
      return m_i;
    }
  };


  HaloProperties( real* pPos, real* pVel, real* pMass, unsigned nPart, real zsnap );

  void SetCosmology( real OmegaM0, real OmegaL0, real H0, real G );

  real ComputeTotalMass( void );
  real ComputeTotalMassExc( std::vector<bool> excluded );

  void ComputeCM( real Boxlength, Vector3& xcm );
  void ComputeCMExc( real Boxlength, Vector3& xcm, std::vector<bool> excluded );
  void ComputeMB( real Boxlength, Vector3& xmb );
  void ComputeMBExc( real Boxlength, Vector3& xmb, std::vector<bool> excluded );

  // assumes extension smaller than 0.25*boxlength
  void CenterOnCM( real Boxlength );
  void CenterOnCMExc( real Boxlength, std::vector<bool> excluded );

  real CenterOnMB( real Boxlength );
  real CenterOnMB( real Boxlength, real &rmbcm, real &alpha );
		
  void CenterOnCoord( real *xc, real Boxlength );
		
  void AverageVel( const unsigned n, const real eps=0.0 );

  void ComputeCenterVel( Vector3& vcm );
  void ComputeCenterVelExc( Vector3& vcm, std::vector<bool> excluded );
  // subtracts the center of mass velocity from all velocities

  void CenterVel( void );
  void CenterVelOn( real *vc );
  void CenterVelExc( std::vector<bool> excluded );

  void ComputeAngularMomentum( real& Jtot );
  real ComputeAngularMomentum( Vector3& J ) const;
  
  void ComputeVelocityDispersion( real& Sigma );

  void GetSortedDistances( std::vector<Distance>& dvec ) const;

  // computed relative to coordinate origin
  void ComputeMajorAxes( Vector3& Eigenvalues, std::vector<Vector3>& Eigenvectors ) const;
  void ComputeMajorAxesReduced( Vector3& Eigenvalues, std::vector<Vector3>& Eigenvectors ) const;

  void ComputeShapeParameters( real& T, real& S ) const;

  // PHYSICAL!!!!
  void ComputeVirialProperties200( real rhobox, real& Mvir, real& Rvir ) const;
  void ComputeVirialProperties200( real rhobox, real& Mvir, real& Rvir, real &Vcmax ) const;
  void ComputeVirialProperties200( real rhobox, real& Mvir, real& Rvir, std::vector<Distance> &dvec ) const;
  void ComputeVirialProperties200( real rhobox, real& Mvir, real& Rvir, real& Vcmax, std::vector<Distance> &dvec ) const;

  void ComputeVirialProperties( real rhobox, real& Mvir, real& Rvir ) const ;
  void ComputeVirialProperties( real rhobox, real& Mvir, real& Rvir, real& Vcmax ) const;
  void ComputeVirialProperties( real rhobox, real& Mvir, real& Rvir, std::vector<Distance> &dvec ) const;
  void ComputeVirialProperties( real rhobox, real& Mvir, real& Rvir, real& Vcmax, std::vector<Distance> &dvec ) const;

  real ComputePoormanConcentration( real rhobox, real &R20, real &R80 );
  real ComputePoormanConcentration( real rhobox );

  void ComputeLambdaPrime( real &lambda, real Mvir, real Rvir );

  real ComputeKineticEnergy( std::vector<real>& kin ) const;
  real ComputeKineticEnergyNoHubble( std::vector<real>& kin ) const;
  void ComputeKineticEnergyNoHubbleExc( std::vector<real>& kin, std::vector<bool> excluded ) const;
  // 		real ComputeKineticEnergyNoHubbleExc( std::vector<real>& kin, std::vector<bool> excluded ) const;

  real ComputePotentialDirect( real plumsoft, std::vector<real>& phi ) const;
  void ComputePotentialDirectExc( real plumsoft, std::vector<real>& phi, std::vector<bool> excluded ) const;
  real ComputePotentialTree( real theta, real plumsoft, std::vector<real>& phi ) const;
  void ComputePotentialTreeExc( real theta, real plumsoft, std::vector<real>& phi, std::vector<bool> excluded ) const;

  std::vector<bool> CleanUnbound( real boxsize );
  void CleanUnboundExc( real boxsize, std::vector<bool>& excluded );
  void CleanUnboundSimple( real boxsize, std::vector<bool>& excluded );

  int ComputeDensityProfile( real rhobox, unsigned n, std::vector<real>& r,
			     std::vector<real>& rho, std::vector<real>& sigma_rho ) const;

  int FitNFWProfile2( real rhobox, unsigned n, real& fdeltac, real& fedeltac,
		      real& frs, real& fers ) const;
  int FitNFWProfile1( real rhobox, unsigned n, real& frs, real& fers ) const;
  int FitNFWProfile200( real rhobox, unsigned n, real& frs, real& fers ) const;
};

#endif
