// HaloProperties (c) 2008 by Oliver Hahn, hahn@phys.ethz.ch
//                (c) 2010    Pascal Steger, psteger@phys.ethz.ch

/**
 * \file HaloProperties.cpp
 * \brief Implementations of the HaloProperties class
 */

#include "HaloProperties.h"

#include <stdexcept>
#include "ANN/ANN.h"
#include "NFWFit.h"

HaloProperties::HaloProperties( real *pPos, real *pVel, real *pMass, unsigned nPart, real zsnap ) {
  m_pPos = pPos;
  m_pVel = pVel;
  m_nPart = nPart;
  m_pMass = pMass;
  m_Cosm.OmegaM0 = 0.25;
  m_Cosm.OmegaL0 = 0.75;
  m_Cosm.H0 = 0.1;
  m_Cosm.G = 43011.7902;
  m_Cosm.asnap = 1.0 / ( 1.0 + zsnap );

  m_imostbound = -1;

  for ( unsigned j = 0; j < 3; ++j ) {
    m_xcm[ j ] = 0.0;
    m_xmb[ j ] = 0.0;
    m_vcm[ j ] = 0.0;
    m_vmb[ j ] = 0.0;
  }
  m_mtot = 0.0;
  for( unsigned i=0; i < m_nPart; ++i)
    m_mtot += m_pMass[i];
}
/**
 * Set internal units and values for cosmological constants
 * @param OmegaM0 fraction of mass density and critical density
 * @param OmegaL0 fraction of false vacuum density and critical density
 * @param H0 present day Hubble constant
 * @param G gravitation constant, in internal units
 */
void 
HaloProperties::SetCosmology( real OmegaM0, real OmegaL0, real H0, real G = 43011.7902 ) {
  m_Cosm.OmegaM0 = OmegaM0;
  m_Cosm.OmegaL0 = OmegaL0;
  m_Cosm.H0 = H0;
  m_Cosm.G = 43011.7902;
}

real 
HaloProperties::ComputeTotalMass( void ) {
  return m_mtot;
}

real 
HaloProperties::ComputeTotalMassExc( std::vector<bool> excluded ) {
  real mtot = 0.0;
  for ( unsigned i = 0; i < m_nPart; ++i )
    if ( !excluded.at( i ) )
      mtot += m_pMass[ i ];
  return mtot;
}

void 
HaloProperties::ComputeCM( real Boxlength, Vector3& xcm ) {
  xcm = 0.0;

  real rfirst[ 3 ];
  for(unsigned j=0; j<3; ++j)
    rfirst[ j ] = m_pPos[ j ];

  real BL05 = Boxlength / 2.0;
  for ( unsigned i = 0, i3=0; i < m_nPart; ++i, i3=3*i ) {
    for ( unsigned j = 0; j < 3; ++j ) {
      real dr = m_pPos[ i3 + j ] - rfirst[ j ];
      if ( dr >= BL05 )
	dr -= Boxlength;
      else if ( dr < -BL05 )
	dr += Boxlength;
      xcm( j ) += dr * m_pMass[ i ];
    }
  }
  assert( m_mtot > 0.0 );
  for ( unsigned j = 0; j < 3; ++j )
    xcm( j ) = fmod( xcm( j ) / m_mtot + rfirst[ j ], Boxlength );
}


void 
HaloProperties::ComputeCMExc( real Boxlength, Vector3& xcm, std::vector<bool> excluded ) {
  xcm = 0.0;

  real rfirst[ 3 ];
  for(unsigned j=0; j<3; ++j)
    rfirst[ j ] = m_pPos[ j ];

  real BL05 = Boxlength / 2.0;
  real mtotexc = 0.0;
  for ( unsigned i = 0, i3=0; i < m_nPart; ++i, i3=3*i ) {
    if ( excluded.at( i ) )
      continue;
    for ( unsigned j = 0; j < 3; ++j ) {
      real dr = m_pPos[ i3 + j ] - rfirst[ j ];
      if ( dr >= BL05 )
	dr -= Boxlength;
      else if ( dr < -BL05 )
	dr += Boxlength;
      xcm( j ) += dr * m_pMass[ i ];
    }
    mtotexc += m_pMass[ i ];
  }
  for ( unsigned j = 0; j < 3; ++j )
    xcm( j ) = fmod(xcm( j ) / mtotexc + rfirst[ j ], Boxlength);
}

/*! centers on the CM obeying periodic boundary conditions. REQUIRES THAT THE EXTENSION
  OF THE OBJECT DO NOT EXCEED 0.25 OF THE BOX LENGTH.
  \param Boxlength the length of the simulation box
  \return nothing, centers the coordinates on the CM and stores the CM
*/
void 
HaloProperties::CenterOnCM( real Boxlength ) {
  m_xcm[ 0 ] = 0.0; m_xcm[ 1 ] = 0.0; m_xcm[ 2 ] = 0.0;

  real rfirst[ 3 ];
  rfirst[ 0 ] = m_pPos[ 0 ]; rfirst[ 1 ] = m_pPos[ 1 ]; rfirst[ 2 ] = m_pPos[ 2 ];

  real BL05 = Boxlength / 2.0;
  for ( unsigned i = 0; i < m_nPart; ++i ) {
    for ( unsigned j = 0; j < 3; ++j ) {
      real dr = m_pPos[ 3 * i + j ] - rfirst[ j ];
      if ( dr >= BL05 )
	dr -= Boxlength;
      else if ( dr < -BL05 )
	dr += Boxlength;
      m_xcm[ j ] += dr * m_pMass[ i ];
    }
  }
  for ( unsigned j = 0; j < 3; ++j )
    m_xcm[ j ] = fmod( m_xcm[ j ] / m_mtot + rfirst[ j ], Boxlength );

  for ( unsigned i = 0, i3 = 0; i < m_nPart; ++i, i3 = 3 * i ) {
    for ( unsigned j = 0; j < 3; ++j ) {
      m_pPos[ i3 + j ] -= m_xcm[ j ];

      if ( m_pPos[ i3 + j ] >= BL05 )
	m_pPos[ i3 + j ] -= Boxlength;
      else if ( m_pPos[ i3 + j ] < -BL05 )
	m_pPos[ i3 + j ] += Boxlength;
    }
  }
}
void 
HaloProperties::CenterOnCMExc( real Boxlength, std::vector<bool> excluded ) {
  real rfirst[ 3 ];
  for(unsigned j=0; j<3; ++j){
    m_xcm[ j ] = 0.0;
    rfirst[ j ] = m_pPos[ j ];
  }

  real BL05 = Boxlength / 2.0;
  real mtotexc = 0.0;
  for ( unsigned i = 0; i < m_nPart; ++i ) {
    if ( excluded.at( i ) )
      continue;
    for ( unsigned j = 0; j < 3; ++j ) {
      real dr = m_pPos[ 3 * i + j ] - rfirst[ j ];
      if ( dr >= BL05 )
	dr -= Boxlength;
      else if ( dr < -BL05 )
	dr += Boxlength;
      m_xcm[ j ] += dr * m_pMass[ i ];
    }
    mtotexc += m_pMass[ i ];
  }
  for ( unsigned j = 0; j < 3; ++j )
    m_xcm[ j ] = fmod( m_xcm[ j ] / mtotexc + rfirst[ j ], Boxlength );

  for ( unsigned i = 0, i3 = 0; i < m_nPart; ++i, i3 = 3 * i ) {
    if ( excluded.at( i ) )
      continue;
    for ( unsigned j = 0; j < 3; ++j ) {
      m_pPos[ i3 + j ] -= m_xcm[ j ];

      if ( m_pPos[ i3 + j ] >= BL05 )
	m_pPos[ i3 + j ] -= Boxlength;
      else if ( m_pPos[ i3 + j ] < -BL05 )
	m_pPos[ i3 + j ] += Boxlength;
    }
  }
}

void 
HaloProperties::ComputeMB( real Boxlength, Vector3& xmb ) {
  std::vector<real> kin;
  std::vector<real> pot;
  real totkin = 0.0;
  real totpot = 0.0;

  totkin = ComputeKineticEnergyNoHubble( kin );

  if ( m_nPart > FORCE_DIRECT_TREE_SPLIT )
    totpot = ComputePotentialTree( 0.75, 0.001, pot );
  else
    totpot = ComputePotentialDirect( 0.001, pot );

  real Emin = 1.0e+99; // big number, not less than maximum
  unsigned imostbound = 0;

  for ( unsigned i = 0; i < m_nPart; ++i ) {
    if ( kin[ i ] + pot[ i ] < Emin ) {
      Emin = kin[ i ] + pot[ i ];
      imostbound = i;
    }
  }

  for(unsigned j=0; j<3; xmb( j ) = m_pPos[ 3 * imostbound + j ], ++j);
}

void 
HaloProperties::ComputeMBExc(real Boxlength, Vector3& xmb, std::vector<bool> excluded ) {

  std::vector<real> kin, pot;
  real totkin = 0.0, totpot = 0.0;

  totkin = ComputeKineticEnergyNoHubble( kin );

  if ( m_nPart > FORCE_DIRECT_TREE_SPLIT )
    totpot = ComputePotentialTree( 0.75, 0.001, pot );
  else
    totpot = ComputePotentialDirect( 0.001, pot );

  real Emin = 1.0e+99; // big number, not less than maximum
  unsigned imostbound = 0;

  for ( unsigned i = 0; i < m_nPart; ++i ) {
    if( excluded.at(i) )
      continue;
    if ( kin[ i ] + pot[ i ] < Emin ) {
      Emin = kin[ i ] + pot[ i ];
      imostbound = i;
    }
  }

  for(unsigned j=0; j<3; xmb( j ) = m_pPos[ 3 * imostbound + j ], ++j);
}

/**
 * centers all coordinates on the MB obeying periodic boundary conditions ... centers on the MB obeying periodic boundary conditions. REQUIRES THAT THE EXTENSION OF THE OBJECT DOES NOT EXCEED 0.25 * BOX LENGTH.
 * @param Boxlength 
 * @return distance of mb from cm, centers the coordinates on the MB and stores the MB
 */
real 
HaloProperties::CenterOnMB( real Boxlength ) {
  real r, a;
  return CenterOnMB( Boxlength, r, a );
}

/**
 * Center on Most Bound particle. CenterOnCM( Boxlength ); CenterVel(); is needed before invoking CenterOnMB
 * @param Boxlength 
 * @param rmbcm 
 * @param alpha 
 * @return 
 */
real 
HaloProperties::CenterOnMB( real Boxlength, real &rmbcm, real &alpha ) {
  std::vector<real> kin, pot;
  real totkin = 0.0, totpot = 0.0;

  totkin = ComputeKineticEnergyNoHubble( kin );

  if ( m_nPart > FORCE_DIRECT_TREE_SPLIT )
    totpot = ComputePotentialTree( 0.75, 0.001, pot );
  else
    totpot = ComputePotentialDirect( 0.001, pot );

  real Emin = 1.0e+99;
  unsigned imostbound = 0;

  for ( unsigned i = 0; i < m_nPart; ++i ) {
    if ( kin[ i ] + pot[ i ] < Emin ) {
      Emin = kin[ i ] + pot[ i ];
      imostbound = i;
    }
  }

  for(unsigned j=0; j<3; m_xmb[ j ] = m_pPos[ 3 * imostbound + j ], ++j);

  m_rmbcm = sqrt( m_xmb[ 0 ] * m_xmb[ 0 ]
		  + m_xmb[ 1 ] * m_xmb[ 1 ]
		  + m_xmb[ 2 ] * m_xmb[ 2 ] );

  m_alpha = 2.0 * totkin / totpot + 1.0;

  CenterOnCoord( m_xmb, Boxlength );

  m_imostbound = imostbound;
  alpha = m_alpha;
  rmbcm = m_rmbcm;
  return m_rmbcm;
}


/**
 * centers all coordinates on specified point
 * @param xc pointer to an array of three reals contain the point to be centered on
 * @param Boxlength size of simulation (cubic) box
 */
void 
HaloProperties::CenterOnCoord( real *xc, real Boxlength ) {
  for ( unsigned i = 0, i3 = 0; i < m_nPart; ++i, i3 = 3 * i )
    for ( unsigned j = 0; j < 3; ++j )
      m_pPos[ i3 + j ] = fmod( m_pPos[ i3 + j ] - xc[ j ], Boxlength );
}


/**
 * average the velocities of particles using ANN
 */
void 
HaloProperties::AverageVel( const unsigned n, const real eps ) {
  ANNpointArray		dataPts;
  ANNpoint 		queryPt;
  ANNidxArray 		nnIdx;
  ANNdistArray 		dists;
  ANNkd_tree* 		kdTree;
	
  queryPt = annAllocPt( 3 );
  dataPts = annAllocPts( m_nPart, 3 );
  nnIdx   = new ANNidx[n];
  dists   = new ANNdist[n];
	
  if( m_nPart <= n ){
    throw std::runtime_error( "in HaloProperties::AverageVel: number of points in averaging >= number of points available" );
  }
	
  std::vector<real> SPHrho(m_nPart,0.0);
  std::vector<real> SPHvel(3*m_nPart,0.0);
	
  for( unsigned i=0,i3=0; i<m_nPart; ++i,i3+=3 )
    for( unsigned j=0; j<3; ++j )
      (dataPts[i])[j] = m_pPos[i3+j];
	
  kdTree = new ANNkd_tree( dataPts, m_nPart, 3 );
	
  for( unsigned i=0,i3=0; i<m_nPart; ++i,i3+=3 ){
		
    for( unsigned j=0; j<3; ++j )
      queryPt[j] = m_pPos[i3+j];
		
    kdTree->annkSearch( queryPt, n, nnIdx, dists, eps );
		
    real h = sqrt(dists[n-1]);

    for( unsigned k=0; k<n; ++k ){
      dists[k] = sqrt(dists[k]);
      SPHrho[i] += m_pMass[nnIdx[k]]*SPHK3D( dists[k], (double)h );
    }
  }
	
  for( unsigned i=0,i3=0; i<m_nPart; ++i,i3+=3 ){
		
    for( unsigned j=0; j<3; ++j )
      queryPt[j] = m_pPos[i3+j];
		
    kdTree->annkSearch( queryPt, n, nnIdx, dists, eps );

    real h = sqrt(dists[n-1]);	
		
    for( unsigned k=0; k<n; ++k ){
      dists[k] = sqrt(dists[k]);
      for( unsigned j=0; j<3; ++j )
	SPHvel[i3+j] += m_pMass[nnIdx[k]]/SPHrho[nnIdx[k]]*m_pVel[3*nnIdx[k]+j]*SPHK3D( dists[k], (double)h );
    }
  }
	
  delete[] nnIdx;
  delete[] dists;
  delete kdTree;
	
  for( unsigned i=0,i3=0; i<m_nPart; ++i,i3+=3 ){
    for( unsigned j=0; j<3; ++j )
      m_pVel[i3+j] = SPHvel[i3+j];
  }
}

void 
HaloProperties::ComputeCenterVel( Vector3& vcm ) {
  for(unsigned j=0; j<3; ++j)
    m_vcm[j]= 0.0;

  for ( unsigned i = 0, i3 = 0; i < m_nPart; ++i, i3 = 3 * i )
    for ( unsigned j = 0; j < 3; ++j )
      m_vcm[ j ] += m_pVel[ i3 + j ] * m_pMass[ i ];

  assert( m_mtot > 0.0 );
  for( unsigned j = 0; j < 3; ++j){
    m_vcm[ j ] /= m_mtot;
    vcm(j)=m_vcm[j];
  }
}

void 
HaloProperties::ComputeCenterVelExc( Vector3& vcm, std::vector<bool> excluded ) {
  real mtotexc = 0.0;
  vcm = 0.0;

  for ( unsigned i = 0, i3=0; i < m_nPart; ++i, i3=3*i ) {
    if ( excluded.at( i ) )
      continue;
    for ( unsigned j = 0; j < 3; ++j ) {
      vcm( j ) += m_pVel[ i3 + j ] * m_pMass[ i ];
    }
    mtotexc += m_pMass[ i ];
  }

  if( ! mtotexc > 0.0 ) {
    vcm = 0.0;
    return;
  }
  vcm /= mtotexc;
}

void 
HaloProperties::CenterVel( void ) {
  // set vcm to 0.0
  for(unsigned j=0; j<3; m_vcm[ j ] = 0.0, ++j);
	
  // add up velocities with weight given by mass
  for ( unsigned i = 0, i3 = 0; i < m_nPart; ++i, i3=3*i ) {
    for ( unsigned j = 0; j < 3; ++j ) {
      m_vcm[ j ] += m_pVel[ i3 + j ] * m_pMass[ i ];
    }
  }

  // normalize velocity and subtract it from all other velocities
  assert( m_mtot > 0.0 );
  for(unsigned j=0; j<3; ++j){
    m_vcm[ j ] /= m_mtot;
    for ( unsigned i = 0, i3=0; i < m_nPart; m_pVel[ i3 + j ] -= m_vcm[ j ], ++i, i3=3*i );
  }
  // set vcm to 0.0
  for( unsigned j=0; j<3; m_vcm[ j ] = 0.0, ++j );
}

void 
HaloProperties::CenterVelOn( real* vc ) {
  for( unsigned j=0; j<3; ++j){
    m_vcm[j]=0.0;
    for ( unsigned i = 0, i3 = 0; i < m_nPart; m_pVel[ i3 + j ] = m_pVel[ i3 + j ] - vc[ j ], ++i, i3 = 3 * i );
  }
}

void 
HaloProperties::CenterVelExc( std::vector<bool> excluded ) {
  real mtotexc = 0.0;
  for(unsigned j = 0; j < 3; ++j)
    m_vcm[ j ] = 0.0;

  for ( unsigned i = 0; i < m_nPart; ++i ) {
    if ( excluded.at( i ) )
      continue;
    for ( unsigned j = 0; j < 3; ++j ) {
      m_vcm[ j ] += m_pVel[ 3 * i + j ] * m_pMass[ i ];
    }
    mtotexc += m_pMass[ i ];
  }

  assert( mtotexc > 0.0 );
  for(unsigned j=0; j<3; ++j){
    m_vcm[ j ] /= mtotexc;
    for ( unsigned i = 0, i3=0; i < m_nPart; ++i, i3=3*i )
      m_pVel[ i3 + j ] -= m_vcm[ j ];
  }
}





/**
 * computes the angular momentum in internal units.
 * assumes that coordinates and velocities are centered on a reference point
 * @param Jtot total angular momentum
 * @return void
 */
void 
HaloProperties::ComputeAngularMomentum( real &Jtot ) {
  Jtot = 0.0;
  real r[ 3 ], v[ 3 ], J[ 3 ];
  //real Hza = m_Cosm.H0 * m_Cosm.asnap * sqrt( m_Cosm.OmegaM0 / ( m_Cosm.asnap * m_Cosm.asnap * m_Cosm.asnap )+ m_Cosm.OmegaL0 );
  real sqrta = sqrt( m_Cosm.asnap );
  J[ 0 ] = J[ 1 ] = J[ 2 ] = 0.0;
  for ( unsigned i = 0; i < m_nPart; ++i ) {
    for ( unsigned j = 0; j < 3; ++j ) {
      v[ j ] = sqrta * m_pVel[ 3 * i + j ];// + Hza * m_pPos[ 3 * i + j ];
      r[ j ] = m_Cosm.asnap * m_pPos[ 3 * i + j ];
    }

    J[ 0 ] += ( r[ 1 ] * v[ 2 ] - r[ 2 ] * v[ 1 ] ) * m_pMass[ i ];
    J[ 1 ] += ( r[ 2 ] * v[ 0 ] - r[ 0 ] * v[ 2 ] ) * m_pMass[ i ];
    J[ 2 ] += ( r[ 0 ] * v[ 1 ] - r[ 1 ] * v[ 0 ] ) * m_pMass[ i ];

  }
  Jtot = sqrt( J[ 0 ] * J[ 0 ] + J[ 1 ] * J[ 1 ] + J[ 2 ] * J[ 2 ] );
}

/**
 * computes the angular momentum in internal units. assumes that coordinates and velocities are centered on a reference point
 * @return real Jtot, total angular momentum
 */
real 
HaloProperties::ComputeAngularMomentum( Vector3& J ) const {
  // real Jtot = 0.0;
  Vector3 r;
  Vector3 v;
  // real Hza = m_Cosm.H0 * m_Cosm.asnap * sqrt( m_Cosm.OmegaM0 / ( m_Cosm.asnap * m_Cosm.asnap * m_Cosm.asnap ) + m_Cosm.OmegaL0 );
  real sqrta = sqrt( m_Cosm.asnap );
  J( 0 ) = J( 1 ) = J( 2 ) = 0.0;
  for ( unsigned i = 0; i < m_nPart; ++i ) {
    for ( unsigned j = 0; j < 3; ++j ) {
      v( j ) = sqrta * m_pVel[ 3 * i + j ];// + Hza * m_pPos[ 3 * i + j ]; //HubbleFlow
      r( j ) = m_Cosm.asnap * m_pPos[ 3 * i + j ];
    }
    J += (r^v) * m_pMass[i];
    // J( 0 ) += ( r[ 1 ] * v[ 2 ] - r[ 2 ] * v[ 1 ] ) * m_pMass[ i ];
    // J( 1 ) += ( r[ 2 ] * v[ 0 ] - r[ 0 ] * v[ 2 ] ) * m_pMass[ i ];
    // J( 2 ) += ( r[ 0 ] * v[ 1 ] - r[ 1 ] * v[ 0 ] ) * m_pMass[ i ];

  }
  //Jtot = sqrt( J( 0 ) * J( 0 ) + J( 1 ) * J( 1 ) + J( 2 ) * J( 2 ) );

  return real(J.norm());
}


/**
 * computes the velocity dispersion (sigma) in internal units
 * assumes that pos and vel be centered
 */
void
HaloProperties::ComputeVelocityDispersion( real& Sigma ){
  // either go on with norm of velocity vector
  // use method minimizing roundoff errors
  real a = 0;
  real q = 0;
  for(unsigned i = 0; i<m_nPart; ++i){
    real vx = m_pVel[3*i+0];
    real vy = m_pVel[3*i+1];
    real vz = m_pVel[3*i+2];
    real v2 = sqrt(vx*vx + vy*vy + vz*vz);

    q = q + (1.0*(i-1))/i*(v2-a);
    a = a + (v2-a)/i;
  }
  Sigma = sqrt(q/m_nPart);

  // or step through all euclidean components
  // or through (r,phi,theta) components
  // or use other definition of sigma
}

void 
HaloProperties::GetSortedDistances( std::vector<Distance>& dvec ) const {
  dvec.clear();

  for ( unsigned i = 0, i3 = 0; i < m_nPart; ++i, i3 = 3*i ) {
    real ds2 = 0.0;
    for ( unsigned j = 0; j < 3 ;++j )
      ds2 += m_pPos[ i3 + j ] * m_pPos[ i3 + j ];
    dvec.push_back( Distance( sqrt( ds2 ), i ) );
  }

  std::sort( dvec.begin(), dvec.end() );
}

// Centering in position is required before calling ComputeMajorAxes.
void 
HaloProperties::ComputeMajorAxes( Vector3& Eigenvalues, std::vector<Vector3>& Eigenvectors ) const{
  Matrix33 I; I = 0.0;
  std::vector<real> d;

  for ( unsigned i = 0, i3 = 0; i < m_nPart; ++i, i3 = 3 * i ) {

    I( 0, 0 ) += m_pMass[i] * ( m_pPos[ i3 + 1 ] * m_pPos[ i3 + 1 ] + m_pPos[ i3 + 2 ] * m_pPos[ i3 + 2 ]);
    I( 0, 1 ) += m_pMass[i] * (-m_pPos[ i3 + 0 ] * m_pPos[ i3 + 1 ]);
    I( 0, 2 ) += m_pMass[i] * (-m_pPos[ i3 + 0 ] * m_pPos[ i3 + 2 ]);

    I( 1, 0 ) += m_pMass[i] * (-m_pPos[ i3 + 1 ] * m_pPos[ i3 + 0 ]);
    I( 1, 1 ) += m_pMass[i] * ( m_pPos[ i3 + 0 ] * m_pPos[ i3 + 0 ] + m_pPos[ i3 + 2 ] * m_pPos[ i3 + 2 ]);
    I( 1, 2 ) += m_pMass[i] * (-m_pPos[ i3 + 1 ] * m_pPos[ i3 + 2 ]);

    I( 2, 0 ) += m_pMass[i] * (-m_pPos[ i3 + 2 ] * m_pPos[ i3 + 0 ]);
    I( 2, 1 ) += m_pMass[i] * (-m_pPos[ i3 + 2 ] * m_pPos[ i3 + 1 ]);
    I( 2, 2 ) += m_pMass[i] * ( m_pPos[ i3 + 0 ] * m_pPos[ i3 + 0 ] + m_pPos[ i3 + 1 ] * m_pPos[ i3 + 1 ]);
  }

  I.Eigen( Eigenvalues, Eigenvectors );

  real la = Eigenvalues( 0 );
  real lb = Eigenvalues( 1 );
  real lc = Eigenvalues( 2 );
  if(m_nPart == 0)
    std::cerr << " No particles in halo! Axes not defined! Returning NaNs for Eigenvalues!" << std::endl;
  //assert(m_mtot>0.0);
  Eigenvalues( 0 ) = sqrt( 5.0 / ( 2.0 * m_mtot ) * ( -la + lb + lc ) );
  Eigenvalues( 1 ) = sqrt( 5.0 / ( 2.0 * m_mtot ) * ( + la - lb + lc ) );
  Eigenvalues( 2 ) = sqrt( 5.0 / ( 2.0 * m_mtot ) * ( + la + lb - lc ) );
}



void 
HaloProperties::ComputeMajorAxesReduced( Vector3& Eigenvalues, std::vector<Vector3>& Eigenvectors ) const {
  Matrix33 I;

  I = 0.0;
  std::vector<real> d;
  real mtot = 0.0;

  for ( unsigned j = 0; j < m_nPart; ++j ) {
    mtot += m_pMass[j];

    real r2 = ( m_pPos[ 3 * j + 0 ] * m_pPos[ 3 * j + 0 ] + m_pPos[ 3 * j + 1 ] * m_pPos[ 3 * j + 1 ] + m_pPos[ 3 * j + 2 ] * m_pPos[ 3 * j + 2 ] );

    if ( j != m_imostbound ) {
      I( 0, 0 ) += m_pMass[j] * (  m_pPos[ 3 * j + 1 ] * m_pPos[ 3 * j + 1 ] + m_pPos[ 3 * j + 2 ] * m_pPos[ 3 * j + 2 ] ) / r2;
      I( 0, 1 ) += m_pMass[j] * ( -m_pPos[ 3 * j + 0 ] * m_pPos[ 3 * j + 1 ] ) / r2;
      I( 0, 2 ) += m_pMass[j] * ( -m_pPos[ 3 * j + 0 ] * m_pPos[ 3 * j + 2 ] ) / r2;

      I( 1, 0 ) += m_pMass[j] * ( -m_pPos[ 3 * j + 1 ] * m_pPos[ 3 * j + 0 ] ) / r2;
      I( 1, 1 ) += m_pMass[j] * (  m_pPos[ 3 * j + 0 ] * m_pPos[ 3 * j + 0 ] + m_pPos[ 3 * j + 2 ] * m_pPos[ 3 * j + 2 ] ) / r2;
      I( 1, 2 ) += m_pMass[j] * ( -m_pPos[ 3 * j + 1 ] * m_pPos[ 3 * j + 2 ] ) / r2;

      I( 2, 0 ) += m_pMass[j] * ( -m_pPos[ 3 * j + 2 ] * m_pPos[ 3 * j + 0 ] ) / r2;
      I( 2, 1 ) += m_pMass[j] * ( -m_pPos[ 3 * j + 2 ] * m_pPos[ 3 * j + 1 ] ) / r2;
      I( 2, 2 ) += m_pMass[j] * (  m_pPos[ 3 * j + 0 ] * m_pPos[ 3 * j + 0 ] + m_pPos[ 3 * j + 1 ] * m_pPos[ 3 * j + 1 ] ) / r2;
    }
  }

  I.Eigen( Eigenvalues, Eigenvectors );

  real la = Eigenvalues( 0 );
  real lb = Eigenvalues( 1 );
  real lc = Eigenvalues( 2 );

  Eigenvalues( 0 ) = sqrt( (5.0 / ( 2.0 * mtot )) * ( - la + lb + lc ) );
  Eigenvalues( 1 ) = sqrt( (5.0 / ( 2.0 * mtot )) * ( + la - lb + lc ) );
  Eigenvalues( 2 ) = sqrt( (5.0 / ( 2.0 * mtot )) * ( + la + lb - lc ) );
}

// computed relative to coordinate origin
void 
HaloProperties::ComputeShapeParameters( real& T, real& S ) const {
  Vector3 EW;
  std::vector<Vector3 > EV;
  ComputeMajorAxes( EW, EV );

  T = ( EW(0) * EW(0) - EW(1) * EW(1) ) / ( EW(0) * EW(0) - EW(2) * EW(2) );
  S = EW(2) / EW(0);
}

void 
HaloProperties::ComputeVirialProperties200( real rhobox, real& Mvir, real& Rvir, real& Vcmax, std::vector<Distance> &dvec ) const {
  GetSortedDistances( dvec );

  real rhomean = 1.0e+99;
  unsigned ninside = 1;
  real R = 0.0;
  real pi43 = 4.0 * M_PI / 3.0;
  real fz =  m_Cosm.OmegaM0 / ( m_Cosm.asnap * m_Cosm.asnap * m_Cosm.asnap ) /
    ( m_Cosm.OmegaM0 / ( m_Cosm.asnap * m_Cosm.asnap * m_Cosm.asnap )
      + m_Cosm.OmegaL0 ) - 1.0;

  real rhocrit = rhobox / m_Cosm.OmegaM0 * m_Cosm.asnap * m_Cosm.asnap * m_Cosm.asnap;
  // mean density / critical density (e.g. Bryan, Norman 1998)
  real delta = ( 18.0 * M_PI * M_PI + 82.0 * fz - 39.0 * fz * fz );
  real rhodelta = delta * rhocrit;

  real Vc;
  Vcmax = -1.0e+99;


  rhodelta = 200.0 * rhobox;
  //    std::cerr << "delta = " << rhodelta << ", overdens = " << rhodelta/rhobox << std::endl;

  real minside = 0.0;

  while ( rhomean > rhodelta && ninside < dvec.size() - 1 ) {
    R = m_Cosm.asnap * dvec[ ninside ].getr();
    minside += m_pMass[ dvec[ ninside ].geti() ];
    rhomean = minside / ( pi43 * R * R * R );

    Vc = minside / R;
    if ( Vc > Vcmax )
      Vcmax = Vc;
    ++ninside;
  }

  Vcmax = Vcmax * m_Cosm.G;
  Mvir = minside;
  Rvir = R;
}

void 
HaloProperties::ComputeVirialProperties200( real rhobox, real& Mvir, real& Rvir, real &Vcmax ) const {
  std::vector<Distance> dvec;
  ComputeVirialProperties200( rhobox, Mvir, Rvir, Vcmax, dvec );
}
void 
HaloProperties::ComputeVirialProperties200( real rhobox, real& Mvir, real& Rvir, std::vector<Distance> &dvec ) const {
  real Vcmax;
  ComputeVirialProperties200( rhobox, Mvir, Rvir, Vcmax, dvec );
}
void 
HaloProperties::ComputeVirialProperties200( real rhobox, real& Mvir, real& Rvir ) const {
  std::vector<Distance> dvec;
  real Vcmax;
  ComputeVirialProperties200( rhobox, Mvir, Rvir, Vcmax, dvec );
}


void 
HaloProperties::ComputeVirialProperties( real rhobox, real& Mvir, real& Rvir, real& Vcmax, std::vector<Distance> &dvec ) const {
  GetSortedDistances( dvec );

  real rhomean = 1.0e+99;
  unsigned ninside = 1;
  real R = 0.0;
  real pi43 = 4.0 * M_PI / 3.0;
  // real fz = m_Cosm.OmegaM0 / ( m_Cosm.asnap * m_Cosm.asnap * m_Cosm.asnap ) / ( m_Cosm.OmegaM0 / ( m_Cosm.asnap * m_Cosm.asnap * m_Cosm.asnap ) + m_Cosm.OmegaL0 ) - 1.0; // see below
  // real rhocrit = rhobox / m_Cosm.OmegaM0 * m_Cosm.asnap * m_Cosm.asnap * m_Cosm.asnap; // see below
  // real delta = ( 18.0 * M_PI * M_PI + 82.0 * fz - 39.0 * fz * fz ); // see below
  // real rhodelta = delta * rhocrit; //used only if overdensity criterion is ON

  real Vc;
  Vcmax = -1.0e+99;
  real minside = 0.0;
#warning Disabled checking for maximum radius based on overdensity criterion
  while ( /*rhomean > rhodelta &&*/ ninside < dvec.size() - 1 ) {

    R = m_Cosm.asnap * dvec[ ninside ].getr();
    minside += m_pMass[ dvec[ ninside ].geti() ];
    rhomean = minside / ( pi43 * R * R * R );

    Vc = minside / R;
    if ( Vc > Vcmax )
      Vcmax = Vc;
    ++ninside;
  }

  Vcmax = Vcmax * m_Cosm.G;
  Mvir = minside;
  Rvir = R;
}
void 
HaloProperties::ComputeVirialProperties( real rhobox, real& Mvir, real& Rvir, std::vector<Distance> &dvec ) const {
  real Vcmax;
  ComputeVirialProperties( rhobox, Mvir, Rvir, Vcmax, dvec );
}
void 
HaloProperties::ComputeVirialProperties( real rhobox, real& Mvir, real& Rvir, real& Vcmax ) const {
  std::vector<Distance> dvec;
  ComputeVirialProperties( rhobox, Mvir, Rvir, Vcmax, dvec );
}
void 
HaloProperties::ComputeVirialProperties( real rhobox, real& Mvir, real& Rvir ) const {
  std::vector<Distance> dvec;
  real Vcmax;
  ComputeVirialProperties( rhobox, Mvir, Rvir, Vcmax, dvec );
}

real 
HaloProperties::ComputePoormanConcentration( real rhobox, real &R20, real &R80 ) {
  std::vector<Distance> dvec;
  real Mvir, Rvir;
  ComputeVirialProperties( rhobox, Mvir, Rvir, dvec );

  R20 = dvec[ ( unsigned ) ( 0.2 * ( real ) dvec.size() + 0.5 ) ].getr() * m_Cosm.asnap;
  R80 = dvec[ ( unsigned ) ( 0.8 * ( real ) dvec.size() + 0.5 ) ].getr() * m_Cosm.asnap;

  return R80 / R20;
}
real 
HaloProperties::ComputePoormanConcentration( real rhobox ) {
  real R20, R80;
  return ComputePoormanConcentration( rhobox, R20, R80 );
}

/**
 * \brief computes the simplified spin parameter lambda prime
 * computes the lambda prime (Bullock et al. 2001). ASSUMES THAT COORDINATES
 * AND VELOCITIES ARE CENTERED ON REFERENCE POINT. so center on center of mass
 * or something.
 * @param lambda storage for halo spin parameter
 * @param Mvir the virial mass is returned here
 * @param Rvir the virial radius is returned here
 * @return void
 */
void 
HaloProperties::ComputeLambdaPrime( real &lambda, real Mvir, real Rvir ) {
  std::vector<Distance> dvec;
  GetSortedDistances(dvec);

  Vector3 J = 0.0;
  Vector3 d = 0.0;
  Vector3 v = 0.0;
  unsigned ninside = 0;
  real minside = 0.0;
  real sqrta = sqrt( m_Cosm.asnap );

  //for( unsigned i=0; i<10 && i<m_nPart; ++i )
  //	std::cout << m_Cosm.asnap*(dvec[i].getr()) << " " << m_pPos[3*dvec[i].geti()] << " " << Rvir << " mass: "<< m_pMass[ dvec[ ninside ].geti() ] << std::endl;
  while ( ninside < m_nPart /*&& m_Cosm.asnap * dvec[ ninside ].getr() < Rvir*/) {
    //std::cout << ninside << std::endl;
    for ( unsigned j = 0; j < 3; ++j ) {
      v( j ) = sqrta * m_pVel[ 3 * dvec[ ninside ].geti() + j ];
      d( j ) = m_Cosm.asnap * m_pPos[ 3 * dvec[ ninside ].geti() + j ];
    }

    //std::cout <<  m_pMass[ dvec[ ninside ].geti() ] << std::endl;
    J += (d^v)*m_pMass[ dvec[ ninside ].geti() ];

    minside += m_pMass[ dvec[ ninside ].geti() ];
    ++ninside;
  }
  std::cerr << " Mvir = " << Mvir << ", Rvir = " << Rvir << std::endl;
  assert( Mvir > 0.0 );
  assert( Rvir > 0.0 );
  if( ninside == 0 )
    lambda = -1.0;

  lambda = real( J.norm() / ( sqrt( 2.0 * m_Cosm.G * Rvir * Mvir ) * minside ) );
}

/**
 * Compute the kinetic energies for all particles
 * @param kin vector to be filled with kin. energies
 * @return total kinetic energy
 */
real 
HaloProperties::ComputeKineticEnergy( std::vector<real>& kin ) const {
  real kintot = 0.0;
  real v[ 3 ];

  real Hza = m_Cosm.H0 * m_Cosm.asnap * sqrt( m_Cosm.OmegaM0 / ( m_Cosm.asnap * m_Cosm.asnap * m_Cosm.asnap ) + m_Cosm.OmegaL0 );
  real sqrta = sqrt( m_Cosm.asnap );

  kin.clear();

  if ( m_pVel != NULL ) {
    for ( unsigned i = 0; i < m_nPart; ++i ) {
      real k = 0;
      for ( unsigned j = 0; j < 3; ++j ) {
	v[ j ] = sqrta * m_pVel[ 3 * i + j ] + Hza * m_pPos[ 3 * i + j ];
	k += v[ j ] * v[ j ];
      }

      k *= 0.5 * m_pMass[ i ];

      kin.push_back( ( real ) k );
      kintot += k;
    }
  } else {
    static bool bvelignore( false );
    if ( !bvelignore ) {
      std::cerr << " - Warning in HaloProperties::ComputeKineticEnergy" << std::endl;
      std::cerr << "     Setting kinetic energy to zero since I do not have velocities available!" << std::endl;
      bvelignore = true;
    }
    for ( unsigned i = 0; i < m_nPart; ++i )
      kin.push_back( 0.0 );
  }
  return real( kintot );
}

/**
 * Computes the kinetic energy of all particles without Hubble flow.
 * @param kin vector with kinetic energies
 * @return total kinetic energy
 */
real 
HaloProperties::ComputeKineticEnergyNoHubble( std::vector<real>& kin ) const {
  real kintot = 0.0;
  real v[ 3 ];

  real sqrta = sqrt( m_Cosm.asnap );
  kin.clear();

  // need for memory
  kin.reserve(m_nPart);
  if ( m_pVel != NULL ) {
    for ( unsigned i = 0; i < m_nPart; ++i ) {
      real k = 0;
      for ( unsigned j = 0; j < 3; ++j ) {
	v[ j ] = sqrta * m_pVel[ 3 * i + j ];
	k += v[ j ] * v[ j ];
      }

      k *= m_pMass[ i ] * 0.5;

      kin.push_back( real(k) );
      kintot += k;
    }
  } else {
    static bool bvelignore( false );
    if ( !bvelignore ) {
      std::cerr << " - Warning in HaloProperties::ComputeKineticEnergyNoHubble" << std::endl;
      std::cerr << "     Setting kinetic energy to zero since I do not have velocities available!" << std::endl;
      bvelignore = true;
    }
    kin.assign(m_nPart,0.0);
    //for ( unsigned i = 0; i < m_nPart; ++i )
    //kin.at( i ) = 0.0;
  }
  return real( kintot );
}

void 
HaloProperties::ComputeKineticEnergyNoHubbleExc( std::vector<real>& kin, std::vector<bool> excluded ) const {
  for ( unsigned i = 0, i3 = 0; i < m_nPart; ++i, i3 = 3 * i ) {
    if ( excluded.at( i ) )
      kin.at(i) = 0.0;
    else
      kin.at( i ) = m_Cosm.asnap * ( m_pVel[ i3 + 0 ] * m_pVel[ i3 + 0 ] + m_pVel[ i3 + 1 ] * m_pVel[ i3 + 1 ] + m_pVel[ i3 + 2 ] * m_pVel[ i3 + 2 ] ) * m_pMass[ i ] / 2.0;
  }
}

/**
 * Computes the potential directly from particle positions and masses. Center on CM before invoking this routine
 * @param plumsoft smoothing length
 * @param phi potential
 * @return total potential energy
 */
real 
HaloProperties::ComputePotentialDirect( real plumsoft, std::vector<real>& phi ) const {
  real phitot = 0.0;
  plumsoft *= m_Cosm.asnap;

  for ( unsigned i = 0, i3 = 0; i < m_nPart; ++i, i3 = 3 * i ) {
    //std::cout << "i = " << i << " of " << m_nPart << std::endl;
    real phii = 0.0;
    for ( unsigned j = i + 1, j3 = 0; j < m_nPart; ++j, j3 = 3 * j ) {
      //for ( unsigned j = 0; j < m_nPart; ++j ) {
      real
	dx = m_Cosm.asnap * ( m_pPos[ i3 + 0 ] - m_pPos[ j3 + 0 ] ),
	dy = m_Cosm.asnap * ( m_pPos[ i3 + 1 ] - m_pPos[ j3 + 1 ] ),
	dz = m_Cosm.asnap * ( m_pPos[ i3 + 2 ] - m_pPos[ j3 + 2 ] );

      real r2 = dx * dx + dy * dy + dz * dz;
      if ( j != i )
	phii += -m_Cosm.G * m_pMass[ i ] * m_pMass[ i ] / ( sqrt( r2 ) + plumsoft );
    }

    phi.push_back( ( real ) phii );
    phitot += phii;
  }
  return real( phitot );
}

void 
HaloProperties::ComputePotentialDirectExc( real plumsoft, std::vector<real>& phi, std::vector<bool> excluded ) const {
  plumsoft *= m_Cosm.asnap;

  for ( unsigned i = 0, i3 = 0; i < m_nPart; ++i, i3 = 3 * i ) {
    //std::cout << "i = " << i << " of " << m_nPart << std::endl;
    //for all particles needed. so don't check for excluded here, but only for j
    // really???
    if( excluded.at(i) )
      continue;
    real dx = 0.0;
    real dy = 0.0;
    real dz = 0.0;
    real phii = 0.0;

    for ( unsigned j = 0, j3 = 0; j < m_nPart; ++j, j3 = 3 * j ) {
      if ( excluded.at( j ) )// || i == j )
	continue;
      j3 = 3 * j;
      dx = m_Cosm.asnap * ( m_pPos[ i3 + 0 ] - m_pPos[ j3 + 0 ] );
      dy = m_Cosm.asnap * ( m_pPos[ i3 + 1 ] - m_pPos[ j3 + 1 ] );
      dz = m_Cosm.asnap * ( m_pPos[ i3 + 2 ] - m_pPos[ j3 + 2 ] );

      phii += -m_Cosm.G * m_pMass[ i ] * m_pMass[ j ] / ( sqrt( dx * dx + dy * dy + dz * dz ) + plumsoft );
    }
    //std::cout << " phi(" << i << ") = " << phii << std::endl;
    phi.at( i ) = real( phii );
  }
}
real 
HaloProperties::ComputePotentialTree( real theta, real plumsoft, std::vector<real>& phi ) const {
  // 0.7, 0,
  phi.clear();
  plumsoft *= m_Cosm.asnap;

  real phitot = 0.0;

  std::vector< PseudoBody > barray;

  real xmin = 1.0e+99; //real xmax = -1.0e+99;
  real ymin = 1.0e+99; //real ymax = -1.0e+99;
  real zmin = 1.0e+99; //real zmax = -1.0e+99;

  for ( unsigned i = 0, i3 = 0; i < m_nPart; ++i, i3 = 3 * i ) {
    xmin = ( m_pPos[ i3 + 0 ] < xmin ) ? m_pPos[ i3 + 0 ] : xmin;
    //xmax = ( m_pPos[ i3 + 0 ] > xmax ) ? m_pPos[ i3 + 0 ] : xmax;
    ymin = ( m_pPos[ i3 + 1 ] < ymin ) ? m_pPos[ i3 + 1 ] : ymin;
    //ymax = ( m_pPos[ i3 + 1 ] > ymax ) ? m_pPos[ i3 + 1 ] : ymax;
    zmin = ( m_pPos[ i3 + 2 ] < zmin ) ? m_pPos[ i3 + 2 ] : zmin;
    //zmax = ( m_pPos[ i3 + 2 ] > zmax ) ? m_pPos[ i3 + 2 ] : zmax;
  }

  real x, y, z;
  for ( unsigned i = 0, i3 = 0; i < m_nPart; ++i, i3 = 3 * i ) {
    x = m_Cosm.asnap * ( m_pPos[ i3 + 0 ] - xmin );
    y = m_Cosm.asnap * ( m_pPos[ i3 + 1 ] - ymin );
    z = m_Cosm.asnap * ( m_pPos[ i3 + 2 ] - zmin );

    barray.push_back( PseudoBody( m_pMass[ i ], x, y, z ) );
  }

  //real length = std::max( std::max(xmax-xmin,ymax-ymin),zmax-zmin );

  //#ifndef OCTTREE
  //	kDTreeNode *BHroot = new kDTreeNode( barray );
  //#else
  OctTreeNode *BHroot = new OctTreeNode( barray, plumsoft ); //, length );
  //#endif


  for ( unsigned i = 0; i < m_nPart; ++i ) {
    real phii = m_Cosm.G * m_pMass[ i ] * BHroot->traverse( barray[ i ].getx(),
							    barray[ i ].gety(),
							    barray[ i ].getz(),
							    theta, plumsoft );
    phitot += phii;
    phi.push_back( real( phii ) );
  }

  delete BHroot;

  return 0.5*phitot;
}
void 
HaloProperties::ComputePotentialTreeExc( real theta, real plumsoft, std::vector<real>& phi, std::vector<bool> excluded ) const {
  plumsoft *= m_Cosm.asnap;

  std::vector< PseudoBody > barray;

  real xmin = 1.0e+99; //real xmax = -1.0e+99;
  real ymin = 1.0e+99; //real ymax = -1.0e+99;
  real zmin = 1.0e+99; //real zmax = -1.0e+99;

  // recenter all particles to minima of coordinates
  for ( unsigned i = 0, i3 = 0; i < m_nPart; ++i, i3 = 3 * i ) {
    xmin = ( m_pPos[ i3 + 0 ] < xmin ) ? m_pPos[ i3 + 0 ] : xmin;
    //xmax = ( m_pPos[ i3 + 0 ] > xmax ) ? m_pPos[ i3 + 0 ] : xmax;
    ymin = ( m_pPos[ i3 + 1 ] < ymin ) ? m_pPos[ i3 + 1 ] : ymin;
    //ymax = ( m_pPos[ i3 + 1 ] > ymax ) ? m_pPos[ i3 + 1 ] : ymax;
    zmin = ( m_pPos[ i3 + 2 ] < zmin ) ? m_pPos[ i3 + 2 ] : zmin;
    //zmax = ( m_pPos[ i3 + 2 ] > zmax ) ? m_pPos[ i3 + 2 ] : zmax;
  }

  // write masses and positions into tree
  real x, y, z;
  for ( unsigned i = 0, i3 = 0; i < m_nPart; ++i, i3 = 3 * i ) {
    x = m_Cosm.asnap * ( m_pPos[ i3 + 0 ] - xmin );
    y = m_Cosm.asnap * ( m_pPos[ i3 + 1 ] - ymin );
    z = m_Cosm.asnap * ( m_pPos[ i3 + 2 ] - zmin );
    if ( excluded.at( i ) )
      barray.push_back( PseudoBody( 0.0, x, y, z ) );
    else
      barray.push_back( PseudoBody( m_pMass[ i ], x, y, z ) );
  }

  //#ifndef OCTTREE
  //	kDTreeNode *BHroot = new kDTreeNode( barray );
  //#else
  OctTreeNode *BHroot = new OctTreeNode( barray, plumsoft ); //, length );
  //#endif

  // compute potential for all particles
  for ( unsigned i = 0; i < m_nPart; ++i ) {
    real phii = m_Cosm.G * m_pMass[ i ] * BHroot->traverse( barray[ i ].getx(),
							    barray[ i ].gety(),
							    barray[ i ].getz(),
							    theta, plumsoft );
    phi.at( i ) = real( phii );
  }
  delete BHroot;
}

std::vector<bool> 
HaloProperties::CleanUnbound( real boxsize ) {
  std::vector<bool> excluded( m_nPart );
  CleanUnboundExc( boxsize, excluded );

  return excluded;
}

void 
HaloProperties::CleanUnboundExc( real boxsize, std::vector<bool>& excluded ) {
  //std::cout << "Cleanup started with " << m_nPart << " particles...";
  //std::cout.flush();

  unsigned pc = 0;
  for ( unsigned j = 0; j < m_nPart; ++j ) {
    if ( m_nPart > 1000 && j % 10 == 1 )
      std::cout << "j = " << j << std::endl;
    std::vector<real> phi( m_nPart );
    std::vector<real> kin( m_nPart );

    //std::cout << "Centering particles...";
    CenterOnCMExc( boxsize, excluded );
    //std::cout << "and velocities...";
    CenterVelExc( excluded );
    //std::cout << "done!" << std::endl;

    //ComputePotentialDirectExc( 0.0, phi, excluded );
    //ComputePotentialTreeExc( theta, plumsoft, phi );
    // theta gives "angle" to leaves
    real theta = 0.7;
    ComputePotentialTreeExc( theta, 1.0, phi, excluded );
    //std::cout << "potential computed" << std::endl;
    ComputeKineticEnergyNoHubbleExc( kin, excluded );
    //std::cout << "kinetic energy computed" << std::endl;

    real emax = -1E99;
    int maxpos = -1;
    //get most unbound particle
    for ( unsigned i = 0; i < m_nPart; ++i ) {
      if ( excluded.at( i ) )
	continue;
      //std::cout << "kin[" << i << "]=" << kin[ i ] << ", phi[" << i << "]=" << phi[ i ] << std::endl;

      if ( phi[ i ] + kin[ i ] < emax )
	continue;
      maxpos = i;
      //std::cout << "maxpos = " << maxpos << std::endl;

      emax = phi[ i ] + kin[ i ];
    }

    if ( emax < 0 || pc + 1 == m_nPart ) {
      // std::cout << "no more unbound particles found" << std::endl;
      break;
    } else {
      //std::cout << "removing particle " << maxpos << " from list..." << std::endl;
      excluded.at( maxpos ) = true;
      ++pc;
    }

    // reinclude bound particles that were excluded previously
    for ( unsigned i = 0; i < m_nPart; ++i ) {
      if ( excluded.at( i ) && phi[ i ] + kin[ i ] < 0 ) {
	excluded.at( i ) = false;
	--pc;
	std::cout << i << " is bound again..." << std::endl;
      }
    }
  }
  std::cout << "	Removed " << pc << " out of " << m_nPart << " particles." << std::endl;
}

void 
HaloProperties::CleanUnboundSimple( real boxsize, std::vector<bool>& excluded ) {
  int pc = 0;
  std::vector<real> phi( m_nPart );
  std::vector<real> kin( m_nPart );

  //std::cout << "Centering particles...";
  CenterOnCM( boxsize );
  //std::cout << "and velocities...";
  CenterVel( );
  //std::cout << "done!" << std::endl;

  ComputePotentialDirectExc( 0.0, phi, excluded );
  //std::cout << "potential computed" << std::endl;
  ComputeKineticEnergyNoHubbleExc( kin, excluded );
  //std::cout << "kinetic energy computed" << std::endl;
  for ( unsigned i = 0; i < m_nPart; ++i ) {
    if ( phi[ i ] + kin[ i ] > 0 ) {
      //std::cout << "removing particle " << maxpos << " from list..." << std::endl;
      excluded.at( i ) = true;
      ++pc;
    }
  }
  std::cout << " * * S:Removed " << pc << " out of " << m_nPart << " particles." << std::endl;
}






#if 1
int 
HaloProperties::ComputeDensityProfile( real rhobox, unsigned n,
				       std::vector<real>& r,
				       std::vector<real>& rho,
				       std::vector<real>& sigma_rho ) const {
  std::vector<real> rr;
  real Mvir, Rvir;
  rho.clear();
  std::vector<Distance> dvec;
  ComputeVirialProperties( rhobox, Mvir, Rvir, dvec );

  if ( dvec.size() < 50 ) {
    std::cerr << "not enough particles. npart=" << m_nPart << std::endl;
    std::cerr << "dvec.size() = " << dvec.size() << std::endl;

    return 0;
  }

  // hmm. either - or?
  unsigned nn = unsigned( ( real ) dvec.size() / ( n + 1.0 ) );
  nn = unsigned( ( real ) dvec.size() / n );

  std::vector<real> m, dV;
  std::vector<unsigned> nshell;


  for ( unsigned i = 1;i < n + 1; ++i ) {
    rr.push_back( dvec[ nn * i ].getr() * m_Cosm.asnap );
    m.push_back( nn * m_pMass[ i ] );
    nshell.push_back( nn );
  }


  dV.push_back( 4.0 * M_PI / 3.0 * rr[ 0 ] * rr[ 0 ] * rr[ 0 ] );
  for ( unsigned i = 1; i < n; ++i )
    dV.push_back( 4.0 * M_PI / 3.0 * ( rr[ i ] * rr[ i ] * rr[ i ] - rr[ i - 1 ] * rr[ i - 1 ] * rr[ i - 1 ] ) );

  for ( unsigned i = 0; i < n; ++i ) {
    real rhonext = m[ i ] / ( dV[ i ] * m_Cosm.asnap * m_Cosm.asnap * m_Cosm.asnap );
    if ( i > 1 && rhonext > 1.25 * rho[ i - 1 ] )
      break;

    rho.push_back( rhonext );
    sigma_rho.push_back( rhonext / sqrt( nshell[ i ] ) );
    r.push_back( rr[ i ] );

  }

  return 1;
}
#else // now 0
int 
HaloProperties::ComputeDensityProfile( real rhobox, unsigned n,
				       std::vector<real>& r,
				       std::vector<real>& rho,
				       std::vector<real>& sigma_rho ) const {
  std::vector<real> rr;
  real Mvir, Rvir;
  rho.clear();
  std::vector<Distance> dvec;
  ComputeVirialProperties( rhobox, Mvir, Rvir, dvec );

  if ( dvec.size() < 50 ) {
    std::cerr << " - not enough particles. npart=" << m_nPart << std::endl;
    std::cerr << "     dvec.size() = " << dvec.size() << std::endl;
    abort();
  }

  logspace( log10( dvec[ 50 ].getr() ), log10( Rvir / m_Cosm.asnap ), n, rr );

  if ( dvec[ 50 ].getr() > Rvir / m_Cosm.asnap ) {
    std::cerr << "  Found strange profile, will skip this step." << std::endl;
    return 0;
  }

  std::vector<real> m;
  std::vector<real> dV;
  std::vector<unsigned> nshell;
  m.assign( n, 0.0 );
  nshell.assign( n, 0 );

  real Rvirp = ( 1.0 + 1.0e-6 ) * Rvir;

  unsigned ishell = 0, ninside = 0;
  while ( dvec[ ninside ].getr() <= Rvirp && ninside < m_nPart ) {
    if ( dvec[ ninside ].getr() > rr[ ishell ] )
      ++ishell;
    m[ ishell ] += m_pMass[ ninside ];
    ++nshell[ ishell ];
    ++ninside;
  }

  dV.push_back( 4.0 * M_PI / 3.0 * rr[ 0 ] * rr[ 0 ] * rr[ 0 ] );
  for ( unsigned i = 1; i < n; ++i )
    dV.push_back( 4.0 * M_PI / 3.0 * ( rr[ i ] * rr[ i ] * rr[ i ] - rr[ i - 1 ] * rr[ i - 1 ] * rr[ i - 1 ] ) );

  for ( unsigned i = 0; i < n; ++i ) {
    real rhonext = m[ i ] / ( dV[ i ] * m_Cosm.asnap * m_Cosm.asnap * m_Cosm.asnap );
    if ( i > 0 && rhonext > 1.1 * rho[ i - 1 ] )
      break;

    rho.push_back( rhonext );
    sigma_rho.push_back( rho[ i ] / sqrt( nshell[ i ] ) );
    r.push_back( rr[ i ] );


  }

  if ( r.size() < n / 2 )
    return 0;

  return 1;
}
#endif

int 
HaloProperties::FitNFWProfile2( real rhobox, unsigned n, real& fdeltac, real& fedeltac, real& frs, real& fers ) const {
  std::vector<real> r, rho, erho;
  real dc0, rs0, dc, edc, rs, ers;

  if ( !ComputeDensityProfile( rhobox, n, r, rho, erho ) ) {
    std::cerr << " - no convergence in FitNFW!" << std::endl;
    return 0;
  }

  dc0 = rho[ 0 ] * 10.0;
  rs0 = r[ r.size() - 1 ] * 0.1;
  NFW_fit( r, 
	   rho, 
	   erho, 
	   dc0, 
	   rs0, 
	   dc, 
	   edc,
	   rs,
	   ers );

  fdeltac = dc;
  frs = rs;

  return 1;
}

int 
HaloProperties::FitNFWProfile1( real rhobox, unsigned n, real& frs, real& fers ) const {
  std::vector<real> r, rho, erho;
  real rs0;
  real rs;
  real ers;
  real mvir;
  real rvir;

  ComputeVirialProperties( rhobox, mvir, rvir );

  if ( !ComputeDensityProfile( rhobox, n, r, rho, erho ) ) {
    std::cerr << " - no convergence in FitNFW!" << std::endl;
    return 0;
  }

  rs0 = r[ r.size() - 1 ] * 0.1;
  NFW_fix_fit( r, rho, erho, rvir, mvir, rs0, rs, ers );

  frs = rs;

  return 1;
}
int 
HaloProperties::FitNFWProfile200( real rhobox, unsigned n, real& frs, real& fers ) const {
  std::vector<real> r, rho, erho;
  real rs0;
  real rs;
  real ers;
  real mvir;
  real rvir;

  ComputeVirialProperties200( rhobox, mvir, rvir );

  if ( !ComputeDensityProfile( rhobox, n, r, rho, erho ) ) {
    std::cerr << " - no convergence in FitNFW!" << std::endl;
    return 0;
  }

  rs0 = r[ r.size() - 1 ] * 0.1;
  NFW_fix_fit( r, rho, erho, rvir, mvir, rs0, rs, ers );

  frs = rs;

  return 1;
}
