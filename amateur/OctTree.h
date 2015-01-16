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
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/**
 * \file OctTree.h
 * \brief Build trees as new class.
 */

#ifndef OCTTREE_H
#define OCTTREE_H

#include <vector>
#include <iostream>
#include <fstream>

#include <stack>
#include <queue>
#include <cmath>

#include "Global.h"

#define EPS 1.0e-8

extern std::ofstream *ofs;
const real g_geomfac = 2.0 / sqrt( 3.0 );

inline void 
SPLINEQ( double invr, double r2, double twoh, double& a, double& b,
	 double& c, double& d ) {
  double u, dih, dir = ( invr );
  if ( ( r2 ) < ( twoh ) * ( twoh ) ) {
    dih = 2.0 / ( twoh );
    u = dih / dir;
    if ( u < 1.0 ) {
      a = dih * ( 7.0 / 5.0 - 2.0 / 3.0 * u * u + 3.0 / 10.0 * u * u * u * u
		  - 1.0 / 10.0 * u * u * u * u * u );
      b = dih * dih * dih * ( 4.0 / 3.0 - 6.0 / 5.0 * u * u + 1.0 / 2.0 * u * u * u );
      c = dih * dih * dih * dih * dih * ( 12.0 / 5.0 - 3.0 / 2.0 * u );
      d = 3.0 / 2.0 * dih * dih * dih * dih * dih * dih * dir;
    } else {
      a = -1.0 / 15.0 * dir + dih * ( 8.0 / 5.0 - 4.0 / 3.0 * u * u + u * u * u
				      - 3.0 / 10.0 * u * u * u * u + 1.0 / 30.0 * u * u * u * u * u );
      b = -1.0 / 15.0 * dir * dir * dir + dih * dih * dih * ( 8.0 / 3.0 - 3.0 * u + 6.0 / 5.0 * u * u - 1.0 / 6.0 * u * u * u );
      c = -1.0 / 5.0 * dir * dir * dir * dir * dir + 3.0 * dih * dih * dih * dih * dir
	+ dih * dih * dih * dih * dih * ( -12.0 / 5.0 + 1.0 / 2.0 * u );
      d = -dir * dir * dir * dir * dir * dir * dir
	+ 3.0 * dih * dih * dih * dih * dir * dir * dir
	- 1.0 / 2.0 * dih * dih * dih * dih * dih * dih * dir;
    }
  } else {
    a = dir;
    b = a * a * a;
    c = 3.0 * b * a * a;
    d = 5.0 * c * a * a;
  }
}


class Body {
 protected:
  real mass;
  real x, y, z;

 public:
  Body( void )
    : mass( 0.0 ), x( 0.0 ), y( 0.0 ), z( 0.0 ) { }

  Body( const std::vector<Body>& Bodies )
    : mass( 0.0 ), x( 0.0 ), y( 0.0 ), z( 0.0 ) {
    unsigned n = Bodies.size();

    for ( unsigned i = 0;i < n;++i ) {
      x += Bodies[ i ].mass * Bodies[ i ].x;
      y += Bodies[ i ].mass * Bodies[ i ].y;
      z += Bodies[ i ].mass * Bodies[ i ].z;
      mass += Bodies[ i ].mass;
    }

    x /= mass;
    y /= mass;
    z /= mass;
  }

  Body( const Body& b )
    : mass( b.mass ), x( b.x ), y( b.y ), z( b.z ) { }

  Body( real mass_, real x_, real y_, real z_ )
    : mass( mass_ ), x( x_ ), y( y_ ), z( z_ ) { }


  inline real getx( void ) const {
    return x;
  }
  inline real gety( void ) const {
    return y;
  }
  inline real getz( void ) const {
    return z;
  }
  inline real getmass( void ) const {
    return mass;
  }
};

class PseudoBody {
 protected:
  real mass;
  real x, y, z;
  //real px, py, pz;
  real xx, yy, zz, xy, xz, yz;

 public:
  PseudoBody( void )
    : mass( 0.0 ), x( 0.0 ), y( 0.0 ), z( 0.0 ), xx( 0.0 ), yy( 0.0 ), zz( 0.0 ), xy( 0.0 ), xz( 0.0 ), yz( 0.0 ) { }

  PseudoBody( real mass_, real x_, real y_, real z_ )
    : mass( mass_ ), x( x_ ), y( y_ ), z( z_ ),
    xx( 0.0 ), yy( 0.0 ), zz( 0.0 ), xy( 0.0 ), xz( 0.0 ), yz( 0.0 ) { }

  PseudoBody( const std::vector<PseudoBody>& Bodies )
    : mass( 0.0 ), x( 0.0 ), y( 0.0 ), z( 0.0 ),
    xx( 0.0 ), yy( 0.0 ), zz( 0.0 ), xy( 0.0 ), xz( 0.0 ), yz( 0.0 ) {
    unsigned n = Bodies.size();
    real dx, dy, dz;

    for ( unsigned i = 0;i < n;++i ) {
      x += Bodies[ i ].getmass() * Bodies[ i ].getx();
      y += Bodies[ i ].getmass() * Bodies[ i ].gety();
      z += Bodies[ i ].getmass() * Bodies[ i ].getz();
      mass += Bodies[ i ].getmass();
    }

    x /= mass;
    y /= mass;
    z /= mass;

    for ( unsigned i = 0;i < n;++i ) {

      dx = x - Bodies[ i ].getx();
      dy = y - Bodies[ i ].gety();
      dz = z - Bodies[ i ].getz();

      //real r2 = dx*dx+dy*dy+dz*dz;

      xx += Bodies[ i ].getmass() * dx * dx;
      yy += Bodies[ i ].getmass() * dy * dy;
      zz += Bodies[ i ].getmass() * dz * dz;
      xy += Bodies[ i ].getmass() * dx * dy;
      xz += Bodies[ i ].getmass() * dx * dz;
      yz += Bodies[ i ].getmass() * dy * dz;
    }

  }



  PseudoBody( const PseudoBody& b )
    : mass( b.mass ), x( b.x ), y( b.y ), z( b.z ),
    xx( b.xx ), yy( b.yy ), zz( b.zz ),
    xy( b.xy ), xz( b.xz ), yz( b.yz ) { }


  inline real getx( void ) const {
    return x;
  }
  inline real gety( void ) const {
    return y;
  }
  inline real getz( void ) const {
    return z;
  }

  inline real getmass( void ) const {
    return mass;
  }

  inline real getxx( void ) const {
    return xx;
  }
  inline real getyy( void ) const {
    return yy;
  }
  inline real getzz( void ) const {
    return zz;
  }
  inline real getxy( void ) const {
    return xy;
  }
  inline real getyx( void ) const {
    return xy;
  }
  inline real getxz( void ) const {
    return xz;
  }
  inline real getzx( void ) const {
    return xz;
  }
  inline real getyz( void ) const {
    return yz;
  }
  inline real getzy( void ) const {
    return yz;
  }
};


class kDTreeNode {
 protected:
  kDTreeNode *m_nleft, *m_nright;

  PseudoBody m_body;
  // real m_lx, m_ly, m_lz, m_l;
  real m_l;
  //real m_xc,m_yc,m_zc;
  bool m_bleaf;
  int m_axis;

  void partition( std::vector<PseudoBody>& All, std::vector<PseudoBody> **bleft, std::vector<PseudoBody> **bright ) {
    unsigned n = All.size();
    *bleft = new std::vector<PseudoBody>;
    *bright = new std::vector<PseudoBody>;

    int splitdim = m_axis;

    if ( splitdim == 0 ) {
      for ( unsigned i = 0; i < n; ++i ) {
	if ( All[ i ].getx() < m_body.getx() )
	  ( *bleft ) ->push_back( All[ i ] );
	else
	  ( *bright ) ->push_back( All[ i ] );
      }
    } else if ( splitdim == 1 ) {
      for ( unsigned i = 0; i < n; ++i ) {
	if ( All[ i ].gety() < m_body.gety() )
	  ( *bleft ) ->push_back( All[ i ] );
	else
	  ( *bright ) ->push_back( All[ i ] );
      }
    } else {
      for ( unsigned i = 0; i < n; ++i ) {
	if ( All[ i ].getz() < m_body.getz() )
	  ( *bleft ) ->push_back( All[ i ] );
	else
	  ( *bright ) ->push_back( All[ i ] );
      }
    }
  }

  void getBoundingBox( const std::vector<PseudoBody>& bodies ) {
    real xmin, xmax, ymin, ymax, zmin, zmax;
    xmin = ymin = zmin = 1.0e+30;
    xmax = ymax = zmax = -1.0e+30;

    real lx( 0.0 ), ly( 0.0 ), lz( 0.0 );
    for ( unsigned i = 0; i < bodies.size(); ++i ) {

      if ( bodies[ i ].getx() < xmin )
	xmin = bodies[ i ].getx();
      if ( bodies[ i ].getx() > xmax )
	xmax = bodies[ i ].getx();

      if ( bodies[ i ].gety() < ymin )
	ymin = bodies[ i ].gety();
      if ( bodies[ i ].gety() > ymax )
	ymax = bodies[ i ].gety();

      if ( bodies[ i ].getz() < zmin )
	zmin = bodies[ i ].getz();
      if ( bodies[ i ].getz() > zmax )
	zmax = bodies[ i ].getz();


    }

    lx = xmax - xmin;
    ly = ymax - ymin;
    lz = zmax - zmin;

    m_l = std::max( std::max( lx, ly ), lz );

  }


 public:

  ~kDTreeNode() {
    if ( m_nleft != NULL )
      delete m_nleft;
    if ( m_nright != NULL )
      delete m_nright;
  }

  kDTreeNode( std::vector<PseudoBody> *Bodies )
    : m_nleft( NULL ),
    m_nright( NULL ),
    m_body( *Bodies ),
    m_bleaf( false ) {

    if ( Bodies->size() > 1 ) {
      getBoundingBox( *Bodies );
    }
  }



  kDTreeNode( std::vector<PseudoBody>& Bodies )
    : m_nleft( NULL ),
    m_nright( NULL ),
    m_body( Bodies ),
    m_bleaf( false ) {

    std::queue< kDTreeNode * > nqueue;
    std::queue< std::vector<PseudoBody>* > nbodyqueue;
    kDTreeNode *nptr;
    std::vector<PseudoBody>* nbodyptr;

    std::vector<PseudoBody> *pbod = new std::vector<PseudoBody>( Bodies );

    nqueue.push( this );
    nbodyqueue.push( pbod );
    getBoundingBox( Bodies );

    m_axis = 0;
    unsigned nodecount( 0 );

    while ( ! nqueue.empty() ) {

      ++nodecount;

      nptr = nqueue.front();
      nqueue.pop();

      nbodyptr = nbodyqueue.front();
      nbodyqueue.pop();

      if ( nbodyptr->size() > 1 ) {

	std::vector<PseudoBody> *LeftBodies, *RightBodies;
	nptr->partition( *nbodyptr, &LeftBodies, &RightBodies );

	nptr->m_nleft = new kDTreeNode( LeftBodies );
	nptr->m_nleft->m_axis = ( nptr->m_axis + 1 ) % 3;
	nqueue.push( nptr->m_nleft );
	nbodyqueue.push( LeftBodies );
	++nodecount;

	nptr->m_nright = new kDTreeNode( RightBodies );
	nptr->m_nright->m_axis = ( nptr->m_axis + 1 ) % 3;
	nqueue.push( nptr->m_nright );
	nbodyqueue.push( RightBodies );
	++nodecount;

      } else {
	nptr->m_nleft = NULL;
	nptr->m_nright = NULL;
	nptr->m_bleaf = true;

      }

      delete nbodyptr;
    }

    // std::cerr << " - kdTree possesses " << nodecount << " nodes, each requiring " << sizeof(kDTreeNode) << " Bytes = " << nodecount*sizeof(kDTreeNode)/1024.0 << "kB\n";
  }


  // do a level order traversal
  real traverse( const real& x, const real& y, const real& z, const real& theta, const real& plumsoft ) {

    std::stack< kDTreeNode * > nqueue;

    nqueue.push( this );
    real totforce = 0.0;
    kDTreeNode *nptr;

    while ( !nqueue.empty() ) {

      nptr = nqueue.top();
      nqueue.pop();

      real
	xx = nptr->m_body.getx() - x,
	yy = nptr->m_body.gety() - y,
	zz = nptr->m_body.getz() - z;

      real r = sqrt( xx * xx + yy * yy + zz * zz );

      if ( nptr->m_bleaf ) {
	if ( r > EPS )
	  totforce += nptr->m_body.getmass() / ( r + plumsoft );

      } else if ( nptr->m_l / r < theta ) {
	double a, b, c, d, A, B;
	SPLINEQ( 1.0 / r, r * r, plumsoft, a, b, c, d );

	A = 0.5 * ( nptr->m_body.getxx() * xx * xx + nptr->m_body.getyy() * yy * yy + nptr->m_body.getzz() * zz * zz
		    + 2.0 * nptr->m_body.getxy() * xx * yy + 2.0 * nptr->m_body.getxz() * xx * zz + 2.0 * nptr->m_body.getyz() * yy * zz );
	B = 0.5 * ( nptr->m_body.getxx() + nptr->m_body.getyy() + nptr->m_body.getzz() );

	totforce -= nptr->m_body.getmass() * a + c * A - b * B;

      } else {
	if ( nptr->m_nleft != NULL )
	  nqueue.push( nptr->m_nleft );
	if ( nptr->m_nright != NULL )
	  nqueue.push( nptr->m_nright );
      }

    }
    return totforce;
  }

  /*real traverse( const real& x, const real& y, const real& z, const real& theta )
    {

    std::queue< kDTreeNode * > nqueue;

    nqueue.push(this);
    real totforce = 0.0;
    kDTreeNode *nptr;

    while( !nqueue.empty() ){

    nptr = nqueue.front();
    nqueue.pop();

    real 
    xx = nptr->m_body.getx()-x,
    yy = nptr->m_body.gety()-y,
    zz = nptr->m_body.getz()-z;

    real r = sqrt(xx*xx+yy*yy+zz*zz);

    if( nptr->m_bleaf ){
    if( r > EPS )
    totforce += nptr->m_body.getmass()/r;
		    
    } else if( nptr->m_l/r < theta ){
    totforce += nptr->m_body.getmass()/r;
    } else{
    if( nptr->m_nleft != NULL )
    nqueue.push( nptr->m_nleft );
    if( nptr->m_nright != NULL )
    nqueue.push( nptr->m_nright );
    }
		   
    }


    return totforce;
    }*/
};






class OctTreeNode {

 protected:
  OctTreeNode **m_subnodes;

  PseudoBody m_body;


  real m_xc, m_yc, m_zc;
  real m_l;
  bool m_bleaf;

  void partition( const std::vector<PseudoBody>& All, std::vector<PseudoBody> ***Bodies ) {
    unsigned n = All.size();

    *Bodies = new std::vector<PseudoBody>*[ 8 ];

    for ( unsigned i = 0; i < 8; ++i )
      ( *Bodies ) [ i ] = new std::vector<PseudoBody>;

    //... assign bodies to octants ...
    for ( unsigned i = 0; i < n; ++i ) {
      unsigned ix = 0, iy = 0, iz = 0;
      if ( All[ i ].getx() - m_xc > 0 )
	ix = 1;
      if ( All[ i ].gety() - m_yc > 0 )
	iy = 1;
      if ( All[ i ].getz() - m_zc > 0 )
	iz = 1;

      ( *Bodies ) [ ix + 2 * ( iy + 2 * iz ) ] ->push_back( All[ i ] );
    }
  }



  void getBoundingBox( const std::vector<PseudoBody>& bodies, real &xc, real &yc, real &zc, real &l ) {
    real xmin, xmax, ymin, ymax, zmin, zmax;
    xmin = ymin = zmin = 1.0e+30;
    xmax = ymax = zmax = -1.0e+30;

    for ( unsigned i = 0; i < bodies.size(); ++i ) {

      if ( bodies[ i ].getx() < xmin )
	xmin = bodies[ i ].getx();
      if ( bodies[ i ].getx() > xmax )
	xmax = bodies[ i ].getx();

      if ( bodies[ i ].gety() < ymin )
	ymin = bodies[ i ].gety();
      if ( bodies[ i ].gety() > ymax )
	ymax = bodies[ i ].gety();

      if ( bodies[ i ].getz() < zmin )
	zmin = bodies[ i ].getz();
      if ( bodies[ i ].getz() > zmax )
	zmax = bodies[ i ].getz();

    }

    l = std::max( std::max( xmax - xmin, ymax - ymin ), zmax - zmin );
    xc = 0.5 * ( xmax + xmin );
    yc = 0.5 * ( ymax + ymin );
    zc = 0.5 * ( zmax + zmin );
  }

 public:

  OctTreeNode( std::vector<PseudoBody> *Bodies )
    : m_subnodes( NULL ),
    m_body( *Bodies ),
    m_bleaf( false ) {
    if ( Bodies->size() > 1 ) {

      getBoundingBox( *Bodies, m_xc, m_yc, m_zc, m_l );
    }
  }

  OctTreeNode( std::vector<PseudoBody>& Bodies, real min_size )  //const std::vector<PseudoBody>& Bodies )
    :
    m_subnodes( NULL ),
    m_body( Bodies ),
    m_xc( 0.0 ),
    m_yc( 0.0 ),
    m_zc( 0.0 ),
    m_l( 0.0 ),
    m_bleaf( false ) {
    std::queue< OctTreeNode * > nqueue;
    std::queue< std::vector<PseudoBody>* > nbodyqueue;
    OctTreeNode *nptr;
    std::vector<PseudoBody>* nbodyptr;

    std::vector<PseudoBody> *pbod = &Bodies; //new std::vector<PseudoBody>( Bodies );

    nqueue.push( this );
    nbodyqueue.push( pbod );
    getBoundingBox( Bodies, m_xc, m_yc, m_zc, m_l );
    //std::cerr << " - building oct-tree from " << Bodies.size() << " particles...\n";
    unsigned nodecount( 0 );
    while ( ! nqueue.empty() ) {
      ++nodecount;
      nptr = nqueue.front();
      nqueue.pop();

      nbodyptr = nbodyqueue.front();
      nbodyqueue.pop();

      if ( nbodyptr->size() > 1 ) {
	std::vector<PseudoBody> **OctantBodies;
	nptr->partition( *nbodyptr, &OctantBodies );
	nptr->m_subnodes = new OctTreeNode * [ 8 ];

	for ( int i = 0; i < 8; ++i ) {
	  if ( OctantBodies[ i ] ->size() > 0 && nptr->m_l > min_size ) {
	    nptr->m_subnodes[ i ] = new OctTreeNode( OctantBodies[ i ] );
	    nqueue.push( nptr->m_subnodes[ i ] );
	    nbodyqueue.push( OctantBodies[ i ] );
	  } else {
	    nptr->m_subnodes[ i ] = NULL;
	    delete OctantBodies[ i ];
	  }

	}

	delete[] OctantBodies;

      } else {
	nptr->m_subnodes = NULL;
	nptr->m_bleaf = true;

      }
      if ( nbodyptr != pbod )
	delete nbodyptr;
    }

    //std::cerr << " - OctTree possesses " << nodecount << " nodes, each requiring " << sizeof(OctTreeNode) << " Bytes = " << nodecount*sizeof(OctTreeNode)/1024.0 << "kB\n";
  }

  ~OctTreeNode() {
    if ( m_subnodes != NULL ) {
      for ( unsigned i = 0; i < 8; ++i )
	if ( m_subnodes[ i ] != NULL ) {
	  delete m_subnodes[ i ];
	}

      delete[] m_subnodes;
    }
  }

  // do a level order traversal
  real traverse( const real& x, const real& y, const real& z, const real& theta, const real& plumsoft ) {

    //std::queue< OctTreeNode * > nqueue;
    std::stack< OctTreeNode * > nqueue;

    nqueue.push( this );
    real totforce( 0.0 );
    OctTreeNode *nptr;

    while ( !nqueue.empty() ) {
      //nptr = nqueue.front();
      nptr = nqueue.top();
      nqueue.pop();

      real
	xx = nptr->m_body.getx() - x,
	yy = nptr->m_body.gety() - y,
	zz = nptr->m_body.getz() - z;

      real r = sqrt( xx * xx + yy * yy + zz * zz );

      if ( nptr->m_bleaf ) {
	if ( r > EPS )
	  totforce -= nptr->m_body.getmass() / ( r + plumsoft );

      } else if ( nptr->m_l / r < theta ) {

	double a, b, c, d, A, B;
	SPLINEQ( 1.0 / r, r * r, plumsoft, a, b, c, d );

	A = 0.5 * ( nptr->m_body.getxx() * xx * xx + nptr->m_body.getyy() * yy * yy + nptr->m_body.getzz() * zz * zz
		    + 2.0 * nptr->m_body.getxy() * xx * yy + 2.0 * nptr->m_body.getxz() * xx * zz + 2.0 * nptr->m_body.getyz() * yy * zz );
	B = 0.5 * ( nptr->m_body.getxx() + nptr->m_body.getyy() + nptr->m_body.getzz() );

	totforce -= nptr->m_body.getmass() * a + c * A - b * B;

      } else
	for ( unsigned i = 0; i < 8; ++i ) {
	  if ( nptr->m_subnodes[ i ] != NULL )
	    nqueue.push( nptr->m_subnodes[ i ] );
	}
    }
    return totforce;
  }
};

#endif
