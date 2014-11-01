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

/**
 * \file ScalarField.cpp
 * \brief Phi
 */

#include "ScalarField.h"
#include "Stencil.h"

ScalarField::ScalarField( void ){
  std::cout << " in ScalarField()..." << std::endl;
  nx_ = 0;
  ny_ = 0;
  nz_ = 0;
  phi_ = new Array3(512,512,512);
  dx_ = 0.0f;
  nStencil_ = 0;
  nStencil2_ = 0;
}

ScalarField::ScalarField( unsigned nx, unsigned ny, unsigned nz, 
			  float dx, int nStencil ){
  std::cout << " in ScalarField(" << nx << "," << ny << "," << nz;
  std::cout << "," << dx << "," << nStencil << ")..."<< std::endl;

  nx_ = nx;
  ny_ = ny;
  nz_ = nz;
  phi_ = new Array3( nx, ny, nz );   //dy_( local_nx, ny, nz ),
  //phi_ = new float(nx*ny*nz);
  dx_ = dx;
  nStencil_ = nStencil;
  nStencil2_ = ( nStencil_ - 1 ) / 2;
}

ScalarField::~ScalarField(){
  //  phi_.~Array3();
}

ScalarField*
ScalarField::operator=(ScalarField* s){
  nx_ = s->getNx();
  ny_ = s->getNy();
  nz_ = s->getNz();
  phi_ = s->getPhi();
  dx_ = s->getDx();
  nStencil_ = s->getNStencil();
  nStencil2_ = s->getNStencil2();
  return this;
}

//float
//ScalarField::operator() ( int ix, int iy, int iz ) {
//  return phi_.getValue(ix,iy,iz);
//}

float
ScalarField::getValue( int ix, int iy, int iz ) {
  return phi_.getValue(ix,iy,iz);
}

void
ScalarField::setValue(int ix, int iy, int iz, float val){
  phi_.setValue(ix,iy,iz,val);
}

void
ScalarField::recreate( unsigned nx, unsigned ny, unsigned nz, 
		       float dx, int nStencil ) {
  nx_ = nx;
  ny_ = ny;
  nz_ = nz;
  dx_ = dx;
  nStencil_ = nStencil;
  nStencil2_ = ( nStencil - 1 ) / 2;
  phi_.recreate( nx_, ny_, nz_ );
}

Array3*
ScalarField::getPhi( void ){
  return &phi_;
}

void
ScalarField::setPhi( Array3* phi ){
  phi_ = phi;
}

int
ScalarField::getNx( void ){
  return nx_;
}

int
ScalarField::getNy( void ){
  return ny_;
}

int
ScalarField::getNz( void ){
  return nz_;
}

float
ScalarField::getDx( void ){
  return dx_;
}

int
ScalarField::getNStencil( void ){
  return nStencil_;
}

int
ScalarField::getNStencil2(void){
  return nStencil2_;
}

void 
ScalarField::lookup( float *Coord, unsigned nCoord, 
		     std::vector< float >& values ) {
  values.clear();
  for ( unsigned i = 0; i < nCoord; ++i ) {
    values.push_back( this->getValue( Coord[3*i],Coord[3*i+1],Coord[3*i+2] ) );
  }
}

void
ScalarField::lookup( float *Coord, unsigned nCoord, float ll, unsigned N, 
		     std::vector< float >& values ) {
  values.clear();
  float factor = 1.0; //(float)N/ll;
  for ( unsigned i = 0; i < nCoord; ++i ) {
    values.push_back( this->getValue( factor*Coord[ 3*i ], 
				      factor*Coord[ 3*i+1 ], 
				      factor*Coord[ 3*i+2 ] ) );
  }
}


float 
ScalarField::get_CIC( float coordx, float coordy, float coordz ) {
  
  int
    ix = ( int ) coordx,
    iy = ( int ) coordy,
    iz = ( int ) coordz,
    ix1 = ( ix + 1 ) % nx_,
    iy1 = ( iy + 1 ) % ny_,
    iz1 = ( iz + 1 ) % nz_;
  
  
  float
    dx = ( float ) ( coordx - ( float ) ix ),
    dy = ( float ) ( coordy - ( float ) iy ),
    dz = ( float ) ( coordz - ( float ) iz ),
    tx = 1.0 - dx,
    ty = 1.0 - dy,
    tz = 1.0 - dz;

  float
    f_xyz = ( *this ).getValue( ix, iy, iz ) * tx * ty * tz,
    f_Xyz = ( *this ).getValue( ix1, iy, iz ) * dx * ty * tz,
    f_xYz = ( *this ).getValue( ix, iy1, iz ) * tx * dy * tz,
    f_xyZ = ( *this ).getValue( ix, iy, iz1 ) * tx * ty * dz,
    f_XYz = ( *this ).getValue( ix1, iy1, iz ) * dx * dy * tz,
    f_XyZ = ( *this ).getValue( ix1, iy, iz1 ) * dx * ty * dz,
    f_xYZ = ( *this ).getValue( ix, iy1, iz1 ) * tx * dy * dz,
    f_XYZ = ( *this ).getValue( ix1, iy1, iz1 ) * dx * dy * dz;
  
  return f_xyz + f_Xyz + f_xYz + f_xyZ + f_XYz + f_XyZ + f_xYZ + f_XYZ;
}

float
ScalarField::diff1( int ix, int iy, int iz, int i ) {
  int iright[ 3 ];
  for ( int k = 0; k < 3; ++k ){
    iright[ k ] = 0;
  }
  iright[ i ] = + 1;
  
  float d1i = 0.0;
  
  for ( int ii = -nStencil2_; ii <= nStencil2_; ++ii ) {
    int iix = ( ix + ii * iright[ 0 ] + nx_ ) % nx_;
    int iiy = ( iy + ii * iright[ 1 ] + ny_ ) % ny_;
    int iiz = ( iz + ii * iright[ 2 ] + nz_ ) % nz_;
    d1i += pStencil1_[ ii + nStencil2_ ] * this->getValue( iix, iiy, iiz );
    }
  
  return d1i;
}


float
ScalarField::diff2( int ix, int iy, int iz, int i, int j ) {
  // trace component of the Hessian?
  if ( i == j ) {
    
    int iright[ 3 ];
    
    for ( int k = 0; k < 3; ++k )
      iright[ k ] = 0;
    iright[ i ] = + 1;
    
    float d2ii = 0.0;
    
    for ( int ii = -nStencil2_; ii <= nStencil2_; ++ii ) {
      int iix = ( ix + ii * iright[ 0 ] + nx_ ) % nx_;
      int iiy = ( iy + ii * iright[ 1 ] + ny_ ) % ny_;
      int iiz = ( iz + ii * iright[ 2 ] + nz_ ) % nz_;
      d2ii += pStencil2_[ ii + nStencil2_ ] * this->getValue( iix, iiy, iiz );
    }
    return d2ii;
  }
  
  // off-trace component of the Hessian
  int jright[ 3 ];
  
  for ( int k = 0; k < 3; ++k )
    jright[ k ] = 0;
  jright[ j ] = + 1;
  
  float d1j = 0;
  
  for ( int ii = -nStencil2_; ii <= nStencil2_; ++ii ) {
    int iix = ( ix + ii * jright[ 0 ] + nx_ ) % nx_;
    int iiy = ( iy + ii * jright[ 1 ] + ny_ ) % ny_;
    int iiz = ( iz + ii * jright[ 2 ] + nz_ ) % nz_;
    d1j += pStencil1_[ ii + nStencil2_ ] 
      * this->diff1( iix, iiy, iiz, i );
  }
  
  return d1j;
}

void
ScalarField::test(){
  std::cout << " in ScalarField::test()..." << std::endl;

  for(unsigned ix=0; ix<nx_; ++ix){
    for(unsigned iy=0; iy<ny_; ++iy){
      for(unsigned iz=0; iz<nz_; ++iz){
	phi_.setValue( ix,iy,iz,-1.0f);
      }
    }
  }
}

/**
 * @brief Put clound in cell
 * @param npart number of particles
 * @param coord positions of particles
 * @param wpar mass ("weight") of particles
 */
void 
ScalarField::put_CIC( unsigned npart, float *coord, float *wpar ) {

  unsigned ix, iy, iz, ix1, iy1, iz1;
  float x, y, z;
  float dx, dy, dz, dyw;
  float tx, ty, tz, tyw;

  static bool bInitialized( false );

  if ( !bInitialized ) {
    bInitialized = true;
    std::cout << " * * initializing density field with -1...";
    for ( ix = 0; ix < nx_; ++ix ){
      for ( iy = 0; iy < ny_; ++iy ){
	for ( iz = 0; iz < nz_; ++iz ){
	  phi_.setValue( ix, iy, iz, -1.0f);
	}
      }
    }
    std::cout << "done!" << std::endl;
  }

  for ( unsigned i = 0; i < npart; ++i ) {
    x = coord[ 3*i+0 ];
    y = coord[ 3*i+1 ];
    z = coord[ 3*i+2 ];

    ix = unsigned ( x );
    iy = unsigned ( y );
    iz = unsigned ( z );

    dx = x - ( float ( ix ) );
    dy = y - ( float ( iy ) );
    dz = z - ( float ( iz ) );

    ix %= nx_;
    iy %= ny_;
    iz %= nz_;

    tx = 1.0 - dx;
    ty = 1.0 - dy;
    tz = 1.0 - dz;

    tyw = ty * wpar[ i ];
    dyw = dy * wpar[ i ];

    ix1 = ( ix + 1 ) % nx_;
    iy1 = ( iy + 1 ) % ny_;
    iz1 = ( iz + 1 ) % nz_;

    phi_.addValue( ix, iy, iz, tz * tx * tyw);
    phi_.addValue( ix, iy1, iz, tz * tx * dyw);
    phi_.addValue( ix, iy, iz1, dz * tx * tyw);
    phi_.addValue( ix, iy1, iz1, dz * tx * dyw);

    phi_.addValue( ix1, iy, iz, tz * dx * tyw);
    phi_.addValue( ix1, iy1, iz, tz * dx * dyw);
    phi_.addValue( ix1, iy, iz1, dz * dx * tyw);
    phi_.addValue( ix1, iy1, iz1, dz * dx * dyw);

  }
}


/**
 * \brief diff2atCIC
 * @param Coord: list of coordinates
 * @param nPoints: number of points in Coord where we want to know the tensor
 * @param t_ij: tensor
 * @return void
 */
void
ScalarField::diff2atCIC( float *Coord, unsigned nPoints,
			 std::vector< Matrix33 > &t_ij ) {

  ChooseStencil();

  for ( unsigned ip = 0; ip < nPoints; ++ip ) {
    int ix = int ( Coord[ 3 * ip + 0 ] );
    int iy = int ( Coord[ 3 * ip + 1 ] );
    int iz = int ( Coord[ 3 * ip + 2 ] );

    if ( ix > int ( nx_ ) ) {
      std::cerr << " * * ERROR: index out of bounds! ";
      std::cerr << Coord[ 3 * ip ] << std::endl;
      ix = nx_;
    }

    int
      ix1 = ( ix + 1 ) % nx_,
      iy1 = ( iy + 1 ) % ny_,
      iz1 = ( iz + 1 ) % nz_;

    float x = ( ( float ) Coord[ 3*ip + 0 ] );
    float y = ( ( float ) Coord[ 3*ip + 1 ] );
    float z = ( ( float ) Coord[ 3*ip + 2 ] );

    float
      cic_dx = x - ( float ) ix;
    float cic_dy = y - ( float ) iy;
    float cic_dz = z - ( float ) iz;
    float cic_tx = 1.0 - cic_dx;
    float cic_ty = 1.0 - cic_dy;
    float cic_tz = 1.0 - cic_dz;

    float xyz = cic_tx * cic_ty * cic_tz;
    float Xyz = cic_dx * cic_ty * cic_tz;
    float xYz = cic_tx * cic_dy * cic_tz;
    float XYz = cic_dx * cic_dy * cic_tz;
    float xyZ = cic_tx * cic_ty * cic_dz;
    float XyZ = cic_dx * cic_ty * cic_dz;
    float xYZ = cic_tx * cic_dy * cic_dz;
    float XYZ = cic_dx * cic_dy * cic_dz;

    Matrix33 t;

    for ( int i = 0; i < 3; ++i ){
      for ( int j = 0; j < 3; ++j ) {

	float
	  f_ijk, f_Ijk, f_iJk, f_IJk,
	  f_ijK, f_IjK, f_iJK, f_IJK;

	float
	  dy_ijk, dy_Ijk, dy_iJk, dy_IJk,
	  dy_ijK, dy_IjK, dy_iJK, dy_IJK;

	f_ijk = this->diff2( ix, iy, iz, i, j );
	f_Ijk = this->diff2( ix1, iy, iz, i, j );
	f_iJk = this->diff2( ix, iy1, iz, i, j );
	f_IJk = this->diff2( ix1, iy1, iz, i, j );
	f_ijK = this->diff2( ix, iy, iz1, i, j );
	f_IjK = this->diff2( ix1, iy, iz1, i, j );
	f_iJK = this->diff2( ix, iy1, iz1, i, j );
	f_IJK = this->diff2( ix1, iy1, iz1, i, j );

	dy_ijk = f_ijk * xyz;
	dy_Ijk = f_Ijk * Xyz;
	dy_iJk = f_iJk * xYz;
	dy_IJk = f_IJk * XYz;
	dy_ijK = f_ijK * xyZ;
	dy_IjK = f_IjK * XyZ;
	dy_iJK = f_iJK * xYZ;
	dy_IJK = f_IJK * XYZ;

	t( i, j ) = dy_ijk + dy_Ijk + dy_iJk + dy_IJk + dy_ijK + dy_IjK + dy_iJK + dy_IJK;
      }
    }
    t_ij.push_back( t );
  }
}



/**
 * @brief Choose a stencil.
 * Choose a stencil
 * @return void
 */
void
ScalarField::ChooseStencil( void ) {
  switch ( nStencil_ ) {

  case 3:
    pStencil1_ = stencil3_1;
    pStencil2_ = stencil3_2;
    break;

  case 5:
    pStencil1_ = stencil5_1;
    pStencil2_ = stencil5_2;
    break;

  case 7:
    pStencil1_ = stencil7_1;
    pStencil2_ = stencil7_2;
    break;

  case 9:
    pStencil1_ = stencil9_1;
    pStencil2_ = stencil9_2;
    break;

  case 11:
    pStencil1_ = stencil11_1;
    pStencil2_ = stencil11_2;
    break;

  default:
    std::cerr << " * * ERROR: unsupported finite difference scheme ";
    std::cerr << "in ScalarField::ChooseStencil!" << std::endl;
    abort();
  }
}
