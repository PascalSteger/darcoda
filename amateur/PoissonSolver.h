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

/*!
 * \file PoissonSolver.h
 * \brief Functions to solve the Poisson problem for rho
 */

#ifndef POISSONSOLVER_H
#define POISSONSOLVER_H

//#include "sfftw.h" // old inclusions for fftw 2.*
//#include "srfftw.h"
#include "fftw3.h"
#include "TidalField.h"
#include "ScalarField.h"

#define fftw_real double
#define POSdouble(ix,iy,iz) (((ix)*(ny) + (iy)) * (2*((nz)/2+1)) + (iz)) 
#define ACC_R(ix,iy,iz) (( (ix) * (ny) + (iy) ) * (nz) + (iz) )
#define ACC_C(ix,iy,iz) (( (ix) * (ny) + (iy) ) * ((nz)/2+1) + (iz))

#ifdef GAUSSIAN_KERNEL
// #warning "FYI: Enabling Gaussian Kernel in PoissonSolver"
#else
// #warning "FYI: Enabling Top Hat Kernel in PoissonSolver"
#endif

/**
 * implementation of Poisson equation solver using Fourier transformation
 * of the cloud-in-cell density field,
 * templated for the number of stencils
 */
template < int nstencil >
class PoissonSolver {
 protected:
  ScalarField* prho_;
  ScalarField* pphi_;

 public:
  PoissonSolver( ScalarField *prho, ScalarField *pphi )
    { prho_ = prho; pphi_ = pphi;}

  /**
   * \brief Gaussian kernel
   */
  inline double Gauss3D( const double &k, const double &sigma ) {
    double v = 1.0 / ( pow( sigma * sigma * 2.0 * M_PI, 1.5 ) );
    v *= exp( -0.5 * k * k / ( sigma * sigma ) );
    return v;
  }

  /**
   * top hat kernel: all points inside radius R are fully accounted for,
   *                 all points outside are not
   */
  inline double TopHat3D( const double &k, const double &R ) {
    if ( k <= R ) {
      double norm = (4.0*M_PI/3.0) * R * R * R;
      return 1.0 / norm;
    }
    return 0.0;
  }

  /**
   * little test routine for access to prho_ and pphi_
   */
  void test(){
    prho_->setValue(0,0,0,0.3f);
    cout << prho_->getValue(0,0,0) << endl;
    cout << pphi_->getValue(0,0,0) << endl;
    for(int i=0; i<512; ++i){
      pphi_->setValue(i,i,i,-i);
    }
    cout << pphi_->getValue(0,0,0) << endl;
    
  }

  /**
   * invoke solver
   */
  // old syntax
  // void solve( ScalarField *pphi, double Boxlength, double sigma ) {
  void solve( double Boxlength, double sigma ){
    cout << " in solve()..." << endl;
    int nx = prho_->getNx();
    int ny = prho_->getNy();
    int nz = prho_->getNz();

    GreensFunction<nstencil> gf( 2.0 * M_PI / ( double ) nx );

    cout << " allocating memory..." << endl;
    fftw_real *data = new fftw_real[ nx * ny * nz ];
    fftw_real *kern = new fftw_real[ nx * ny * nz ];

    cout << " copying data..." << endl;
    // copy data without ghost cells to array data for FFTW
    for ( int ix = 0; ix < nx; ++ix ) {
      int iix = (ix>=nx/2?ix:ix-nx/2);
      for ( int iy = 0; iy < ny; ++iy ) {
	int iiy = (iy>=ny/2?iy:iy-ny/2);
	for ( int iz = 0; iz < nz; ++iz ) {
	  int iiz = (iz>=nz/2?iz:iz-nz/2);
	  double k = sqrt( iix*iix + iiy*iiy + iiz*iiz );

#ifdef TOPHAT_KERNEL
	  kern[ ACC_R( ix, iy, iz ) ] = TopHat3D( k, sigma );
#else
	  kern[ ACC_R( ix, iy, iz ) ] = Gauss3D( k, sigma );
#endif
	  data[ ACC_R( ix, iy, iz ) ] = (double) prho_->getValue( ix, iy, iz );
	}
      }
    }

    cout << " creating cdata..." << endl;
    fftw_complex *cdata = new fftw_complex[ nx * ny * ( nz / 2 + 1 ) ];
    cout << " creating plan..." << endl;
    fftw_plan plan   = fftw_plan_dft_r2c_3d( nx, ny, nz, &data[0], &cdata[0],
					     FFTW_ESTIMATE ); 
    cout << " executing plan..." << endl;
    fftw_execute( plan );
    cout << " destroying plan..." << endl;
    fftw_destroy_plan(plan);

    cout << " creating ckern..." << endl;
    fftw_complex *ckern = new fftw_complex[ nx * ny * ( nz / 2 + 1 ) ];
    cout << " creating kplan..." << endl;
    fftw_plan kplan  = fftw_plan_dft_r2c_3d( nx, ny, nz, &kern[0], &ckern[0],
					     FFTW_ESTIMATE );
    cout << " executing kplan..." << endl;
    fftw_execute( kplan );
    cout << " destroying kplan..." << endl;
    fftw_destroy_plan(kplan);

    double Omegam0 = 0.25;
    double a = 1.0;
    double fftnorm = 1.0 / ( nx * ny * nz );
    double mGpi = -1.5 * Omegam0 / a * fftnorm;

    cout << " copying data..." << endl;
    for ( int ix = 0; ix < nx; ++ix ) {
      int iix = (ix<=nx/2?ix:ix-nx);
      for ( int iy = 0; iy < ny; ++iy ) {
	int iiy = (iy<=ny/2?iy:iy-ny);
	for ( int iz = 0; iz < ( nz / 2 + 1 ); ++iz ) {

	  double fac = mGpi/(gf.apply( iix )+gf.apply( iiy )+gf.apply( iz ) );
	  double dre = cdata[ ACC_C( ix, iy, iz ) ][0];//real part 
	  double dim = cdata[ ACC_C( ix, iy, iz ) ][1];//imag part
#if 0
	  //ckern[ijk][0] = dre*fftnorm;
	  //ckern[ijk][1] = dim*fftnorm;
	  //cdata[ijk][0] *= fac;
	  //cdata[ijk][1] *= fac;

	  ckern[ ACC_C( ix, iy, iz ) ][0] = dre * fftnorm;//real
	  ckern[ ACC_C( ix, iy, iz ) ][1] = dim * fftnorm;//imag
	  //std::cerr << cdata[ACC_C(ix,iy,iz)][0] << ", ";

	  cdata[ ACC_C( ix, iy, iz ) ][0] *= fac;
	  cdata[ ACC_C( ix, iy, iz ) ][1] *= fac;
#else
	  double kre = ckern[ ACC_C( ix, iy, iz ) ][0];//real part
	  double kim = ckern[ ACC_C( ix, iy, iz ) ][1];//imag part

	  ckern[ ACC_C( ix, iy, iz ) ][0] = (dre*kre - dim*kim) * fftnorm;
	  ckern[ ACC_C( ix, iy, iz ) ][1] = (dre*kim + dim*kre) * fftnorm;

	  cdata[ ACC_C( ix, iy, iz ) ][0] = (dre*kre - dim*kim) * fac;
	  cdata[ ACC_C( ix, iy, iz ) ][1] = (dre*kim + dim*kre) * fac;
#endif
	}
      }
    }

    cdata[ 0 ][0] = 0.0;
    cdata[ 0 ][1] = 0.0;

    cout << " creating iplan..." << endl;
    fftw_plan iplan  = fftw_plan_dft_c2r_3d( nx, ny, nz, &cdata[0], &data[0],
					     FFTW_ESTIMATE ); 
    cout << " executing iplan..." << endl;
    fftw_execute( iplan );
    cout << " destroying iplan..." << endl;
    fftw_destroy_plan(iplan);

    cout << " creating ikplan..." << endl;
    fftw_plan ikplan = fftw_plan_dft_c2r_3d( nx, ny, nz, &ckern[0], &kern[0],
					     FFTW_ESTIMATE ); 
    cout << " executing ikplan..." << endl;
    fftw_execute( ikplan );
    cout << " destroying ikplan..." << endl;
    fftw_destroy_plan(ikplan);
    cout << " deleting ckern..." << endl;
    delete[] ckern;

    cout << " setting prho_..." << endl;
    for ( int ix = 0; ix < nx; ++ix ){
      for ( int iy = 0; iy < ny; ++iy ){
	for ( int iz = 0; iz < nz; ++iz ){
	  prho_->setValue( ix, iy, iz, kern[ ACC_R( ix, iy, iz ) ]); 
	  //kern[ix][iy][iz];
	}
      }
    }
    
    cout << " deleting kern..." << endl;
    delete[] kern;
    

#if 0

    if ( true ) { //MPI::COMM_WORLD.Get_rank() == 0 ){
      std::ofstream ofs( "dump_rho.dat" );
      for ( int ix = 0; ix < nx; ++ix ) {
	for ( int iy = 0; iy < ny; ++iy ) {
	  for ( int iz = 0; iz < nz; ++iz ){
	    ofs << std::setw( 16 ) << ( *prho_ ) ( ix, iy, iz ) << " ";
	  //ofs << std::setw(16) << data[(ix*ny + iy) * (2*(nz/2+1)) + iz];
	  }
	  ofs << std::endl;
	}
	ofs << std::endl;
      }
      ofs.close();
    }
#endif

#if 0
    if ( true ) { //MPI::COMM_WORLD.Get_rank() == 0 ){
      std::ofstream ofs( "dump_phi.dat" );
      for ( int ix = 0; ix < nx; ++ix ) {
	for ( int iy = 0; iy < ny; ++iy ) {
	  for ( int iz = 0; iz < nz; ++iz ) {
	    //ofs << std::setw(16) << pphi_->getValue(ix,iy,iz) << " ";
	    ofs << std::setw( 16 ) << data[ ACC_R( ix, iy, iz ) ] << " ";
	    //data[(ix*ny + iy) * (2*(nz/2+1)) + iz] << " ";
	  }
	  ofs << std::endl;
	}
	ofs << std::endl;
      }
      ofs.close();
    }
#endif
    
#if 0
    if ( true ) { //MPI::COMM_WORLD.Get_rank() == 0 ){
      std::ofstream ofs( "dump_DelPhiDirect.dat" );
      for ( int ix = 0; ix < nx; ++ix ) {
	for ( int iy = 0; iy < ny; ++iy ) {
	  for ( int iz = 0; iz < nz; ++iz ) {
	    //ofs << std::setw(16) << pphi->getValue(ix,iy,iz) << " ";
	    //double f0   = data[POSdouble(ix,iy,iz)];
	    //double fxm2 = data[POSdouble((ix+nx-1)%nx,iy,iz)];
	    //double fxp2 = data[POSdouble((ix+1)%nx,iy,iz)];
	    //double fym2 = data[POSdouble(ix,(iy+ny-1)%ny,iz)];
	    //double fyp2 = data[POSdouble(ix,(iy+1)%ny,iz)];
	    //double fzm2 = data[POSdouble(ix,iy,(iz+nz-1)%nz )];
	    //double fzp2 = data[POSdouble(ix,iy,(iz+1)%nz)];
	    
	    double f0 = data[ ACC_R( ix, iy, iz ) ];
	    double fxm2 = data[ ACC_R( (ix + nx - 1 ) % nx, iy, iz) ];
	    double fxp2 = data[ ACC_R( (ix + 1 ) % nx, iy, iz) ];
	    double fym2 = data[ ACC_R( ix, ( iy + ny - 1 ) % ny, iz )];
	    double fyp2 = data[ ACC_R( ix, ( iy + 1 ) % ny,iz) ];
	    double fzm2 = data[ ACC_R( ix, iy, ( iz + nz - 1 ) % nz) ];
	    double fzp2 = data[ ACC_R( ix, iy, ( iz + 1 ) % nz )];
	    
	    double D2f = 0.25 * ( fxm2 + fxp2 + fym2 + fyp2 + fzm2 + fzp2
				  - 6.0 * f0 );
	    
	    ofs << std::setw( 16 ) << ( D2f ) * 10.66667 - 1 << " ";
	  }
	  ofs << std::endl;
	}
	ofs << std::endl;
      }
      ofs.close();
    }
#endif
    
    // recreate pphi to ensure correct size
    //pphi_->recreate( nx, ny, nz, prho_->getDx(), prho_->getNStencil() );
    
    // copy data to obey ghost cell structure
    cout << " copying to pphi_..." << endl;
    for ( int ix = 0; ix < nx; ++ix ) {
      for ( int iy = 0; iy < ny; ++iy ) {
    	for ( int iz = 0; iz < nz; ++iz ) {
	  pphi_->setValue( ix, iy, iz, data[ ACC_R( ix, iy, iz ) ]); 
	  //data[ix][iy][iz];//data[(ix*ny + iy) * (2*(nz/2+1)) + iz];
	}
      }
    }

    
    delete[] data;
    delete[] cdata;
  }
  
};


#endif
