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
 * \file PoissonSolver_new.h
 * \brief Functions to solve the Poisson problem for rho
 * The _new suffix means that the new fftw3.h is used
 */

#ifndef POISSONSOLVER_H
#define POISSONSOLVER_H

#include <iostream>
#include "fftw3/fftw3.h"
#include "TidalField.h"

#define POSdouble(ix,iy,iz) (((ix)*(ny) + (iy)) * (2*((nz)/2+1)) + (iz)) 
#define ACC_R(ix,iy,iz) (((ix)*(ny) + (iy)) * ((nz) + (iz)))
#define ACC_C(ix,iy,iz) (((ix)*(ny) + (iy)) * (((nz)/2+1) + (iz)))

//#define ACC_R(ix,iy,iz) (((ix)*(ny) + (iy)) * (nz) + (iz) )
//#define ACC_C(ix,iy,iz) (((ix)*(ny) + (iy)) * ((nz)/2+1) + (iz))

#define c_re(c) ((c)[0])
#define c_im(c) ((c)[1])


//int POSdouble(int ix, int iy, int iz ){
//  return (ix*ny + iy) * (2*(nz/2+1)) + iz;

#ifdef GAUSSIAN_KERNEL
#warning "FYI: Enabling Gaussian Kernel in PoissonSolver"
#else
#warning "FYI: Enabling Top Hat Kernel in PoissonSolver"
#endif

template < int nstencil >
class PoissonSolver {
 protected:
  ScalarField *m_prho;

 public:
  PoissonSolver( ScalarField *prho )
    : m_prho( prho ) { }

    inline double Gauss3D( const double &k, const double &sigma ) {
      double v = 1.0 / ( pow( sigma * sigma * 2.0 * M_PI, 1.5 ) );
      v *= exp( -0.5 * k * k / ( sigma * sigma ) );
      return v;
    }

    inline double TopHat3D( const double &k, const double &R ) {
      if ( k <= R ) {
	double norm = 4.0 / 3.0 * M_PI * R * R * R;
	return 1.0 / norm;
      }
      return 0.0;
    }



    void solve( ScalarField *pphi, double Boxlength, double sigma ) {
      int nx = m_prho->get_xsize();
      int ny = m_prho->get_ysize();
      int nz = m_prho->get_zsize(); 

      GreensFunction<nstencil> gf( 2.0 * M_PI / ( double ) m_prho->get_ysize() );


      double *data = new double[ nx * ny * nz ];
      double *kern = new double[ nx * ny * nz ];


      // copy data without ghost cells to array data for FFTW
      for ( int ix = 0; ix < nx; ++ix )
	for ( int iy = 0; iy < ny; ++iy )
	  for ( int iz = 0; iz < nz; ++iz ) {

	    int iix( ix ), iiy( iy ), iiz( iz );

	    if ( iix > nx / 2 )
	      iix -= nx;
	    if ( iiy > ny / 2 )
	      iiy -= ny;
	    if ( iiz > nz / 2 )
	      iiz -= nz;

	    double k = sqrt( iix * iix + iiy * iiy + iiz * iiz );



#ifdef TOPHAT_KERNEL
	    kern[ ACC_R( ix, iy, iz ) ] = TopHat3D( k, sigma );
#else
	    kern[ ACC_R( ix, iy, iz ) ] = Gauss3D( k, sigma );
						
#endif
	    data[ ACC_R( ix, iy, iz ) ] = ( double ) ( *m_prho ) ( ix, iy, iz );
	  }

      fftw_complex *cdata = new fftw_complex[ nx * ny * ( nz / 2 + 1 ) ];
      //fftw_plan plan  = rfftw3d_create_plan( nx, ny, nz,
      //				 FFTW_REAL_TO_COMPLEX,
      //				 FFTW_ESTIMATE ); 
      //rfftwnd_one_real_to_complex( plan, &data[ 0 ], &cdata[ 0 ] );
      fftw_plan plan = fftw_plan_dft_r2c_3d( nx, ny, nz,
					     &data[0], &cdata[0],
					     FFTW_ESTIMATE);
      fftw_execute( plan );
      fftw_destroy_plan( plan ); // FFTW3 uses one method to destroy all different plans!
      std::cout << c_re(cdata[0]) << std::endl;
    
      fftw_complex *ckern = new fftw_complex[ nx * ny * ( nz / 2 + 1 ) ];
      //fftw_plan kplan = rfftw3d_create_plan( nx, ny, nz,
      //				 FFTW_REAL_TO_COMPLEX,
      //				 FFTW_ESTIMATE );
      //rfftwnd_one_real_to_complex( kplan, &kern[ 0 ], &ckern[ 0 ] );
      fftw_plan kplan = fftw_plan_dft_r2c_3d( nx, ny, nz,
					      &kern[0], &ckern[0],
					      FFTW_ESTIMATE);
      fftw_execute( kplan );
      fftw_destroy_plan( kplan );    

      double Omegam0 = 0.25;
      double a = 1.0;

      double fftnorm = 1.0 / ( nx * ny * nz );
      double mGpi = -1.5 * Omegam0 / a * fftnorm;

      for ( int ix = 0; ix < nx; ++ix )
	for ( int iy = 0; iy < ny; ++iy )
	  for ( int iz = 0; iz < ( nz / 2 + 1 ); ++iz ) {

	    int iix( ix ), iiy( iy ), iiz( iz );

	    if ( iix > nx / 2 )
	      iix -= nx;
	    if ( iiy > ny / 2 )
	      iiy -= ny;

	    double fac = mGpi / ( gf.apply( iix ) + gf.apply( iiy ) + gf.apply( iiz ) );

	    double dre = c_re(cdata[ ACC_C( ix, iy, iz ) ]); 
	    double dim = c_im(cdata[ ACC_C( ix, iy, iz ) ]);


#if 0
	    //c_re(ckern[ijk]) = dre*fftnorm;
	    //c_im(ckern[ijk]) = dim*fftnorm;

	    //c_re(cdata[ijk]) *= fac;
	    //c_im(cdata[ijk]) *= fac;

	    c_re(ckern[ ix ][ iy ][ iz ]) = dre * fftnorm;
	    c_im(ckern[ ix ][ iy ][ iz ]) = dim * fftnorm;

	    //std::cerr << c_re(cdata[ix][iy][iz]) << ", ";

	    c_re(cdata[ ix ][ iy ][ iz ]) *= fac;
	    c_im(cdata[ ix ][ iy ][ iz ]) *= fac;

#else
	    double kre = c_re(ckern[ ACC_C( ix, iy, iz ) ]);
	    double kim = c_im(ckern[ ACC_C( ix, iy, iz ) ]); 

	    c_re(ckern[ ACC_C( ix, iy, iz ) ]) = ( dre * kre - dim * kim ) * fftnorm;
	    c_im(ckern[ ACC_C( ix, iy, iz ) ]) = ( dre * kim + dim * kre ) * fftnorm;

	    c_re(cdata[ ACC_C( ix, iy, iz ) ]) = ( dre * kre - dim * kim ) * fac;
	    c_im(cdata[ ACC_C( ix, iy, iz ) ]) = ( dre * kim + dim * kre ) * fac;

#endif

	  }

      c_re(cdata[ 0 ]) = 0.0;
      c_im(cdata[ 0 ]) = 0.0;

      //fftw_plan iplan = rfftw3d_create_plan( nx, ny, nz,
      //				 FFTW_COMPLEX_TO_REAL,
      //				 FFTW_ESTIMATE ); 
      //rfftwnd_one_complex_to_real( iplan, &cdata[ 0 ], &data[ 0 ] );
      fftw_plan iplan = fftw_plan_dft_c2r_3d( nx, ny, nz,
					      &cdata[0], &data[0],
					      FFTW_ESTIMATE );
      fftw_execute( iplan );
      fftw_destroy_plan( iplan ); // FFTW3 uses one method to destroy all different plans!
      delete[] cdata;

      //fftw_plan ikplan = rfftw3d_create_plan( nx, ny, nz,
      //				 FFTW_COMPLEX_TO_REAL,
      //				 FFTW_ESTIMATE );
      //rfftwnd_one_complex_to_real( ikplan, &ckern[ 0 ], &kern[ 0 ] );
      fftw_plan ikplan = fftw_plan_dft_c2r_3d( nx, ny, nz,
					       &ckern[0], &kern[0],
					       FFTW_ESTIMATE );
      fftw_execute( ikplan );
      fftw_destroy_plan( ikplan );
      delete[] ckern;

      for ( int ix = 0; ix < nx; ++ix )
	for ( int iy = 0; iy < ny; ++iy )
	  for ( int iz = 0; iz < nz; ++iz )
	    ( *m_prho ) ( ix, iy, iz ) = kern[ ACC_R( ix, iy, iz ) ]; //kern[ix][iy][iz];

      delete[] kern;

      // recreate pphi to ensure correct size
      pphi->recreate( nx, ny, nz, m_prho->get_dx(), m_prho->get_stencil() );

      // copy data to obey ghost cell structure
      for ( int ix = 0; ix < nx; ++ix )
	for ( int iy = 0; iy < ny; ++iy )
	  for ( int iz = 0; iz < nz; ++iz ) {
	    ( *pphi ) ( ix, iy, iz ) = data[ ACC_R( ix, iy, iz ) ]; //data[ix][iy][iz];//data[(ix*ny + iy) * (2*(nz/2+1)) + iz];
	  }

      delete[] data;

    }

};


#endif
