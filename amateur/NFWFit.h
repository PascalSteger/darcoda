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
 * \file NFWFit.h
 * \brief Fit of a Navarro-Frenk-White potential
 */

#ifndef NFWFIT_H
#define NFWFIT_H

#include <cassert>
#include <iostream>
#include <vector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_errno.h>

#include "Global.h"

struct fitdata {

  const size_t n;
  const real *x;
  const real *y;
  const real *sigma;

  fitdata( const size_t _n, const real *_x, const real *_y, const real *_sigma )
    : n( _n ), x( _x ), y( _y ), sigma( _sigma ) {}
}
;

struct fitdata_fix {
  const size_t n;
  const real *x;
  const real *y;
  const real *sigma;
  const real rvir;
  const real mvir;

  fitdata_fix( const size_t _n, const real *_x, const real *_y,
	       const real *_sigma, const real _rvir, const real _mvir )
    : n( _n ), x( _x ), y( _y ), sigma( _sigma ), rvir( _rvir ), mvir( _mvir ) {}
}
;


void 
print_state ( size_t iter, gsl_multifit_fdfsolver * s ) {
  /*std::cerr << "iter: " <<iter << "x= " <<gsl_vector_get (s->x, 0) //<< " " <<gsl_vector_get (s->x, 1)
    << "|f(x)|=" << gsl_blas_dnrm2 (s->f) << "\n";*/
}

int 
NFW_f( const gsl_vector *p, void* data, gsl_vector *f ) {
  struct fitdata * dataptr = ( struct fitdata* ) data;
  size_t n( dataptr->n );
  const real *x( dataptr->x );
  const real *y( dataptr->y );
  const real *sigma( dataptr->sigma );

  real deltac = gsl_vector_get( p, 0 );
  real rs = gsl_vector_get( p, 1 );

  for ( size_t i = 0; i < n; ++i ) {
    real r = x[ i ];
    real Yi = deltac / ( ( r / rs ) * ( 1.0 + r / rs ) * ( 1.0 + r / rs ) );
    gsl_vector_set( f, i, ( Yi - y[ i ] ) / sigma[ i ] );
  }

  return GSL_SUCCESS;
}

int 
NFW_fix_f( const gsl_vector *p, void* data, gsl_vector *f ) {
  struct fitdata_fix * dataptr = ( struct fitdata_fix* ) data;
  size_t n( dataptr->n );
  const real *x( dataptr->x );
  const real *y( dataptr->y );
  const real *sigma( dataptr->sigma );
  const real rvir( dataptr->rvir );
  const real mvir( dataptr->mvir );

  //  real deltac = gsl_vector_get(p,0);
  real rs = gsl_vector_get( p, 0 );

  real fac = mvir / ( 4.0 * M_PI * rs * rs * rs ) / ( log( ( rvir + rs ) / rs ) - rvir / ( rvir + rs ) );

  for ( size_t i = 0; i < n; ++i ) {
    real r = x[ i ];
    real Yi = fac / ( ( r / rs ) * ( 1.0 + r / rs ) * ( 1.0 + r / rs ) );
    gsl_vector_set( f, i, ( Yi - y[ i ] ) / sigma[ i ] );
  }

  return GSL_SUCCESS;
}

int 
NFW_df( const gsl_vector *p, void *data, gsl_matrix *J ) {
  struct fitdata * dataptr = ( struct fitdata* ) data;
  size_t n( dataptr->n );
  const real *x( dataptr->x );
  const real *sigma( dataptr->sigma );

  real deltac = gsl_vector_get( p, 0 );
  real rs = gsl_vector_get( p, 1 );

  for ( size_t i = 0; i < n; ++i ) {
    real x1 = x[ i ] / rs;
    real x1p1 = 1.0 + x[ i ] / rs;

    real df_dp1 = 1.0 / ( x1 * x1p1 * x1p1 );
    real df_dp2 = deltac / ( x[ i ] * x1p1 * x1p1 ) + 2.0 * deltac / ( rs * x1p1 * x1p1 * x1p1 );

    gsl_matrix_set( J, i, 0, df_dp1 / sigma[ i ] );
    gsl_matrix_set( J, i, 1, df_dp2 / sigma[ i ] );
  }

  return GSL_SUCCESS;
}

/*int 
  NFW_fix_df( const gsl_vector *p, void *data, gsl_matrix *J )
  {
  struct fitdata_fix* dataptr = (struct fitdata_fix*)data;
  size_t n( dataptr->n );
  const real *x( dataptr->x );
  const real *sigma( dataptr->sigma );
  const real rvir( dataptr->rvir );
  const real mvir( dataptr->mvir );
 
  real rs     = gsl_vector_get(p,0);
  real fac1 = mvir/(4.0*M_PI*rs)/pow(log((rvir+rs)/rs)*(rvir+rs)-rvir,2);
  real fac2 = log((rvir+rs)/rs);
 
  for( size_t i=0; i<n; ++i )
  {
  real r    = x[i];///rs;
  //      real x1p1 = 1.0+x[i]/rs;
 
  real df_dp1 = fac1/(r*(rs+r)*(rs+r)*(rs+r))*(fac2*(-2.0*rvir*rvir*rs-4.0*rvir*rs*rs-2.0*rs*rs*rs)+3.0*rvir*rs+rvir*rvir*r+2*rs*rs*rvir);;
  //      real df_dp2 = deltac/(x[i]*x1p1*x1p1)+2.0*deltac/(rs*x1p1*x1p1*x1p1);
 
  gsl_matrix_set(J,i,0,df_dp1/sigma[i]);
  //      gsl_matrix_set(J,i,1,df_dp2/sigma[i]);
  }
 
  return GSL_SUCCESS;
  }*/


int 
NFW_fix_df( const gsl_vector *p, void *data, gsl_matrix *J ) {
  struct fitdata_fix * dataptr = ( struct fitdata_fix* ) data;
  size_t n( dataptr->n );
  const real *x( dataptr->x );
  const real *sigma( dataptr->sigma );
  const real rvir( dataptr->rvir );
  const real mvir( dataptr->mvir );

  real rs = gsl_vector_get( p, 0 );

  real a1 = mvir / ( 4.0 * M_PI * rs ) / pow( log( ( rvir + rs ) / rs ) * ( rvir + rs ) - rvir, 2 );
  real a2 = log( ( rvir + rs ) / rs ) * ( -2.0 * rvir * rvir * rs - 4.0 * rvir * rs * rs - 2.0 * rs * rs * rs ) + 3.0 * rvir * rvir * rs + 2.0 * rvir * rs * rs;

  for ( size_t i = 0; i < n; ++i ) {
    real r = x[ i ];
    real df_dp1 = a1 * ( a2 + rvir * rvir * r ) / ( r * ( rs + r ) * ( rs + r ) * ( rs + r ) );
    gsl_matrix_set( J, i, 0, df_dp1 / sigma[ i ] );
  }

  return GSL_SUCCESS;
}

int 
NFW_fdf( const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J ) {
  NFW_f( x, data, f );
  NFW_df( x, data, J );
  return GSL_SUCCESS;
}

int 
NFW_fix_fdf( const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J ) {
  NFW_fix_f( x, data, f );
  NFW_fix_df( x, data, J );
  return GSL_SUCCESS;
}

int 
NFW_fit( const std::vector<real>& r, const std::vector<real>& rho,
	 const std::vector<real>& sigma_rho,
	 real deltac0, real rs0,
	 real& deltac, real& edeltac,
	 real& rs, real& ers ) {
  assert( ( r.size() == rho.size() ) && ( r.size() == sigma_rho.size() ) );

  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  unsigned iter( 0 );
  const size_t n = r.size();
  const size_t np = 2;

  gsl_matrix *covar = gsl_matrix_alloc( np, np );

  //  struct fitdata d = { n, &r[0], &rho[0], &sigma_rho[0] };
  struct fitdata d( n, &r[ 0 ], &rho[ 0 ], &sigma_rho[ 0 ] );

  gsl_multifit_function_fdf f;
  double x_init[ 2 ] = { deltac0, rs0 };

  gsl_vector_view x = gsl_vector_view_array( x_init, np );

  f.f = &NFW_f;
  f.df = &NFW_df;
  f.fdf = &NFW_fdf;
  f.n = n;
  f.p = np;
  f.params = &d;

  gsl_set_error_handler_off();

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc( T, n, np );
  gsl_multifit_fdfsolver_set( s, &f, &x.vector );

  print_state( iter, s );

  int status;
  do {
    ++iter;
    status = gsl_multifit_fdfsolver_iterate( s );
    print_state( iter, s );

    if ( status )
      break;

    status = gsl_multifit_test_delta( s->dx, s->x, 1e-4, 1e-4 );
  } while ( status == GSL_CONTINUE && iter < 500 );

  gsl_multifit_covar( s->J, 0.0, covar );

  //#define FIT(i) gsl_vector_get(s->x, i )
  //#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  real redchisq;
  real chi = gsl_blas_dnrm2( s->f );
  real dof = n - np;
  real c = GSL_MAX_DBL( 1, chi / sqrt( dof ) );

  redchisq = chi * chi / dof;

  gsl_multifit_fdfsolver_free( s );
  gsl_matrix_free( covar );

  if ( iter < 500 ) {
    deltac = gsl_vector_get( s->x, 0 );
    rs = gsl_vector_get( s->x, 1 );

    edeltac = c * sqrt( gsl_matrix_get( covar, 0, 0 ) );
    ers = c * sqrt( gsl_matrix_get( covar, 1, 1 ) );

    return 1;
  }
  return 0;

}


int 
NFW_fix_fit( const std::vector<real>& r, const std::vector<real>& rho,
	     const std::vector<real>& sigma_rho,
	     real rvir, real mvir,
	     real rs0, real& rs, real& ers ) {
  assert( ( r.size() == rho.size() ) && ( r.size() == sigma_rho.size() ) );

  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  unsigned iter( 0 );
  const size_t n = r.size();
  const size_t np = 1;

  gsl_matrix *covar = gsl_matrix_alloc( np, np );

  //struct fitdata_fix d = { n, &r[0], &rho[0], &sigma_rho[0], rvir, mvir  };
  struct fitdata_fix d( n, &r[ 0 ], &rho[ 0 ], &sigma_rho[ 0 ], rvir, mvir );

  gsl_multifit_function_fdf f;
  double x_init[ np ] = { rs0 };

  gsl_vector_view x = gsl_vector_view_array( x_init, np );

  f.f = &NFW_fix_f;
  f.df = &NFW_fix_df;
  f.fdf = &NFW_fix_fdf;
  f.n = n;
  f.p = np;
  f.params = &d;

  gsl_set_error_handler_off();

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc( T, n, np );
  gsl_multifit_fdfsolver_set( s, &f, &x.vector );

  print_state( iter, s );

  int status;
  do {
    ++iter;
    status = gsl_multifit_fdfsolver_iterate( s );
    print_state( iter, s );

    if ( status )
      break;

    status = gsl_multifit_test_delta( s->dx, s->x, 1e-4, 1e-4 );
  } while ( status == GSL_CONTINUE && iter < 500 );

  gsl_multifit_covar( s->J, 0.0, covar );

  //#define FIT(i) gsl_vector_get(s->x, i )
  //#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  real redchisq;
  real chi = gsl_blas_dnrm2( s->f );
  real dof = n - np;
  real c = GSL_MAX_DBL( 1, chi / sqrt( dof ) );

  redchisq = chi * chi / dof;

  gsl_multifit_fdfsolver_free( s );
  gsl_matrix_free( covar );

  if ( iter < 500 ) {
    rs = gsl_vector_get( s->x, 0 );
    //rs      = gsl_vector_get(s->x,1);

    //    edeltac = c*sqrt(gsl_matrix_get(covar,0,0));
    ers = c * sqrt( gsl_matrix_get( covar, 0, 0 ) );

    return 1;
  }

  return 0;

}

/*
 
template< class Model >
class GSLMultiFit{
protected:
Model m_Model;
 
const gsl_multifit_fdfsolver_type *m_pType;
gsl_multifit_fdfsolver *m_pSolver;
 
public:
GSLMultiFit( void )
: m_Model()
{}
  
  
void fit( const std::vector<real>& x,
const std::vector<real>& y,
const std::vector<real>& guess,
std::vector<real>& fitparam )
{
  std::assert( x.size() == y.size() );
 
  const size_t np = x.size();
  const size_t p  = m_Model.getnpar();
 
    
    
};
*/

#endif
