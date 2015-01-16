//copyright 2008 Oliver Hahn, 2010 Pascal Steger

/**
 * \file TidalField.cpp
 * \brief Implements -t option: Computation of the tidal field.
 */

//#define WITH_MPI
//#include "mpi++.h"

#include <unistd.h>
#include <sstream>

#include "TidalField.h"
#include "HDFIO.h"
#include "PoissonSolver.h"
#include "Eigensystem3.h"
#include "Global.h"


//in kpc, not Mpc
const float G_const = 43011.7902;


void TidalField( unsigned isnap, int N, double Radius ) {
  std::cout << " in TidalField()..." << std::endl;
  //TBD: change to 11 for different smoothing scales
  for(unsigned i=1; i<2; ++i){
    std::string foldsnap = Folder("output") + Snaps( isnap );
    std::string snaphdf5 = foldsnap + "/snapshot.hdf5";
    std::string extra    = foldsnap +"/Hprop.txt";
    std::string tensorhdf5 = foldsnap +"/tensor_"+IntToString(i)+".hdf5";
    std::cout << snaphdf5 << std::endl;
    std::cout << extra << std::endl;
    std::cout << tensorhdf5 << std::endl;
    TidalField( snaphdf5, N, i*Radius, extra, tensorhdf5 );
  }
}

/**
 * @brief Computes the tidal field from a given mass distribution.
 * @param snaphdf5: a snapshot file, converted by AnalyzeSnap
 * @param N: the length of the density mesh
 * @param Radius: the length scale used for smoothing
 * @param extra: file containing the -z output of the halos 
 *               (used for npart, COM)
 * @param tensorhdf5: output file for the tensor components
 * @return void
 */
void TidalField( std::string snaphdf5, unsigned N, double Radius, 
		 std::string extra, std::string tensorhdf5 ) {
  std::cout << "TidalField v1.2 - "
	    << "by O. Hahn, P. Steger (psteger@phys.ethz.ch)" << std::endl;
  std::cout << "--------------------------------------------"
	    << "------------------" << std::endl << std::endl;

  unsigned nParticles = 0;
  float boxsize = 0.0f;
  std::vector<float> pmass;
  std::vector<float> ppos;
  write(" * Reading particle properties...");
  // attention: what type of particles?
  HDFReadGroupAttribute( snaphdf5, "Header", "NumParticles", nParticles );
  HDFReadGroupAttribute( snaphdf5, "Header", "Boxsize", boxsize );
  HDFReadDataset( snaphdf5, "Data_Mass", pmass );
  HDFReadVector( snaphdf5, "Data_Pos", ppos );
  writedone();
  //for(unsigned i=0; i<nParticles; ++i)
  //	std::cout<<pmass.at(i)<<std::endl;

  float dx = boxsize / N;
  float Omegam0 = 0.25f;
  float a = 1.0f;
  float H0 = 100.0;
  float rho0 = 3.0f * H0 * H0 / 8.0 / M_PI / G_const * Omegam0;
  // what's the difference to dx?
  float r0 = boxsize / N;
  float t0 = 1.0f / H0;
  float v0 = r0 / t0;
  float phi0 = v0 * v0;

  // filter size is directly specified, Gaussian
  float sigma = Radius / r0;

  std::cout << " * # particles     = " << nParticles << std::endl;
  std::cout << " * boxsize         = " << boxsize << "kpc/h" << std::endl;
  std::cout << " * smooth @ Radius = " << Radius << "kpc/h" << std::endl;
  std::cout << " * sigma           = " << sigma << std::endl;
  std::cout << " * gridlength N    = " << N << std::endl;
  std::cout << " * r0              = " << r0 << std::endl;
  std::cout << " * dx              = " << dx << std::endl;

  unsigned nx=N;
  unsigned ny=N;
  unsigned nz=N;

  const int nstencil = 3;

  for ( unsigned i = 0; i < ppos.size(); ++i ){
    ppos[i] /= r0;
  }

  std::cout << " * Putting Cloud in Cell..." << std::endl;
  ScalarField rho( nx, ny, nz, dx, nstencil );
  rho.put_CIC( nParticles, &ppos[ 0 ], &pmass[ 0 ] );
  std::cout << " * done!" << std::endl;

  ScalarField phi( nx, ny, nz, dx, nstencil );
  PoissonSolver<nstencil> psolver( &rho, &phi );
  std::cout << " * Invoking FFT Poisson solver..." << std::endl;
  psolver.solve( boxsize, sigma );
  writedone();

  // file input:
  std::ifstream file(extra.c_str());
  write( " * Reading number of Halos..." );
  int idummy;
  float fdummy;
  unsigned nhalos;
  file >> idummy >> fdummy >> nhalos;
  std::cout << " gives " << nhalos << ", done!" << std::endl;

  std::cout << " * Reading positions of COM..." << std::flush;
  std::vector<float> xcm;//( size 3 * nhalos at the end);
  for(unsigned i=0; i<nhalos; ++i){
    file >>idummy>>idummy>>fdummy>>fdummy>>fdummy>>fdummy;
    for(int j=0; j<3; ++j){
      file >> fdummy;
      //      cout << fdummy << " ";
      xcm.push_back(fdummy);
    }
    //cout << endl;
    //{int a; cin >> a;}
    for(int j=0; j<9; ++j){
      file>>fdummy;
    }
    for(int j=0; j<5; ++j){
      file>>idummy;
      for(int k=0; k<18; ++k){
	file>>fdummy;
      }
    }
  }
  file.close();
  writelndone();
  for( unsigned i = 0; i<nhalos; ++i){
    cout << xcm.at(i) << endl;
  }
  {int a; cin >> a;}

  write( " * Converting coordinates to grid units...");
  for ( unsigned i = 0; i < 3*nhalos; ++i )
    xcm.at( i ) /= r0;
  writelndone();
  for( unsigned i = 0; i<nhalos; ++i){
    cout << xcm.at(i) << endl;
  }
  {int a; cin >> a;}

  HDFCreateFile( tensorhdf5 );
  write( " * Computing tensor components at given positions...");
  std::vector< Matrix33 > tij_cm;
  // dphi.diff_at( tij_cm, xcm, nhalos, boxsize );
  phi.diff2atCIC( &xcm[ 0 ], nhalos, tij_cm );
  writelndone();

  write( " * Converting back to physical units...");
  float tnorm = 1.0 / ( 1.5 * Omegam0 / a );
  for ( unsigned i = 0; i < nhalos; ++i ){
    for ( int j = 0; j < 3; ++j ){
      for ( int k = 0; k < 3; ++k ){
	( tij_cm[ i ] ) ( j, k ) *= tnorm;
      }
    }
  }
  writelndone();

  StoreTensorField( tensorhdf5, nhalos, xcm, tij_cm, sigma );

  write( " * Evaluating density field at given positions..." );
  std::vector< float > datadens;

  for ( unsigned i = 0; i < nhalos; ++i ){
    datadens.push_back( 
		     ( rho.get_CIC( xcm.at( 3*i+0 ),
				    xcm.at( 3*i+1 ),
				    xcm.at( 3*i+2 ) ) + 1.0 ) 
		     * rho0 / ( a * a * a ));
  }

  HDFWriteDataset( tensorhdf5, "Density", datadens );
  datadens.clear();

  std::vector<float> dataov;
  for ( unsigned i = 0; i < nhalos; ++i ){
    dataov.push_back( 
		     rho.get_CIC( xcm.at( 3*i+0 ),
				  xcm.at( 3*i+1 ),
				  xcm.at( 3*i+2 ) ));
  }

  HDFWriteDataset( tensorhdf5, "Overdensity", dataov );
  dataov.clear();

  std::vector<float> datapot;
  for ( unsigned i = 0; i < nhalos; ++i ){
    datapot.push_back( 
		      rho.get_CIC( xcm.at( 3*i+0 ),
				   xcm.at( 3*i+1 ),
				   xcm.at( 3*i+2 ) ) * phi0);
  }

  HDFWriteDataset( tensorhdf5, "Potential", datapot );
  datapot.clear();
}




/**
 * @brief Store tensor field.
 * store tensor field to HDF5
 * @param filename: path to output file
 * @param nPoints: number of points
 * @param coordinates: evaluation points
 * @param aTensorField: tensor field that should be evaluated
 * @param SmoothingMassScale: smoothing scale
 */
void StoreTensorField( std::string filename, unsigned nPoints,
		       const std::vector<float>& coordinates,
		       const std::vector< Matrix33 > &aTensorField,
		       float SmoothingMassScale ) {
  std::cout << " in StoreTensorField..." << std::endl;
  std::vector< std::vector<float> > tensor;

  std::cout << " * * Writing data for " << nPoints ;
  std::cout << " points to file " << filename.c_str();
  std::cout << "..." << std::endl;
  std::vector<float> eigenvalues,eigenvector1,eigenvector2,eigenvector3;
  for ( unsigned i = 0; i < nPoints; ++i ) {
    // tensor and trace-free tensor
    std::vector<float> tensortemp, evaltemp, evecttemp;
    for ( unsigned k = 0; k < 3; ++k ){
      for ( unsigned l = 0; l <= k; ++l ){
	tensortemp.push_back( aTensorField[i] ( l, k ) );
      }
    }

    //std::cout << " running eigensystem calculation..." << std::endl;
    Matrix33 A;
    A = aTensorField[ i ];

    Eigensystem3 E(A);
    //std::vector<float> lambda;
    //std::vector<Vector3> V;
    //A.Eigen( lambda, V );

    //std::cout << " extracting eigenvectors..." << std::endl;
    for ( unsigned j = 0; j < 3; ++j ) {
      eigenvalues.push_back( E.getEval(j) );
      eigenvector1.push_back( (E.getEvec(0))(j) );
      eigenvector2.push_back( (E.getEvec(1))(j) );
      eigenvector3.push_back( (E.getEvec(2))(j) );
    }

    tensor.push_back( tensortemp );
  }

  std::cout << " writing coordinates and eigensystem for tidal field..."; 
  std::cout << std::endl;
  HDFWriteDatasetVector( filename, "Coordinates", coordinates );
  HDFWriteDataset2D( filename, "TensorField", tensor );
  HDFWriteDatasetVector( filename, "Eigenvector_1", eigenvector1 );
  HDFWriteDatasetVector( filename, "Eigenvector_2", eigenvector2 );
  HDFWriteDatasetVector( filename, "Eigenvector_3", eigenvector3 );
  HDFWriteDatasetVector( filename, "Eigenvalues", eigenvalues );
}
