#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern struct
{
   int ii, jj, kk;
} ijk_;

int printsample_(int* ndim){
  printf("function call by name worked!\n");
  printf(" - with parameter %d\n", *ndim);
  return(0);
}


int doubleijk_(char *cc, int ll)
{
   cc[ll--] = '\0';  // NULL terminate the string
   printf("Inside doubleIJK: %s\n",cc);
   ijk_.ii *=2;
   ijk_.jj *=2;
   ijk_.kk *=2;
   return(0);
}

/****** loglikelihood routine *************************/
// Now an example, sample an egg box likelihood
// Input arguments
// ndim 						= dimensionality (total number of free parameters) of the problem
// npars 						= total number of free plus derived parameters
// context						void pointer, any additional information
//
// Input/Output arguments
// Cube[npars] 						= on entry has the ndim parameters in unit-hypercube
//	 						on exit, the physical parameters plus copy any derived parameters you want to store with the free parameters
//
// Output arguments
// lnew 						= loglikelihood
void* loglike_(double *Cube, int *ndim, int *npars, double *lnew, void *context)
{
	double chi = 1.0;
	int i;
	for(i = 0; i < *ndim; i++)
	{
		double x = Cube[i]*10.0*M_PI;
		chi *= cos(x/2.0);
		Cube[i] = x;
	}
	*lnew = powf(chi + 2.0, 5.0); // value changed in here
}


/******** dumper routine ********************************/
// The dumper routine will be called every updInt*10 iterations
// MultiNest does not need the user to do anything. Use the arguments any way you want.
// nSamples 						= total number of samples in posterior distribution
// nlive 						= total number of live points
// nPar 						= total number of parameters (free + derived)
// physLive[1][nlive * (nPar + 1)] 			= 2D array containing the last set of live points (physical parameters plus derived parameters) along with their loglikelihood values
// posterior[1][nSamples * (nPar + 2)] 			= posterior distribution containing nSamples points. Each sample has nPar parameters (physical + derived) along with the their loglike value & posterior probability
// paramConstr[1][4*nPar]:
// paramConstr[0][0] to paramConstr[0][nPar - 1] 	= mean values of the parameters
// paramConstr[0][nPar] to paramConstr[0][2*nPar - 1] 	= standard deviation of the parameters
// paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1] = best-fit (maxlike) parameters
// paramConstr[0][nPar*4] to paramConstr[0][4*nPar - 1] = MAP (maximum-a-posteriori) parameters
// maxLogLike						= maximum loglikelihood value
// logZ							= log evidence value from the default (non-INS) mode
// INSlogZ						= log evidence value from the INS mode
// logZerr						= error on log evidence value
// context						void pointer, any additional information
void dumper_(int *nSamples, int *nlive, int *nPar, double **physLive,
             double **posterior, double **paramConstr, double *maxLogLike,
             double *logZ, double *INSlogZ, double *logZerr, void *context)
{
	// convert the 2D Fortran arrays to C arrays
	// the posterior distribution
	// postdist will have nPar parameters in the first nPar columns
    // & loglike value & the posterior probability in the last two columns
	int i, j;
	double postdist[*nSamples][*nPar + 2];
	for( i = 0; i < *nPar + 2; i++ )
		for( j = 0; j < *nSamples; j++ )
			postdist[j][i] = posterior[0][i * (*nSamples) + j];

	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns
    // & loglike value in the last column
	double pLivePts[*nlive][*nPar + 1];
	for( i = 0; i < *nPar + 1; i++ )
		for( j = 0; j < *nlive; j++ )
			pLivePts[j][i] = physLive[0][i * (*nlive) + j];
}
