

#ifndef SRHDPACK_HEADER_INCLUDED
#define SRHDPACK_HEADER_INCLUDED
#include "cow.h"

#define SRHDPACK_VELOCITY_GAMMA          -42
#define SRHDPACK_VELOCITY_BETA           -43
#define SRHDPACK_VELOCITY_GAMMABETA      -44
#define SRHDPACK_VELOCITY_DUMUDXMU       -45
#define SRHDPACK_VELOCITY_DUMUDUMU       -46
#define SRHDPACK_SEPARATION_LAB          -47
#define SRHDPACK_SEPARATION_PROPER       -48
#define SRHDPACK_PROJECTION_NONE         -49
#define SRHDPACK_PROJECTION_TRANSVERSE   -50
#define SRHDPACK_PROJECTION_LONGITUDINAL -51

typedef struct srhdpack_samplemode
{
  double exponent; // exponent value: p
  int velmode; // SRHDPACK_VELOCITY
  int sepmode; // SRHDPACK_SEPARATION
  int projmode; // SRHDPACK_PROJECTION
} srhdpack_samplemode;


char *srhdpack_onepointpdfs(cow_dfield *prim,
			    char *which,
			    char *h5fname,
			    char *h5gname,
			    double *interval);

void srhdpack_shelevequescaling(cow_dfield *vel,
				cow_histogram *hist,
				int velmode,
				int sepmode,
				int projmode,
				int nbatch,
				int nperbatch,
				int seed,
				double exponent);

void srhdpack_collectpairs(cow_dfield *vel,
			   srhdpack_samplemode *modes,
			   int num_modes,
			   int num_pairs,
			   int num_samps,
			   double *samploc,
			   double *outbufx,
			   double *outbufy);

#endif // SRHDPACK_HEADER_INCLUDED
