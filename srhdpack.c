
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cow.h"
#include "srhdpack.h"
#define MODULE "srhdpack"

static void sl94(cow_dfield *vel, cow_histogram *hist,
		 int velmode,
		 int sepmode,
		 int projmode,
		 int exponent,
		 int N);
static void boost(double u[4], double x[4], double xp[4]);
static double len3(double *x);
static double gamm(double *x);
static double dot3(double *u, double *v);
static void project3(double *dx, double *u, double *ulat, double *ulon);

void srhdpack_shelevequescaling(cow_dfield *vel,
				cow_histogram *hist,
				int velmode,
				int sepmode,
				int projmode,
				int exponent,
				int nbatch,
				int nperbatch,
				int seed)
// -----------------------------------------------------------------------------
// This function computes the pairwise relative Lorentz factor between many
// points in the 3-velocity field `vel`. The user needs to supply two
// half-initialized histograms, which have not yet been committed. This function
// will commit, populate, and seal the histograms before returning them. The
// supplies the fields like in the example below, all other will be
// over-written.
//
//  cow_histogram_setnbins(hist, 0, 256);
//  cow_histogram_setspacing(hist, COW_HIST_SPACING_LINEAR); // or LOG
//  cow_histogram_setnickname(hist, "myhist"); // optional
//  cow_histogram_setlower(histpro, 0, 0.1); // good values for most cases
//  cow_histogram_setupper(histpro, 0, 1.5);
//
// The histogram `histpro` bins the pairs in terms of their proper space-like
// separation, while `histlab` uses their lab-frame separations.
//
// NOTE: this function generates its own sample points by calling the math.h
// rand() function. If running in parallel, the user is responsible for making
// sure to seed the ranks differently, unless `seed` is set to 1, in which case
// srand will be called with the domain's cartesian communicator rank.
//
// Modes:
//
// Velocity options: GAMMA, BETA, GAMMABETA
// Separation vector: LAB, PROPER
// Projection options: LONGITUDINAL, TRANSVERSE
// Exponent options: p=1,2,...
// -----------------------------------------------------------------------------
{
  cow_domain *d = cow_dfield_getdomain(vel);
  if (seed) {
    printf("[%s] seed=%d, seeding the random number generator\n", MODULE, seed);
    srand(cow_domain_getcartrank(d));
  }
  else {
    printf("[%s] seed=%d\n", MODULE, seed);
  }
  cow_histogram_setdomaincomm(hist, d);
  cow_histogram_setbinmode(hist, COW_HIST_BINMODE_AVERAGE);
  cow_histogram_commit(hist);

  for (int n=0; n<nbatch; ++n) {
    sl94(vel, hist, velmode, sepmode, projmode, exponent, nperbatch);
    printf("[%s] running batch %d/%d of size %d\n", MODULE, n, nbatch,
	   nperbatch);
  }
  cow_histogram_seal(hist);
}


void sl94(cow_dfield *vel, cow_histogram *hist,
	  int velmode,
	  int sepmode,
	  int projmode,
	  int exponent,
	  int N)
{
  int npair = N;
  int nsamp = N*2;

  double *x = (double*) malloc(nsamp * 3 * sizeof(double));
  double *v;

  for (int n=0; n<nsamp; ++n) {
    x[3*n + 0] = (double) rand() / RAND_MAX;
    x[3*n + 1] = (double) rand() / RAND_MAX;
    x[3*n + 2] = (double) rand() / RAND_MAX;
  }

  cow_dfield_setsamplemode(vel, COW_SAMPLE_LINEAR);
  int err = cow_dfield_setsamplecoords(vel, x, nsamp, 3);
  if (err) {
    printf("[%s] Error! setsamplecoords returned %d\n", MODULE, err);
  }
  free(x);

  int nout1, nout2;
  cow_dfield_sampleexecute(vel);
  cow_dfield_getsamplecoords(vel, &x, &nout1, NULL);
  cow_dfield_getsampleresult(vel, &v, &nout2, NULL);

  double x1[3];
  double x2[3];

  for (int n=0; n<npair; ++n) {
    int i1 = rand() % nsamp;
    int i2 = rand() % nsamp;
    memcpy(x1, &x[3*i1], 3*sizeof(double));
    memcpy(x2, &x[3*i2], 3*sizeof(double));
    // -------------------------------------------------------------------------
    // With periodic BC's on the unit cube, points are never actually more than
    // 1/2 away from one another along a given axis. Of the two possible
    // orderings of points x1 and x2 along axis `d`, we choose the one which
    // keeps them closer together.
    // -------------------------------------------------------------------------
    for (int d=0; d<3; ++d) {
      if (x1[d] - x2[d] > 0.5) {
	x1[d] -= 1.0;
      }
      else if (x1[d] - x2[d] < -0.5) {
	x1[d] += 1.0;
      }
    }
    double *v1 = &v[3*i1];
    double *v2 = &v[3*i2];
    double g1 = gamm(v1);
    double g2 = gamm(v2);
    double umu1[4] = { g1, g1*v1[0], g1*v1[1], g1*v1[2] };
    double umu2[4] = { g2, g2*v2[0], g2*v2[1], g2*v2[2] };
    double dxlab[4] = { 0.0, x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2] };
    double dvlab[4] = { 0.0, v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2] };
    double dxpro[4];
    double dupro[4];

    boost(umu1, umu2, dupro);
    boost(umu1, dxlab, dxpro);

    double drlab = len3(&dxlab[1]);
    double drpro = len3(&dxpro[1]);
    double dvlat[4];
    double dvlon[4];
    double dulat[4];
    double dulon[4];

    project3(dxlab, &dvlab[1], &dvlat[1], &dvlon[1]);
    project3(dxlab, &dupro[1], &dulat[1], &dulon[1]);

    dvlat[0] = len3(&dvlat[1]);
    dvlon[0] = len3(&dvlon[1]);
    dulat[0] = gamm(&dulat[1]);
    dulon[0] = gamm(&dulon[1]);

    double xvalue;
    double yvalue;
    double gammarel;
    double betarel;
    double gammabetarel;

    switch (sepmode) {
    case SRHDPACK_SEPARATION_PROPER:
      xvalue = drpro;
      break;
    case SRHDPACK_SEPARATION_LAB:
      xvalue = drlab;
      break;
    default:
      printf("[%s] Error! invalid argument: sepmode\n", MODULE);
      return;
    }

    switch (projmode) {
    case SRHDPACK_PROJECTION_NONE:
      gammarel = dupro[0];
      betarel = len3(&dvlab[1]);
      gammabetarel = len3(&dupro[1]);
      break;
    case SRHDPACK_PROJECTION_TRANSVERSE:
      gammarel = dulat[0];
      betarel = len3(&dvlat[1]);
      gammabetarel = len3(&dulat[1]);
      break;
    case SRHDPACK_PROJECTION_LONGITUDINAL:
      gammarel = dulon[0];
      betarel = len3(&dvlon[1]);
      gammabetarel = len3(&dulon[1]);
      break;
    default:
      printf("[%s] Error! invalid argument: projmode\n", MODULE);
      return;
    }

    switch (velmode) {
    case SRHDPACK_VELOCITY_GAMMA:
      yvalue = gammarel;
      break;
    case SRHDPACK_VELOCITY_BETA:
      yvalue = betarel;
      break;
    case SRHDPACK_VELOCITY_GAMMABETA:
      yvalue = gammabetarel;
      break;
    default:
      printf("[%s] Error! invalid argument: velmode\n", MODULE);
      return;
    }
    cow_histogram_addsample1(hist, xvalue, pow(yvalue, exponent));
  }
}

double len3(double *x)
{
  return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}
double dot3(double *u, double *v)
{
  return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}
double gamm(double *x)
{
  return 1.0 / sqrt(1.0 - len3(x));
}
void project3(double *dx, double *u, double *ulat, double *ulon)
// -----------------------------------------------------------------------------
// Returns the projection of the vector u onto the vector dx in ulat, and its
// projection onto the plane normal to dx in ulon.
// -----------------------------------------------------------------------------
{
  double dxlen = len3(dx); if (dxlen<1e-12) dxlen = 1.0; // don't divide by zero
  double dxhat[3] = { dx[0] / dxlen, dx[1] / dxlen, dx[2] / dxlen };
  double ulonlen = dot3(dxhat, ulon);
  ulon[0] = ulonlen * dxhat[0];
  ulon[1] = ulonlen * dxhat[1];
  ulon[2] = ulonlen * dxhat[2];
  ulat[0] = u[0] - ulon[0];
  ulat[1] = u[1] - ulon[1];
  ulat[2] = u[2] - ulon[2];
}
void boost(double u[4], double x[4], double xp[4])
// -----------------------------------------------------------------------------
// Maps the 4-vector x into xp by the boost formed from the 4-velocity u.
//
// see:
// http://en.wikipedia.org/wiki/Lorentz_transformation#Boost_in_any_direction
// -----------------------------------------------------------------------------
{
  double u2 = u[1]*u[1] + u[2]*u[2] + u[3]*u[3];
  double gm = sqrt(1.0 + u2);
  double L[4][4];
  L[0][0] = gm;
  L[0][1] = L[1][0] = -u[1];
  L[0][2] = L[2][0] = -u[2];
  L[0][3] = L[3][0] = -u[3];
  for (int i=1; i<4; ++i) {
    for (int j=1; j<4; ++j) {
      L[i][j] = (gm - 1) * u[i]*u[j] / u2 + (i == j);
    }
  }
  for (int i=0; i<4; ++i) {
    xp[i] = 0.0;
    for (int j=0; j<4; ++j) {
      xp[i] += L[i][j] * x[j];
    }
  }
}
