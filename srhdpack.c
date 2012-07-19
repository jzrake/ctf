
#include <stdio.h>
#include <math.h>
#include "cow.h"
#define MODULE "srhdpack"

static void _dobatch(cow_dfield *vel, cow_histogram *hist, int N, char mode);
static void _boost(double u[4], double x[4], double xp[4]);

void srhdpack_relativelorentzpairs(cow_dfield *vel,
				   cow_histogram *histpro,
				   cow_histogram *histlab,
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
//  cow_histogram_setnickname(hist, "mypspec"); // optional
//  cow_histogram_setlower(histpro, 0, 0.1); // good values for most cases
//  cow_histogram_setupper(histpro, 0, 1.5);
//
// The histogram `histpro` bins the pairs in terms of their proper space-like
// separation, while `histlab` uses their lab-frame separations.
//
// NOTE: this function generates its own sample points by calling rand(). If
// running in parallel, the user is responsible for making sure to seed the
// ranks differently, unless `seed` is set to 1, in which case srand will be
// called with the domain's cartesian communicator rank.
// -----------------------------------------------------------------------------
{
  if (seed) {
    srand(cow_domain_getcartrank(cow_dfield_getdomain(vel)));
  }

  cow_histogram_setbinmode(histpro, COW_HIST_BINMODE_AVERAGE);
  cow_histogram_setbinmode(histlab, COW_HIST_BINMODE_AVERAGE);
  cow_histogram_commit(histpro);
  cow_histogram_commit(histlab);

  for (int n=0; n<nbatch; ++n) {
    _dobatch(vel, histlab, nperbatch, 'l');
    _dobatch(vel, histpro, nperbatch, 'p');
    printf("[%s] running batch %d/%d\n", MODULE, n, nbatch);
  }
  cow_histogram_seal(histpro);
  cow_histogram_seal(histlab);
}


void _boost(double u[4], double x[4], double xp[4])
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


void _dobatch(cow_dfield *vel, cow_histogram *hist, int N, char mode)
{
  double *x = (double*) malloc(N * 3 * sizeof(double));
  double *v;

  for (int n=0; n<N; ++n) {
    x[3*n + 0] = (double) rand() / RAND_MAX;
    x[3*n + 1] = (double) rand() / RAND_MAX;
    x[3*n + 2] = (double) rand() / RAND_MAX;
  }

  cow_dfield_setsamplemode(vel, COW_SAMPLE_LINEAR);
  int err = cow_dfield_setsamplecoords(vel, x, N, 3);
  if (err) {
    printf("error on setsamplecoords: %d\n", err);
  }
  free(x);

  int nout1, nout2;
  cow_dfield_sampleexecute(vel);
  cow_dfield_getsamplecoords(vel, &x, &nout1, NULL);
  cow_dfield_getsampleresult(vel, &v, &nout2, NULL);

  for (int n=0; n<nout1/2; ++n) {
    double *x1 = &x[6*n + 0];
    double *x2 = &x[6*n + 3];
    double *v1 = &v[6*n + 0];
    double *v2 = &v[6*n + 3];
    double g1 = 1.0 / sqrt(1.0 - (v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]));
    double g2 = 1.0 / sqrt(1.0 - (v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]));
    double umu1[4] = { g1, g1*v1[0], g1*v1[1], g1*v1[2] };
    double umu2[4] = { g2, g2*v2[0], g2*v2[1], g2*v2[2] };
    double xmu1[4] = { 0.0, x1[0], x1[1], x1[2] };
    double xmu2[4] = { 0.0, x2[0], x2[1], x2[2] };
    double xmu1p[4];
    double xmu2p[4];
    double xrel[4];
    double urel[4];
    _boost(umu1, xmu1, xmu1p);
    _boost(umu1, xmu2, xmu2p);
    _boost(umu1, umu2, urel);
    xrel[0] = xmu2p[0] - xmu1p[0];
    xrel[1] = xmu2p[1] - xmu1p[1];
    xrel[2] = xmu2p[2] - xmu1p[2];
    xrel[3] = xmu2p[3] - xmu1p[3];
    double dx[3] = { x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2] };
    double drlab = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
    double drprop = sqrt(xrel[1]*xrel[1] + xrel[2]*xrel[2] + xrel[3]*xrel[3]);
    double gammarel = urel[0];

    if (mode == 'p') { // proper
      cow_histogram_addsample1(hist, drprop, gammarel);
    }
    else if (mode == 'l') {
      cow_histogram_addsample1(hist, drlab, gammarel);
    }
  }
}
