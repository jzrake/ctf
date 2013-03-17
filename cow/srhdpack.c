
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cow.h"
#include "srhdpack.h"
#define MODULE "srhdpack"

static void boost(double u[4], double x[4], double xp[4]);
static double len3(double *x);
static double gamm(double *x);
static double dot3(double *u, double *v);
static double dot4(double *u, double *v);
static void project3(double *dx, double *u, double *ulat, double *ulon);

#define rho args[0][0]
#define pre args[0][1]
#define  vx args[0][2]
#define  vy args[0][3]
#define  vz args[0][4]

static void reduce_gamma_rho(double *result, double **args, int **strides,
			     void *udata)
{
  double v2 = vx*vx + vy*vy + vz*vz;
  *result = rho / sqrt(1.0 - v2);
}
static void reduce_proper_rho(double *result, double **args, int **strides,
			      void *udata)
{
  *result = rho;
}
static void reduce_gamma_beta(double *result, double **args, int **strides,
			       void *udata)
{
  double v2 = vx*vx + vy*vy + vz*vz;
  *result = sqrt(v2 / (1.0 - v2));
}

#undef rho
#undef pre
#undef vx
#undef vy
#undef vz



char *srhdpack_onepointpdfs(cow_dfield *prim,
			    char *which,
			    char *h5fname,
			    char *h5gname,
			    double *interval)
{
  int nbin = 4096;
  double mms[3]; // min, max, sum
  cow_transform op;

  if (strcmp(which, "proper-rho") == 0) {
    op = reduce_proper_rho;
  }
  else if (strcmp(which, "gamma-rho") == 0) {
    op = reduce_gamma_rho;
  }
  else if (strcmp(which, "gamma-beta") == 0) {
    op = reduce_gamma_beta;
  }
  else {
    return "'which' not in [proper-rho, gamma-rho, gamma-beta]";
  }

  cow_dfield_settransform(prim, op);

  if (interval) { // use manual interval if given
    mms[0] = interval[0];
    mms[1] = interval[1];
  }
  else { // otherwise set the interval automatically
    cow_dfield_reduce(prim, mms);
  }

  cow_histogram *h = cow_histogram_new();
  cow_histogram_setnbins(h, 0, nbin);
  cow_histogram_setlower(h, 0, mms[0]);
  cow_histogram_setupper(h, 0, mms[1]);
  cow_histogram_setspacing(h, COW_HIST_SPACING_LOG);
  cow_histogram_setbinmode(h, COW_HIST_BINMODE_COUNTS);
  cow_histogram_commit(h);
  cow_histogram_populate(h, prim, op);
  cow_histogram_seal(h);
  cow_histogram_setnickname(h, which);
  cow_histogram_dumphdf5(h, h5fname, h5gname);
  cow_histogram_del(h);

  return NULL;
}

void srhdpack_collectpairs(cow_dfield *vel,
			   srhdpack_samplemode *modes,
			   int num_modes,
			   int num_pairs,
			   int num_samps,
			   double *samploc,
			   double *outbufx,
			   double *outbufy)
/* -----------------------------------------------------------------------------
 * DESCRIPTION:
 *
 * This function helps compute pairwise structure functions by sampling the data
 * field `vel` at the locations in `samploc`, and processing them in as many
 * ways as there are `modes`. For each pair, the function loops over `modes`
 * processing it and appending the x (separation) and y (sample weight) result
 * to `outbufx` and `outbufy` respectively.
 *
 * ARGUMENTS:
 *
 * vel       ... I 3-member, 3 dimensional cow_dfield of 3-velocity
 * modes     ... I array of modes for processing pair samples
 * num_modes ... I size of `modes` array
 * num_pairs ... I number of pairs to select and process
 * num_samps ... I size of one-point pool of samples from which to select pairs
 * samploc   ... I buffer (num_samps x 3) of sample coordinates
 * outbufx   ... O buffer (num_pairs x num_modes) to hold separations 
 * outbufy   ... O buffer (num_pairs x num_modes) to hold sample weights
 * -----------------------------------------------------------------------------
 */
{
  int ntotmeas = 0;
  double *x = (double*) malloc(num_samps * 3 * sizeof(double));
  double *v;

  for (int n=0; n<num_samps; ++n) {
    x[3*n + 0] = samploc[3*n + 0];
    x[3*n + 1] = samploc[3*n + 1];
    x[3*n + 2] = samploc[3*n + 2];
  }

  cow_dfield_setsamplemode(vel, COW_SAMPLE_LINEAR);
  int err = cow_dfield_setsamplecoords(vel, x, num_samps, 3);
  if (err) {
    printf("[%s] Error! setsamplecoords returned %d\n", MODULE, err);
    return;
  }
  free(x);

  int nout1, nout2;
  cow_dfield_sampleexecute(vel);
  cow_dfield_getsamplecoords(vel, &x, &nout1, NULL);
  cow_dfield_getsampleresult(vel, &v, &nout2, NULL);

  double x1[3];
  double x2[3];

  for (int n=0; n<num_pairs; ++n) {
    int i1 = rand() % num_samps;
    int i2 = rand() % num_samps;
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
    double dumu[4] = { umu2[0] - umu1[0],
		       umu2[1] - umu1[1],
		       umu2[2] - umu1[2],
		       umu2[3] - umu1[3] };
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
    project3(dxpro, &dupro[1], &dulat[1], &dulon[1]);

    dvlat[0] = len3(&dvlat[1]);
    dvlon[0] = len3(&dvlon[1]);
    dulat[0] = gamm(&dulat[1]);
    dulon[0] = gamm(&dulon[1]);

    double xvalue;
    double yvalue;
    double gammarel;
    double betarel;
    double gammabetarel;

    for (int m=0; m<num_modes; ++m) {

      switch (modes[m].sepmode) {
      case SRHDPACK_SEPARATION_PROPER:
        xvalue = drpro;
        break;
      case SRHDPACK_SEPARATION_LAB:
        xvalue = drlab;
        break;
      default:
        printf("[%s] Error! invalid argument: sepmode %d\n", MODULE,
	       modes[m].sepmode);
        return;
      }

      switch (modes[m].projmode) {
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
        printf("[%s] Error! invalid argument: projmode %d\n", MODULE,
	       modes[m].projmode);
        return;
      }

      switch (modes[m].velmode) {
      case SRHDPACK_VELOCITY_GAMMA:
        yvalue = gammarel;
        break;
      case SRHDPACK_VELOCITY_BETA:
        yvalue = betarel;
        break;
      case SRHDPACK_VELOCITY_GAMMABETA:
        yvalue = gammabetarel;
        break;
      case SRHDPACK_VELOCITY_DUMUDXMU:
	yvalue = dot4(dumu, dxlab) / sqrt(dot4(dxlab, dxlab));
	break;
      case SRHDPACK_VELOCITY_DUMUDUMU:
	yvalue = dot4(dumu, dumu);
	break;
      default:
        printf("[%s] Error! invalid argument: velmode %d\n", MODULE,
	       modes[m].velmode);
        return;
      }

      outbufx[ntotmeas] = xvalue;
      outbufy[ntotmeas] = pow(fabs(yvalue), modes[m].exponent);
      ntotmeas += 1;
    }
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
double dot4(double *u, double *v)
{
  return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2] + u[3]*v[3];
}
double gamm(double *x)
{
  return 1.0 / sqrt(1.0 - dot3(x, x));
}
void project3(double *dx, double *u, double *ulat, double *ulon)
// -----------------------------------------------------------------------------
// Returns the projection of the vector u onto the vector dx in ulat, and its
// projection onto the plane normal to dx in ulon.
// -----------------------------------------------------------------------------
{
  double dxlen = len3(dx); if (dxlen<1e-12) dxlen = 1.0; // don't divide by zero
  double dxhat[3] = { dx[0] / dxlen, dx[1] / dxlen, dx[2] / dxlen };
  double ulonlen = dot3(dxhat, u);
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
