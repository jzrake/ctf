
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#define COW_PRIVATE_DEFS
#include "cow.h"

#define MODULE "sampling"
#define EPS 1e-12

static void _sample1(cow_dfield *f, double *x, double *P, int mode);
static void _sample2(cow_dfield *f, double *x, double *P, int mode);
static void _sample3(cow_dfield *f, double *x, double *P, int mode);
static void _loc(cow_dfield *f, double *Ri, int Nsamp, double *Ro, double *Po,
		 int mode);
static void _rem(cow_dfield *f, double *Ri, int Nsamp, double *Ro, double *Po,
                 int mode);

int cow_dfield_setsamplecoords(cow_dfield *f, double *x, int ns, int nd)
{
  double *gx0 = f->domain->glb_lower;
  double *gx1 = f->domain->glb_upper;
  if (nd != 3) {
    return COW_SAMPLE_ERROR_WRONGD; // nd must be 3
  }
  for (int m=0; m<ns; ++m) {
    double *r1 = &x[3*m];
    for (int d=0; d<f->domain->n_dims; ++d) {
      if (gx0[d] > r1[d] - EPS || r1[d] + EPS > gx1[d]) {
	printf("[%s] warning: sample is off-bounds (%f %f %f)\n",
	       MODULE, r1[0], r1[1], r1[2]);
	return COW_SAMPLE_ERROR_OUT; // sample out of bounds
      }
    }
  }
  int m = f->n_members;
  f->samplecoords = (double*) realloc(f->samplecoords, ns * 3 * sizeof(double));
  f->sampleresult = (double*) realloc(f->sampleresult, ns * m * sizeof(double));
  f->samplecoordslen = ns;
  memcpy(f->samplecoords, x, ns * 3 * sizeof(double));
  return 0;
}
void cow_dfield_getsamplecoords(cow_dfield *f, double **x, int *ns, int *nd)
{
  if (ns) *ns = f->samplecoordslen;
  if (nd) *nd = 3;
  if (x) *x = f->samplecoords;
}
void cow_dfield_getsampleresult(cow_dfield *f, double **P, int *ns, int *nd)
{
  if (ns) *ns = f->samplecoordslen;
  if (nd) *nd = f->n_members;
  if (P) *P = f->sampleresult;
}
void cow_dfield_setsamplemode(cow_dfield *f, int mode)
{
  f->samplemode = mode;
}
void cow_dfield_sampleexecute(cow_dfield *f)
// -----------------------------------------------------------------------------
// Samples `N` points on the global domain. `N` may vary between processes, or
// be equal to zero. Points are output permuted relative to the input
// coordinates, so take care to use `xout` as the coordinates associated with
// samples in `P`.
//
// f:    IN   data field instance to sample
// x:    IN   list of input coordinates at which to sample f's data (N x 3)
// N:    IN   number of points to sample
// xout: OUT  locations of returned samples, permutation of x (N x 3)
// P:    OUT  list of filled samples (N x Q) where Q = f->n_members
// -----------------------------------------------------------------------------
{
  if (cow_domain_getguard(f->domain) < 1 &&
      f->samplemode == COW_SAMPLE_LINEAR) {
    printf("[%s] warning: sample mode COW_SAMPLE_LINEAR requires at least one"
	   " guard zone, sampleexecute function will have no effect\n", MODULE);
    return;
  }
  double *xout = (double*) malloc(f->samplecoordslen * 3 * sizeof(double));
  double *xin = f->samplecoords;
  int N = f->samplecoordslen;
  double *P = f->sampleresult;
  int mode = f->samplemode;
  if (cow_mpirunning()) {
    _rem(f, xin, N, xout, P, mode);
  }
  else {
    _loc(f, xin, N, xout, P, mode);
  }
  memcpy(f->samplecoords, xout, f->samplecoordslen * 3 * sizeof(double));
  free(xout);
}

void cow_dfield_sampleglobalind(cow_dfield *f, int i, int j, int k, double **x,
				int *n0)
{
  int ng = f->domain->n_ghst;
  double xin[3];
  i -= cow_domain_getglobalstartindex(f->domain, 0);
  j -= cow_domain_getglobalstartindex(f->domain, 1);
  k -= cow_domain_getglobalstartindex(f->domain, 2);
  xin[0] = cow_domain_positionatindex(f->domain, 0, i + ng);
  xin[1] = cow_domain_positionatindex(f->domain, 1, j + ng);
  xin[2] = cow_domain_positionatindex(f->domain, 2, k + ng);
  cow_dfield_setsamplecoords(f, xin, 1, 3);
  cow_dfield_setsamplemode(f, COW_SAMPLE_NEAREST);
  cow_dfield_sampleexecute(f);
  cow_dfield_getsampleresult(f, x, NULL, n0);
}

/* -----------------------------------------------------------------------------
 *
 *
 *               Description of method used by COW_SAMPLE_LINEAR
 *
 *
 *     @@@@@@@@@@@@@@@@@@@@@@@         Sampling use a centered gradient, so
 *     @ x := zone center    @         points sampled exactly at zone centers
 *     @ | := zone interface @         will in general not give the exact value
 *     @ o := target point   @         associated with that zone.
 *     @@@@@@@@@@@@@@@@@@@@@@@
 *
 *     { <---------  2*dx[0]  ----------> }
 *     ====================================
 *              |               |
 *     x        |       x    o  |         x
 *     x        |       x    o  |         x
 *              |               |
 *     ====================================
 *     ^                ^    ^           ^
 *     x[i-1]           i    x[0]        x[xi+1]
 *     P0                                P1
 *
 * -----------------------------------------------------------------------------
 */
void _sample1(cow_dfield *f, double *x, double *P, int mode)
{
#define M(i) ((i)*s[0])
  cow_domain *d = f->domain;
  int *s = f->stride;
  int i = cow_domain_indexatposition(d, 0, x[0]);
  double *A = (double*) f->data;
  if (mode == COW_SAMPLE_NEAREST) {
    memcpy(P, A + M(i), f->n_members * sizeof(double));
  }
  else if (mode == COW_SAMPLE_LINEAR) {
    double x0 = cow_domain_positionatindex(d, 0, i-1);
    double *P0 = &A[M(i-1)];
    double *P1 = &A[M(i+1)];
    double delx[1] = { 0.5 * (x[0] - x0) / d->dx[0] };
    for (int q=0; q<f->n_members; ++q) {
      double b1 = P0[q];
      double b2 = P1[q] - P0[q];
      P[q] = b1 + b2*delx[0];
    }
  }
#undef M
}
void _sample2(cow_dfield *f, double *x, double *P, int mode)
{
#define M(i,j) ((i)*s[0] + (j)*s[1])
  cow_domain *d = f->domain;
  int *s = f->stride;
  int i = cow_domain_indexatposition(d, 0, x[0]);
  int j = cow_domain_indexatposition(d, 1, x[1]);
  double *A = (double*) f->data;
  if (mode == COW_SAMPLE_NEAREST) {
    memcpy(P, A + M(i,j), f->n_members * sizeof(double));
  }
  else if (mode == COW_SAMPLE_LINEAR) {
    double x0 = cow_domain_positionatindex(d, 0, i-1);
    double y0 = cow_domain_positionatindex(d, 1, j-1);
    double *P00 = &A[M(i-1,j-1)];
    double *P01 = &A[M(i-1,j+1)];
    double *P10 = &A[M(i+1,j-1)];
    double *P11 = &A[M(i+1,j+1)];
    double delx[2] = {
      0.5 * (x[0] - x0) / d->dx[0],
      0.5 * (x[1] - y0) / d->dx[1] };
    for (int q=0; q<f->n_members; ++q) {
      double b1 = P00[q];
      double b2 = P10[q] - P00[q];
      double b3 = P01[q] - P00[q];
      double b4 = P00[q] - P10[q] - P01[q] + P11[q];
      P[q] = b1 + b2*delx[0] + b3*delx[1] + b4*delx[0]*delx[1];
    }
  }
#undef M
}
void _sample3(cow_dfield *f, double *x, double *P, int mode)
{
#define M(i,j,k) ((i)*s[0] + (j)*s[1] + (k)*s[2])
  cow_domain *d = f->domain;
  int *s = f->stride;
  int i = cow_domain_indexatposition(d, 0, x[0]);
  int j = cow_domain_indexatposition(d, 1, x[1]);
  int k = cow_domain_indexatposition(d, 2, x[2]);
  double *A = (double*) f->data;
  if (mode == COW_SAMPLE_NEAREST) {
    memcpy(P, A + M(i,j,k), f->n_members * sizeof(double));
  }
  else if (mode == COW_SAMPLE_LINEAR) {
    double x0 = cow_domain_positionatindex(d, 0, i-1);
    double y0 = cow_domain_positionatindex(d, 1, j-1);
    double z0 = cow_domain_positionatindex(d, 2, k-1);
    double *P000 = &A[M(i-1,j-1,k-1)];
    double *P001 = &A[M(i-1,j-1,k+1)];
    double *P010 = &A[M(i-1,j+1,k-1)];
    double *P011 = &A[M(i-1,j+1,k+1)];
    double *P100 = &A[M(i+1,j-1,k-1)];
    double *P101 = &A[M(i+1,j-1,k+1)];
    double *P110 = &A[M(i+1,j+1,k-1)];
    double *P111 = &A[M(i+1,j+1,k+1)];
    double delx[3] = {
      0.5 * (x[0] - x0) / d->dx[0],
      0.5 * (x[1] - y0) / d->dx[1],
      0.5 * (x[2] - z0) / d->dx[2] };
    // -------------------------------------------------------------------------
    // See http://en.wikipedia.org/wiki/Trilinear_interpolation
    // -------------------------------------------------------------------------
    for (int q=0; q<f->n_members; ++q) {
      double i1 = P000[q] * (1.0 - delx[2]) + P001[q] * delx[2];
      double i2 = P010[q] * (1.0 - delx[2]) + P011[q] * delx[2];
      double j1 = P100[q] * (1.0 - delx[2]) + P101[q] * delx[2];
      double j2 = P110[q] * (1.0 - delx[2]) + P111[q] * delx[2];
      double w1 = i1 * (1.0 - delx[1]) + i2 * delx[1];
      double w2 = j1 * (1.0 - delx[1]) + j2 * delx[1];
      P[q] = w1 * (1.0 - delx[0]) + w2 * delx[0];
    }
  }
#undef M
}

void _loc(cow_dfield *f, double *Ri, int Nsamp, double *Ro, double *Po,
          int mode)
{
  int Q = f->n_members;
  memcpy(Ro, Ri, Nsamp * 3 * sizeof(double));
  for (int n=0; n<Nsamp; ++n) {
    switch (f->domain->n_dims) {
    case 1: _sample1(f, &Ri[3*n], &Po[Q*n], mode); break;
    case 2: _sample2(f, &Ri[3*n], &Po[Q*n], mode); break;
    case 3: _sample3(f, &Ri[3*n], &Po[Q*n], mode); break;
    default: break;
    }
  }
}

void _rem(cow_dfield *f, double *Ri, int Nsamp, double *Ro, double *Po,
          int mode)
{
#if (COW_MPI)
  int Q = f->n_members;
  int rank = f->domain->cart_rank;
  int size = f->domain->cart_size;
  int Nd = f->domain->n_dims;
  double **remote_r1 = (double**) malloc(size * sizeof(double*));
  double **remote_P1 = (double**) malloc(size * sizeof(double*));
  int *remote_r1_size = (int*) malloc(size * sizeof(int));
  int *remote_P1_size = (int*) malloc(size * sizeof(int));
  for (int n=0; n<size; ++n) {
    remote_r1[n] = NULL;
    remote_P1[n] = NULL;
    remote_r1_size[n] = 0;
    remote_P1_size[n] = 0;
  }
  for (int m=0; m<Nsamp; ++m) {
    double *r1 = &Ri[3*m];
    int remote = cow_domain_subgridatposition(f->domain, r1[0], r1[1], r1[2]);
    double **remR = &remote_r1[remote];
    double **remP = &remote_P1[remote];
    for (int d=0; d<3; ++d) {
      int N = (remote_r1_size[remote] += 1);
      remR[0] = (double*) realloc(remR[0], N * sizeof(double));
      remR[0][N - 1] = r1[d];
    }
    for (int q=0; q<Q; ++q) {
      int N = (remote_P1_size[remote] += 1);
      remP[0] = (double*) realloc(remP[0], N * sizeof(double));
      remP[0][N - 1] = 0.0;
    }
  }
  // ---------------------------------------------------------------------------
  // We will be sampling pairs between ourselves and the remote process, rank +
  // dn. That process is referred to as 'lawyer' because they will work for us,
  // obtaining the records we have determined live on their domain. Similarly,
  // we will be the lawyer for process rank-dn, so that process is called
  // 'client'.
  // ---------------------------------------------------------------------------
  MPI_Status status;
  MPI_Comm comm = f->domain->mpi_cart;
  int queries_satisfied = 0;
  for (int dn=0; dn<size; ++dn) {
    // -------------------------------------------------------------------------
    // This loop contains three pairs of matching Send/Recv's. For the first
    // two, the send is place to the lawyer process, and the receive comes from
    // the client. We call our lawyer and let him know to expect a message of
    // length 'numlawyer' double[3]'s. Those are the positions of the remote
    // points we have chosen, and then determined live on his domain. At the
    // same time, we receive a message from our client, which contains the
    // length 'numclient' of double[3]'s he will be asking us to fetch. The
    // next pair of Send/Recv's is the list of coordinates themselves. The last
    // one is the list of primitive quantities at those locations.
    // -------------------------------------------------------------------------
    int lawyer = (rank + size + dn) % size;
    int client = (rank + size - dn) % size;
    int numlawyer = remote_r1_size[lawyer] / 3;
    int numclient;
    MPI_Sendrecv(&numlawyer, 1, MPI_INT, lawyer, 123,
                 &numclient, 1, MPI_INT, client, 123, comm, &status);
    double *you_get_for_me_r = remote_r1[lawyer];
    double *you_get_for_me_P = remote_P1[lawyer];
    double *I_find_for_you_r = (double*) malloc(numclient * 3 * sizeof(double));
    double *I_find_for_you_P = (double*) malloc(numclient * Q * sizeof(double));
    MPI_Sendrecv(you_get_for_me_r, numlawyer * 3, MPI_DOUBLE, lawyer, 123,
                 I_find_for_you_r, numclient * 3, MPI_DOUBLE, client, 123,
                 comm, &status);
    for (int s=0; s<numclient; ++s) {
      double *r_query = &I_find_for_you_r[3*s];
      double *Panswer = (double*) malloc(Q * sizeof(double));
      switch (Nd) {
      case 1: _sample1(f, r_query, Panswer, mode);
      case 2: _sample2(f, r_query, Panswer, mode);
      case 3: _sample3(f, r_query, Panswer, mode);
      }
      memcpy(&I_find_for_you_P[s*Q], &Panswer[0], Q*sizeof(double));
      free(Panswer);
    }
    MPI_Sendrecv(I_find_for_you_P, numclient*Q, MPI_DOUBLE, client, 123,
                 you_get_for_me_P, numlawyer*Q, MPI_DOUBLE, lawyer, 123,
                 comm, &status);
    memcpy(&Ro[queries_satisfied * 3], remote_r1[lawyer],
           remote_r1_size[lawyer] * sizeof(double));
    memcpy(&Po[queries_satisfied * Q], remote_P1[lawyer],
           remote_P1_size[lawyer] * sizeof(double));
    queries_satisfied += numlawyer;
    free(I_find_for_you_r);
    free(I_find_for_you_P);
  }
  for (int n=0; n<size; ++n) {
    free(remote_r1[n]);
    free(remote_P1[n]);
  }
  free(remote_r1);
  free(remote_P1);
  free(remote_r1_size);
  free(remote_P1_size);
#endif
}
