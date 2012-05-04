
// -----------------------------------------------------------------------------
//
// TODO:
//
// Verify Cees and Dees, Betas had a typo in Shasha's paper
//
// -----------------------------------------------------------------------------


#include "plm-split.hpp"

#define MAXNQ 8 // Used for static array initialization
#define SQU(x) ((x)*(x))


typedef Weno5RiemannMethodOfLinesSplit Deriv;

static const double CeesA2C[3][3] = { {23./24.,  1./12.,  -1./24.},
                                      {-1./24., 13./12.,  -1./24.},
                                      {-1./24.,  1./12.,  23./24.} };
static const double CeesC2A[3][3] = { {25./24., -1./12.,   1./24.},
                                      { 1./24., 11./12.,   1./24.},
                                      { 1./24., -1./12.,  25./24.} };
static const double CeesC2L[3][3] = { {15./8., -5./4.,  3./8.},
                                      { 3./8.,  3./4., -1./8.},
                                      {-1./8.,  3./4.,  3./8.} };
static const double CeesC2R[3][3] = { { 3./8., 3./4.,  -1./8.},
                                      {-1./8., 3./4.,   3./8.},
                                      { 3./8.,-5./4.,  15./8.} };

static const double CeesC2R_A[3][3] = { { 1./3.,  5./6., -1./6. },
					{-1./6.,  5./6.,  1./3. },
					{ 1./3., -7./6., 11./6. } };
static const double CeesC2L_A[3][3] = { {11./6., -7./6.,  1./3. },
					{ 1./3.,  5./6., -1./6. },
					{-1./6.,  5./6.,  1./3. } };

static const double DeesC2L[3] = {   1./ 16.,   5./  8.,   5./ 16. };
static const double DeesC2R[3] = {   5./ 16.,   5./  8.,   1./ 16. };
static const double DeesA2C[3] = {  -9./ 80.,  49./ 40.,  -9./ 80. };
static const double DeesC2A[3] = { -17./240., 137./120., -17./240. };

static const double DeesC2L_A[3] = { 0.1, 0.6, 3.0 };
static const double DeesC2R_A[3] = { 0.3, 0.6, 1.0 };


static inline double weno5(const double *v, const double c[3][3], const double d[3])
{
  const double eps = 1e-16; // 1e-6 -> 1e-4

  const double B[3] = {
    (13.0/12.0)*SQU(  v[ 0] - 2*v[ 1] +   v[ 2]) +
    ( 1.0/ 4.0)*SQU(3*v[ 0] - 4*v[ 1] +   v[ 2]),

    (13.0/12.0)*SQU(  v[-1] - 2*v[ 0] +   v[ 1]) +
    ( 1.0/ 4.0)*SQU(  v[-1] - 0*v[ 0] -   v[ 1]),

    (13.0/12.0)*SQU(  v[-2] - 2*v[-1] +   v[ 0]) +
    ( 1.0/ 4.0)*SQU(  v[-2] - 4*v[-1] + 3*v[ 0])
  };

  const double vs[3] = {
    c[0][0]*v[ 0] + c[0][1]*v[ 1] + c[0][2]*v[2],
    c[1][0]*v[-1] + c[1][1]*v[ 0] + c[1][2]*v[1],
    c[2][0]*v[-2] + c[2][1]*v[-1] + c[2][2]*v[0],
  };
  const double w[3] = {
    d[0] / SQU(eps + B[0]),
    d[1] / SQU(eps + B[1]),
    d[2] / SQU(eps + B[2])
  };

  const double wtot = w[0] + w[1] + w[2];
  return (w[0]*vs[0] + w[1]*vs[1] + w[2]*vs[2])/wtot;
}

std::valarray<double> Deriv::dUdt(const std::valarray<double> &Uin)
{
  this->prepare_integration();

  std::valarray<double> U = Uin;
  std::valarray<double> L(U.size());
  std::valarray<double> &P = Mara->PrimitiveArray;

  ConsToPrim(U, P);
  DriveSweeps(U, L);

  return L;
}
void Deriv::DriveSweeps(const std::valarray<double> &U,
                        std::valarray<double> &L)
{
  switch (ND) {
  case 1: drive_sweeps_1d(&U[0], &L[0]); break;
  case 2: drive_sweeps_2d(&U[0], &L[0]); break;
  case 3: drive_sweeps_3d(&U[0], &L[0]); break;
  }
}
#include <iostream>
void Deriv::intercell_flux_sweep(const double *U_avg, double *F, int dim)
// -----------------------------------------------------------------------------
// U_avg is the volume-averaged conserved array. It has to be WENO-interpolated
// to the point-centered value, once for each direction. ConsToPrim is then run
// on the point-centered value of U.
// -----------------------------------------------------------------------------
{
  int i,S=stride[dim];
  double Pl[MAXNQ], Pr[MAXNQ];

  std::valarray<double> G(stride[0]); // is the point-valued HLL flux at i+1/2
  std::valarray<double> U(stride[0]); // point-valued conserved state for this axis
  std::valarray<double> u(stride[0]);
  std::valarray<double> P = Mara->PrimitiveArray;


  for (i=2*S; i<stride[0]-2*S; i+=NQ) {
    for (int q=0; q<NQ; ++q) {
      const int m = i + q;
      const double *U0 = &U_avg[m];
      double v[5] = { U0[-2*S], U0[-S], U0[0], U0[S], U0[2*S] };
      U[m] = U_avg[m];
      //      U[m] = weno5(v+2, CeesA2C, DeesA2C);
    }

    int error = Mara->fluid->ConsToPrim(&U[i], &P[i]);
    if (error) {
      printf("got c2p error on zone %d!\n", i/NQ);
      std::cout << Mara->fluid->PrintCons(&U_avg[i]) << std::endl;
      std::cout << Mara->fluid->PrintCons(&U[i]) << std::endl;
      std::cout << Mara->fluid->PrintPrim(&P[i]) << std::endl;
      exit(2);
    }
  }
  for (i=4*S; i<stride[0]-5*S; i+=NQ) {
    for (int q=0; q<NQ; ++q) {
      const int m = i + q;
      double vm[5] = { P[m-2*S], P[m-S], P[m+0], P[m+1*S], P[m+2*S] };
      double vp[5] = { P[m-1*S], P[m-0], P[m+S], P[m+2*S], P[m+3*S] };

      Pl[q] = weno5(vm+2, CeesC2R, DeesC2R);
      Pr[q] = weno5(vp+2, CeesC2L, DeesC2L);
    }

    int error = Mara->riemann->IntercellFlux(Pl, Pr, 0, &G[i], 0.0, dim);
    if (error) printf("got riemann solver error!\n");
  }

  for (i=6*S; i<stride[0]-7*S; i+=NQ) {
    for (int q=0; q<NQ; ++q) {
      const int m = i + q;
      double v[5] = { G[m-2*S], G[m-S], G[m+0], G[m+1*S], G[m+2*S] };
      F[m] = G[m];
      //      F[m] = weno5(v+2, CeesC2A, DeesC2A);
    }
  }
}

void Deriv::drive_sweeps_1d(const double *U, double *L)
{
  double *F = (double*) malloc(stride[0]*sizeof(double));

  int i,sx=stride[1];
  intercell_flux_sweep(U,F,1);

  for (i=sx; i<stride[0]; ++i) {
    L[i] = -(F[i]-F[i-sx])/dx;
  }
  free(F);
}
void Deriv::drive_sweeps_2d(const double *U, double *L)
{
  double *F = (double*) malloc(stride[0]*sizeof(double));
  double *G = (double*) malloc(stride[0]*sizeof(double));

  int i,sx=stride[1],sy=stride[2];
  intercell_flux_sweep(U,F,1);
  intercell_flux_sweep(U,G,2);

  Mara->fluid->ConstrainedTransport2d(F,G,stride);

  for (i=sx; i<stride[0]; ++i) {
    L[i] = -(F[i]-F[i-sx])/dx - (G[i]-G[i-sy])/dy;
  }
  free(F); free(G);
}
void Deriv::drive_sweeps_3d(const double *U, double *L)
{
  double *F = (double*) malloc(stride[0]*sizeof(double));
  double *G = (double*) malloc(stride[0]*sizeof(double));
  double *H = (double*) malloc(stride[0]*sizeof(double));

  int i,sx=stride[1],sy=stride[2],sz=stride[3];
  intercell_flux_sweep(U,F,1);
  intercell_flux_sweep(U,G,2);
  intercell_flux_sweep(U,H,3);

  Mara->fluid->ConstrainedTransport3d(F,G,H,stride);

  for (i=sx; i<stride[0]; ++i) {
    L[i] = -(F[i]-F[i-sx])/dx - (G[i]-G[i-sy])/dy - (H[i]-H[i-sz])/dz;
  }
  free(F); free(G); free(H);
}
