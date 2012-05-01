
// -----------------------------------------------------------------------------
//
// TODO:
//
// Verify Cees and Dees, Betas had a typo in Shasha's paper
//
// -----------------------------------------------------------------------------


#include "plm-split.hpp"

#define MAXNQ 8 // Used for static array initialization

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

static const double DeesC2L[3] = {   1./ 16.,   5./  8.,   5./ 16. };
static const double DeesC2R[3] = {   5./ 16.,   5./  8.,   1./ 16. };
static const double DeesA2C[3] = {  -9./ 80.,  49./ 40.,  -9./ 80. };
static const double DeesC2A[3] = { -17./240., 137./120., -17./240. };


static inline double weno5(const double *v, const double c[3][3], const double d[3])
{
  const double eps = 1e-5; // 1e-6 -> 1e-4
  const int k = 3;
  double w[3], vs[3];

  const double B[3] = {
    (13.0/12.0)*pow(  v[ 0] - 2*v[ 1] +   v[ 2], 2.0) +
    ( 1.0/ 4.0)*pow(3*v[ 0] - 4*v[ 1] +   v[ 2], 2.0),

    (13.0/12.0)*pow(  v[-1] - 2*v[ 0] +   v[ 1], 2.0) +
    ( 1.0/ 4.0)*pow(  v[-1] - 0*v[ 0] -   v[ 1], 2.0),

    (13.0/12.0)*pow(  v[-2] - 2*v[-1] +   v[ 0], 2.0) +
    ( 1.0/ 4.0)*pow(  v[-2] - 4*v[-1] + 3*v[ 0], 2.0)
  };

  double wtot = 0.0;
  for (int r=0; r<k; ++r) {
    vs[r] = 0.0;
    for (int j=0; j<k; ++j) {
      vs[r] += c[r][j] * v[j-r];
    }
    w[r] = d[r] / pow(eps + B[r], 2);
    wtot += w[r];
  }

  double viph = 0.0;
  for (int r=0; r<k; ++r) {
    viph += (w[r] / wtot) * vs[r];
  }
  return viph;
}


std::valarray<double> Deriv::dUdt(const std::valarray<double> &Uin)
{
  this->prepare_integration();

  std::valarray<double> U = Uin;
  std::valarray<double> L(U.size());
  std::valarray<double> &P = Mara->PrimitiveArray;

  ConsToPrim(U, P);
  DriveSweeps(P, L);

  return L;
}
void Deriv::DriveSweeps(const std::valarray<double> &P,
                        std::valarray<double> &L)
{
  switch (ND) {
  case 1: drive_sweeps_1d(&P[0], &L[0]); break;
  case 2: drive_sweeps_2d(&P[0], &L[0]); break;
  case 3: drive_sweeps_3d(&P[0], &L[0]); break;
  }
}
void Deriv::intercell_flux_sweep(const double *P, double *F, int dim)
{
  int i,S=stride[dim];
  double Pl[MAXNQ], Pr[MAXNQ];

  for (i=2*S; i<stride[0]-3*S; i+=NQ) {

    for (int q=0; q<NQ; ++q) {
      const int m = i + q;
      double vm[5] = { P[m-2*S], P[m-S], P[m+0], P[m+1*S], P[m+2*S] };
      double vp[5] = { P[m-1*S], P[m-0], P[m+S], P[m+2*S], P[m+3*S] };

      //      Pl[q] = P[m];
      //      Pr[q] = P[m+S];

      Pl[q] = weno5(vm+2, CeesC2R, DeesC2R);
      Pr[q] = weno5(vp+2, CeesC2L, DeesC2L);
    }

    Mara->riemann->IntercellFlux(Pl, Pr, 0, &F[i], 0.0, dim);
    // remember to catch errors here
  }
}

void Deriv::drive_sweeps_1d(const double *P, double *L)
{
  double *F = (double*) malloc(stride[0]*sizeof(double));

  int i,sx=stride[1];
  intercell_flux_sweep(P,F,1);

  for (i=sx; i<stride[0]; ++i) {
    L[i] = -(F[i]-F[i-sx])/dx;
  }
  free(F);
}
void Deriv::drive_sweeps_2d(const double *P, double *L)
{
  double *F = (double*) malloc(stride[0]*sizeof(double));
  double *G = (double*) malloc(stride[0]*sizeof(double));

  int i,sx=stride[1],sy=stride[2];
  intercell_flux_sweep(P,F,1);
  intercell_flux_sweep(P,G,2);

  Mara->fluid->ConstrainedTransport2d(F,G,stride);

  for (i=sx; i<stride[0]; ++i) {
    L[i] = -(F[i]-F[i-sx])/dx - (G[i]-G[i-sy])/dy;
  }
  free(F); free(G);
}
void Deriv::drive_sweeps_3d(const double *P, double *L)
{
  double *F = (double*) malloc(stride[0]*sizeof(double));
  double *G = (double*) malloc(stride[0]*sizeof(double));
  double *H = (double*) malloc(stride[0]*sizeof(double));

  int i,sx=stride[1],sy=stride[2],sz=stride[3];
  intercell_flux_sweep(P,F,1);
  intercell_flux_sweep(P,G,2);
  intercell_flux_sweep(P,H,3);

  Mara->fluid->ConstrainedTransport3d(F,G,H,stride);

  for (i=sx; i<stride[0]; ++i) {
    L[i] = -(F[i]-F[i-sx])/dx - (G[i]-G[i-sy])/dy - (H[i]-H[i-sz])/dz;
  }
  free(F); free(G); free(H);
}
