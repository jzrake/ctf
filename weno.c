

#include <math.h>
#include "weno.h"

#define SQU(x) ((x)*(x))
#define SMOOTHNESS_ZHA 1


static const double CeesA2C_FV[3][3] = { {23./24.,  1./12.,  -1./24.},
                                         {-1./24., 13./12.,  -1./24.},
                                         {-1./24.,  1./12.,  23./24.} };
static const double CeesC2A_FV[3][3] = { {25./24., -1./12.,   1./24.},
                                         { 1./24., 11./12.,   1./24.},
                                         { 1./24., -1./12.,  25./24.} };
static const double CeesC2L_FV[3][3] = { {15./8., -5./4.,  3./8.},
                                         { 3./8.,  3./4., -1./8.},
                                         {-1./8.,  3./4.,  3./8.} };
static const double CeesC2R_FV[3][3] = { { 3./8., 3./4.,  -1./8.},
                                         {-1./8., 3./4.,   3./8.},
                                         { 3./8.,-5./4.,  15./8.} };
static const double DeesC2L_FV[3] = {   1./ 16.,   5./  8.,   5./ 16. };
static const double DeesC2R_FV[3] = {   5./ 16.,   5./  8.,   1./ 16. };
static const double DeesA2C_FV[3] = {  -9./ 80.,  49./ 40.,  -9./ 80. };
static const double DeesC2A_FV[3] = { -17./240., 137./120., -17./240. };



static const double CeesC2R_FD[3][3] = { { 1./3.,  5./6., -1./6. },
                                         {-1./6.,  5./6.,  1./3. },
                                         { 1./3., -7./6., 11./6. } };
static const double CeesC2L_FD[3][3] = { {11./6., -7./6.,  1./3. },
                                         { 1./3.,  5./6., -1./6. },
                                         {-1./6.,  5./6.,  1./3. } };
static const double DeesC2L_FD[3] = { 0.1, 0.6, 0.3 };
static const double DeesC2R_FD[3] = { 0.3, 0.6, 0.1 };
static double __weno5(const double *v, const double c[3][3], const double d[3]);

double weno5(const double *v, enum WenoOperation type)
{
  switch (type) {
  case WENO5_FD_C2L: return __weno5(v, CeesC2L_FD, DeesC2L_FD);
  case WENO5_FD_C2R: return __weno5(v, CeesC2R_FD, DeesC2R_FD);
  case WENO5_FV_C2L: return __weno5(v, CeesC2L_FV, DeesC2L_FV);
  case WENO5_FV_C2R: return __weno5(v, CeesC2R_FV, DeesC2R_FV);
  case WENO5_FV_A2C: return __weno5(v, CeesA2C_FV, DeesA2C_FV);
  case WENO5_FV_C2A: return __weno5(v, CeesC2A_FV, DeesC2A_FV);
  default: return 0.0;
  }
}


static inline double min3(const double *x)
{
  double x01 = x[0] < x[1] ? x[0] : x[1];
  return x01 < x[2] ? x01 : x[2];
}
static inline double max3(const double *x)
{
  double x01 = x[0] > x[1] ? x[0] : x[1];
  return x01 > x[2] ? x01 : x[2];
}
double __weno5(const double *v, const double c[3][3], const double d[3])
{
  const double eps = 1e-20;
  const double eps_prime = 1e-10;

  double B[3] = {
    (13.0/12.0)*SQU(  v[ 0] - 2*v[ 1] +   v[ 2]) +
    ( 1.0/ 4.0)*SQU(3*v[ 0] - 4*v[ 1] +   v[ 2]),

    (13.0/12.0)*SQU(  v[-1] - 2*v[ 0] +   v[ 1]) +
    ( 1.0/ 4.0)*SQU(  v[-1] - 0*v[ 0] -   v[ 1]),

    (13.0/12.0)*SQU(  v[-2] - 2*v[-1] +   v[ 0]) +
    ( 1.0/ 4.0)*SQU(  v[-2] - 4*v[-1] + 3*v[ 0])
  };

  if (SMOOTHNESS_ZHA) {
    const double A = 10.0;
    const double tau5 = fabs(B[0] - B[2]);
    B[0] = (B[0] + eps) / (B[0] + tau5 + eps);
    B[1] = (B[1] + eps) / (B[1] + tau5 + eps);
    B[2] = (B[2] + eps) / (B[2] + tau5 + eps);

    const double R0 = min3(B) / (max3(B) + eps_prime);
    B[0] = R0*A*min3(B) + B[0];
    B[1] = R0*A*min3(B) + B[1];
    B[2] = R0*A*min3(B) + B[2];
  }

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

