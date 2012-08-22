

#include <math.h>
#include "weno.h"

#define SQU(x) ((x)*(x))
#define SGN(x) (((x)>0)-((x)<0))


static double plm_theta = 2.0;    // [1 -> 2 (most aggressive)]
static double shenzha10_A = 50.0; // [0 -> ~100 (most aggressive)]
static enum SmoothnessIndicator IS_mode = OriginalJiangShu96;

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
static double __plm(const double *v, double sgn);

void reconstruct_set_plm_theta(double theta)
{
  plm_theta = theta;
}
void reconstruct_set_smoothness_indicator(enum SmoothnessIndicator IS)
{
  IS_mode = IS;
}
void reconstruct_set_shenzha10_A(double A)
{
  shenzha10_A = A; 
}
double reconstruct(const double *v, enum ReconstructOperation type)
{
  switch (type) {
  case WENO5_FD_C2L: return __weno5(v, CeesC2L_FD, DeesC2L_FD);
  case WENO5_FD_C2R: return __weno5(v, CeesC2R_FD, DeesC2R_FD);
  case WENO5_FV_C2L: return __weno5(v, CeesC2L_FV, DeesC2L_FV);
  case WENO5_FV_C2R: return __weno5(v, CeesC2R_FV, DeesC2R_FV);
  case WENO5_FV_A2C: return __weno5(v, CeesA2C_FV, DeesA2C_FV);
  case WENO5_FV_C2A: return __weno5(v, CeesC2A_FV, DeesC2A_FV);
  case PLM_C2L: return __plm(v, -1.0);
  case PLM_C2R: return __plm(v, +1.0);
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

static inline double __plm_minmod(double ul, double u0, double ur)
{
  const double a = plm_theta * (u0 - ul);
  const double b =     0.5   * (ur - ul);
  const double c = plm_theta * (ur - u0);
  const double fabc[3] = { fabs(a), fabs(b), fabs(c) };
  return 0.25*fabs(SGN(a)+SGN(b))*(SGN(a)+SGN(c))*min3(fabc);
}
double __plm(const double *v, double sgn)
{
  return v[0] + sgn*0.5*__plm_minmod(v[-1], v[0], v[1]);
}
double __weno5(const double *v, const double c[3][3], const double d[3])
// -----------------------------------------------------------------------------
//
// 
//
// Improvement of the WENO scheme smoothness estimator, Shen & Zha (2010)
// -----------------------------------------------------------------------------
{
  double eps = 1e-6;
  double eps_prime = 1e-6;

  const double vs[3] = {
    c[0][0]*v[+0] + c[0][1]*v[+1] + c[0][2]*v[2],
    c[1][0]*v[-1] + c[1][1]*v[+0] + c[1][2]*v[1],
    c[2][0]*v[-2] + c[2][1]*v[-1] + c[2][2]*v[0],
  };
  double B[3] = { // smoothness indicators
    (13./12.)*SQU(1*v[+0] - 2*v[+1] + 1*v[+2]) +
    ( 1./ 4.)*SQU(3*v[+0] - 4*v[+1] + 1*v[+2]),
    (13./12.)*SQU(1*v[-1] - 2*v[+0] + 1*v[+1]) +
    ( 1./ 4.)*SQU(1*v[-1] - 0*v[+0] - 1*v[+1]),
    (13./12.)*SQU(1*v[-2] - 2*v[-1] + 1*v[+0]) +
    ( 1./ 4.)*SQU(1*v[-2] - 4*v[-1] + 3*v[+0])
  };

  double w[3];
  if (IS_mode == ImprovedBorges08) {
    eps = eps_prime = 1e-14; // Borges uses 1e-40, but has Matlab
    const double tau5 = fabs(B[0] - B[2]);

    // Calculate my weights with new smoothness indicators accoding to Borges
    w[0] = d[0] * (1.0 + (tau5 / (B[0] + eps)));
    w[1] = d[1] * (1.0 + (tau5 / (B[1] + eps)));
    w[2] = d[2] * (1.0 + (tau5 / (B[2] + eps)));
  }
  else if (IS_mode == ImprovedShenZha10) {
    eps = 1e-6;
    eps_prime = 1e-10;
    const double A = shenzha10_A; // [0 (less aggressive) -> ~100 (more aggressive)]
    const double minB = min3(B), maxB = max3(B);
    const double R0 = minB / (maxB + eps_prime);
    B[0] = R0*A*minB + B[0];
    B[1] = R0*A*minB + B[1];
    B[2] = R0*A*minB + B[2];
    w[0] = d[0] / SQU(eps_prime + B[0]);
    w[1] = d[1] / SQU(eps_prime + B[1]);
    w[2] = d[2] / SQU(eps_prime + B[2]);
  }
  else { // Use OriginalJiangShu96
    eps = eps_prime = 1e-6; // recommended value by Jiang and Shu
    w[0] = d[0] / SQU(eps_prime + B[0]);
    w[1] = d[1] / SQU(eps_prime + B[1]);
    w[2] = d[2] / SQU(eps_prime + B[2]);
  }

  const double wtot = w[0] + w[1] + w[2];
  return (w[0]*vs[0] + w[1]*vs[1] + w[2]*vs[2])/wtot;
}

