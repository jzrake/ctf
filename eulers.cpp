
/*------------------------------------------------------------------------------
 * FILE: euler.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * REFERENCES:
 *
 *------------------------------------------------------------------------------
 */

#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "eulers.hpp"
#include "eos.hpp"
#include "logging.hpp"
#include "matrix.h"



typedef AdiabaticIdealEulers Eulers;


std::vector<std::string> Eulers::GetPrimNames() const
{
  std::string vars[5] = { "rho", "pre",
                          "vx" , "vy", "vz" };
  return std::vector<std::string>(vars, vars+5);
}
std::string Eulers::PrintPrim(const double *P) const
{
  char str[512];
  sprintf(str, "P = { %+6.4e %+6.4e %+6.4e %+6.4e %+6.4e }",
          P[0], P[1], P[2], P[3], P[4]);
  return std::string(str);
}
std::string Eulers::PrintCons(const double *U) const
{
  char str[512];
  sprintf(str, "U = { %+6.4e %+6.4e %+6.4e %+6.4e %+6.4e }",
          U[0], U[1], U[2], U[3], U[4]);
  return std::string(str);
}
int Eulers::ConsToPrim(const double *U, double *P) const
{
  const double gm1 = Mara->GetEos<AdiabaticEos>().Gamma - 1.0;

  P[rho] = U[rho];
  P[pre] =(U[nrg] - 0.5*(U[px]*U[px] + U[py]*U[py] + U[pz]*U[pz])/U[rho])*gm1;
  P[vx ] = U[px ] / U[rho];
  P[vy ] = U[py ] / U[rho];
  P[vz ] = U[pz ] / U[rho];

  if (P[rho] < 0.0) {
    DebugLog.Error(__FUNCTION__) << "Got negative density." << std::endl;
    DebugLog.Error() << PrintPrim(P) << std::endl
                     << PrintCons(U) << std::endl;
    return 1;
  }
  if (P[pre] < 0.0) {
    DebugLog.Error(__FUNCTION__) << "Got negative pressure." << std::endl;
    DebugLog.Error() << PrintPrim(P) << std::endl
                     << PrintCons(U) << std::endl;
    return 1;
  }
  if (U[nrg] < 0.0) {
    DebugLog.Error(__FUNCTION__) << "Got negative energy." << std::endl;
    DebugLog.Error() << PrintPrim(P) << std::endl
                     << PrintCons(U) << std::endl;
    return 1;
  }

  return 0;
}
int Eulers::PrimToCons(const double *P, double *U) const
{
  const double gm1 = Mara->GetEos<AdiabaticEos>().Gamma - 1.0;

  U[rho] = P[rho];
  U[px]  = P[rho] * P[vx];
  U[py]  = P[rho] * P[vy];
  U[pz]  = P[rho] * P[vz];
  U[nrg] = P[rho] * 0.5*(P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz]) + P[pre]/gm1;

  return 0;
}
void Eulers::FluxAndEigenvalues(const double *U,
                                const double *P, double *F,
                                double *ap, double *am, int dimension) const
{
  switch (dimension) {
  case 1:
    F[rho]  =  U[rho] * P[vx];
    F[nrg]  = (U[nrg] + P[pre])*P[vx];
    F[px]   =  U[px]  * P[vx] + P[pre];
    F[py]   =  U[py]  * P[vx];
    F[pz]   =  U[pz]  * P[vx];
    break;
  case 2:
    F[rho]  =  U[rho] * P[vy];
    F[nrg]  = (U[nrg] + P[pre])*P[vy];
    F[px]   =  U[px]  * P[vy];
    F[py]   =  U[py]  * P[vy] + P[pre];
    F[pz]   =  U[pz]  * P[vy];
    break;
  case 3:
    F[rho]  =  U[rho] * P[vz];
    F[nrg]  = (U[nrg] + P[pre])*P[vz];
    F[px]   =  U[px]  * P[vz];
    F[py]   =  U[py]  * P[vz];
    F[pz]   =  U[pz]  * P[vz] + P[pre];
    break;
  }
  if (ap == 0 || am == 0) return; // User may skip eigenvalue calculation

  const double gm = Mara->GetEos<AdiabaticEos>().Gamma;

  switch (dimension) {
  case 1:
    *ap = P[vx] + sqrt(gm*P[pre] / U[rho]);
    *am = P[vx] - sqrt(gm*P[pre] / U[rho]);
    break;
  case 2:
    *ap = P[vy] + sqrt(gm*P[pre] / U[rho]);
    *am = P[vy] - sqrt(gm*P[pre] / U[rho]);
    break;
  case 3:
    *ap = P[vz] + sqrt(gm*P[pre] / U[rho]);
    *am = P[vz] - sqrt(gm*P[pre] / U[rho]);
    break;
  }
}

#define EIGEN2

#ifdef EIGEN1
void Eulers::Eigensystem(const double *U, const double *P,
                         double *L, double *R, double *lam, int dim) const
{
  const double gm = Mara->GetEos<AdiabaticEos>().Gamma;
  const double g1 = gm - 1.0;
  const double u = P[vx];
  const double v = P[vy];
  const double w = P[vz];
  const double ek = 0.5*(u*u + v*v + w*w);
  const double a2 = gm * P[pre] / P[rho];
  const double a  = sqrt(a2);
  const double h0 = a2/g1 + ek;

  //  const double n0 = sqrt(nhat[0]*nhat[0] + nhat[1]*nhat[1] + nhat[2]*nhat[2]);
  const double nx = dim==1;//nhat[0]/n0;
  const double ny = dim==2;//nhat[1]/n0;
  const double nz = dim==3;//nhat[2]/n0;
  const double vn = u*nx + v*ny + w*nz;


  if (fabs(nx) >= fabs(ny) && fabs(nx) >= fabs(nz)) {
    const double R_[5][5] =
      { {         1,  1,         1,           0,   0         },
        { h0 - a*vn, ek, h0 + a*vn, u*ny - v*nx, w*nx - u*nz },
        {  u - a*nx,  u,  u + a*nx,          ny, -nz         },
        {  v - a*ny,  v,  v + a*ny,         -nx,   0         },
        {  w - a*nz,  w,  w + a*nz,           0,  nx         } };

    const double L_[5][5] =
      { { (g1*ek + a*vn)/(2*a2), g1/(2*a2), -(g1*u + a*nx)/(2*a2), -(g1*v + a*ny)/(2*a2), -(g1*w + a*nz)/(2*a2) },
        { (a2 - g1*ek) / a2, -g1/a2, g1*u/a2, g1*v/a2, g1*w/a2 },
        { (g1*ek - a*vn)/(2*a2), g1/(2*a2), -(g1*u - a*nx)/(2*a2), -(g1*v - a*ny)/(2*a2), -(g1*w - a*nz)/(2*a2) },
        { (v - vn*ny) / nx, 0.0, ny, (ny*ny - 1.0) / nx, ny*nz/nx },
        {-(w - vn*nz) / nx, 0.0,-nz,-ny*nz/nx,-(nz*nz - 1.0) / nx } };

    std::memcpy(R, R_, 5*5*sizeof(double));
    std::memcpy(L, L_, 5*5*sizeof(double));
  }

  else if (fabs(ny) >= fabs(nz) && fabs(ny) >= fabs(nx)) {
    const double R_[5][5] =
      { {         1,  1,         1,           0,   0         },
        { h0 - a*vn, ek, h0 + a*vn, u*ny - v*nx, v*nz - w*ny },
        {  u - a*nx,  u,  u + a*nx,          ny,   0         },
        {  v - a*ny,  v,  v + a*ny,         -nx,  nz         },
        {  w - a*nz,  w,  w + a*nz,           0, -ny         } };

    const double L_[5][5] =
      { { (g1*ek + a*vn)/(2*a2), g1/(2*a2), -(g1*u + a*nx)/(2*a2), -(g1*v + a*ny)/(2*a2), -(g1*w + a*nz)/(2*a2) },
        { (a2 - g1*ek) / a2, -g1/a2, g1*u/a2, g1*v/a2, g1*w/a2 },
        { (g1*ek - a*vn)/(2*a2), g1/(2*a2), -(g1*u - a*nx)/(2*a2), -(g1*v - a*ny)/(2*a2), -(g1*w - a*nz)/(2*a2) },
        {-(u - vn*nx) / ny, 0.0,-(nx*nx - 1.0) / ny,-nx,-nx*nz/ny },
        { (w - vn*nz) / ny, 0.0, nx*nz/ny, nz, (nz*nz - 1.0) / ny } };

    std::memcpy(R, R_, 5*5*sizeof(double));
    std::memcpy(L, L_, 5*5*sizeof(double));
  }

  else if (fabs(nz) >= fabs(nx) && fabs(nz) >= fabs(ny)) {
    const double R_[5][5] =
      { {         1,  1,         1,           0,   0         },
        { h0 - a*vn, ek, h0 + a*vn, w*nx - u*nz, v*nz - w*ny },
        {  u - a*nx,  u,  u + a*nx,         -nz,   0         },
        {  v - a*ny,  v,  v + a*ny,           0,  nz         },
        {  w - a*nz,  w,  w + a*nz,          nx, -ny         } };

    const double L_[5][5] =
      { { (g1*ek + a*vn)/(2*a2), g1/(2*a2), -(g1*u + a*nx)/(2*a2), -(g1*v + a*ny)/(2*a2), -(g1*w + a*nz)/(2*a2) },
        { (a2 - g1*ek) / a2, -g1/a2, g1*u/a2, g1*v/a2, g1*w/a2 },
        { (g1*ek - a*vn)/(2*a2), g1/(2*a2), -(g1*u - a*nx)/(2*a2), -(g1*v - a*ny)/(2*a2), -(g1*w - a*nz)/(2*a2) },
        { (u - vn*nx) / nz, 0.0, (nx*nx - 1.0) / nz, nx*ny/nz, nx },
        {-(v - vn*ny) / nz, 0.0,-nx*ny/nz,-(ny*ny - 1.0) / nz,-ny } };

    std::memcpy(R, R_, 5*5*sizeof(double));
    std::memcpy(L, L_, 5*5*sizeof(double));
  }

  lam[0] = vn-a;
  lam[1] = vn;
  lam[2] = vn+a;
  lam[3] = vn;
  lam[4] = vn;
}
#endif

#ifdef EIGEN2
void Eulers::Eigensystem(const double *U, const double *P_,
                         double *L, double *R, double *lam, int dim) const
// -----------------------------------------------------------------------------
// For 1d only
// -----------------------------------------------------------------------------
{
  const double T[3][5][5] =
  // Tx
    {{{1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 1, 0, 0},
      {0, 0, 0, 1, 0},
      {0, 0, 0, 0, 1}},
  // Ty
     {{1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 0, 1, 0},
      {0, 0,-1, 0, 0},
      {0, 0, 0, 0, 1}},
  // Tz
     {{1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 0, 0, 1},
      {0, 0, 0, 1, 0},
      {0, 0,-1, 0, 0}}};

  const double S[3][5][5] = // T^{-1}
  // Sx
    {{{1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 1, 0, 0},
      {0, 0, 0, 1, 0},
      {0, 0, 0, 0, 1}},
  // Sy
     {{1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 0,-1, 0},
      {0, 0, 1, 0, 0},
      {0, 0, 0, 0, 1}},
  // Sz
     {{1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 0, 0,-1},
      {0, 0, 0, 1, 0},
      {0, 0, 1, 0, 0}}};

  double P[5];
  matrix_vector_product(T[dim-1][0], P_, P, 5, 5);

  const double gm = Mara->GetEos<AdiabaticEos>().Gamma;
  const double gm1 = gm - 1.0;
  const double u = P[vx];
  const double v = P[vy];
  const double w = P[vz];
  const double V2 = u*u + v*v + w*w;
  const double a = sqrt(gm * P[pre] / P[rho]);
  const double H = (U[nrg] + P[pre]) / P[rho];

  // --------------------------------------------------------------------------
  // Toro Equation 3.82 (rows are permuted to deal with Mara's convention on
  // the conserved quantities)
  // --------------------------------------------------------------------------
  const double RR[5][5] =
    {{       1,       1,      1,      0,      0 },  // rho
     { H - u*a, H + u*a, 0.5*V2,      v,      w },  // nrg
     {     u-a,     u+a,      u,      0,      0 },  // px
     {       v,       v,      v,      1,      0 },  // py
     {       w,       w,      w,      0,      1 }}; // pz
  // --------------------------------------------------------------------------
  // Toro Equation 3.83 up to (gam - 1) / (2*a^2)
  // --------------------------------------------------------------------------
  const double LL[5][5] =
    {{    H + (a/gm1)*(u-a),   1, -(u+a/gm1),        -v,        -w },
     {    H - (a/gm1)*(u+a),   1, -(u-a/gm1),        -v,        -w },
     { -2*H + (4/gm1)*(a*a),  -2,        2*u,       2*v,       2*w },
     {         -2*v*a*a/gm1,   0,          0, 2*a*a/gm1,         0 },
     {         -2*w*a*a/gm1,   0,          0,         0, 2*a*a/gm1 }};
  // --------------------------------------------------------------------------
  //                    rho, nrg,         px,        py,        pz
  // --------------------------------------------------------------------------

  matrix_matrix_product(S[dim-1][0], RR[0], R, 5, 5, 5);
  matrix_matrix_product(LL[0], T[dim-1][0], L, 5, 5, 5);

  // Replace the term in eqn 3.83 : (gam - 1) / (2*a^2)
  // ---------------------------------------------------------------------------
  const double norm = gm1 / (2*a*a);
  for (int i=0; i<25; ++i) {
    L[i] *= norm;
  }

  lam[0] = u-a;
  lam[1] = u+a;
  lam[2] = u;
  lam[3] = u;
  lam[4] = u;
}
#endif
