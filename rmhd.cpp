
/*------------------------------------------------------------------------------
 * FILE: rmhd.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * REFERENCES:
 *
 *------------------------------------------------------------------------------
 */

#include <iostream>
#include <cstdio>
#include <cmath>
#include "rmhd.hpp"
#include "rmhd-c2p.h"
#include "eos.hpp"
#include "quartic.hpp"
#include "nrsolver.hpp"
#include "logging.hpp"

#define EOS (Mara->eos)

typedef AdiabaticIdealRmhd Rmhd;



std::vector<std::string> Rmhd::GetPrimNames() const
{
  std::string vars[8] = { "rho", "pre",
                          "vx" , "vy", "vz",
                          "Bx" , "By", "Bz" };
  return std::vector<std::string>(vars, vars+8);
}
std::string Rmhd::PrintPrim(const double *P) const
{
  char str[512];
  sprintf(str, "P = { %+6.4e, %+6.4e, %+6.4e, %+6.4e, %+6.4e, %+6.4e, %+6.4e, %+6.4e }",
          P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7]);
  return std::string(str);
}
std::string Rmhd::PrintCons(const double *U) const
{
  char str[512];
  sprintf(str, "U = { %+6.4e, %+6.4e, %+6.4e, %+6.4e, %+6.4e, %+6.4e, %+6.4e, %+6.4e }",
          U[0], U[1], U[2], U[3], U[4], U[5], U[6], U[7]);
  return std::string(str);
}

int Rmhd::PrimCheck(const double *P) const
{
  return rmhd_c2p_check_prim(P);
}
int Rmhd::ConsCheck(const double *U) const
{
  return rmhd_c2p_check_cons(U);
}

int Rmhd::ConsToPrim(const double *U, double *P) const
{
  int error = 1;

  // This piece of code drives cons to prim inversions for an arbitrary equation
  // of state.
  // ---------------------------------------------------------------------------
  if (typeid(*EOS) != typeid(AdiabaticEos)) {

    rmhd_c2p_eos_set_eos(EOS);
    rmhd_c2p_eos_new_state(U);

    //    std::cout << PrintPrim(P) << std::endl;
    //    std::cout << PrintCons(U) << std::endl;

    if (error) {
      rmhd_c2p_eos_set_starting_prim(P);
      error = rmhd_c2p_eos_solve_duffell3d(P);
    }
    if (error) {

      // NOTE: disregarding further c2p trials for debugging purposes
      return error;
      // ------------------------------------------------------------

      rmhd_c2p_eos_set_starting_prim(P);
      error = rmhd_c2p_eos_solve_noble2dzt(P);
    }
    if (error) {
      rmhd_c2p_eos_estimate_from_cons();
      error = rmhd_c2p_eos_solve_noble2dzt(P);
    }

    return error;
  }

  // This piece of code drives cons to prim inversions for a gamma-law equation
  // of state.
  // ---------------------------------------------------------------------------
  rmhd_c2p_set_gamma(Mara->GetEos<AdiabaticEos>().Gamma);
  rmhd_c2p_new_state(U);

  if (error) {
    rmhd_c2p_set_starting_prim(P);
    error = rmhd_c2p_solve_anton2dzw(P);
  }
  if (error) {
    rmhd_c2p_estimate_from_cons();
    error = rmhd_c2p_solve_anton2dzw(P);
  }
  if (error) {
    rmhd_c2p_set_starting_prim(P);
    error = rmhd_c2p_solve_noble1dw(P);
  }
  if (error) {
    rmhd_c2p_estimate_from_cons();
    error = rmhd_c2p_solve_noble1dw(P);
  }

  return error;
}
int Rmhd::PrimToCons(const double *P, double *U) const
{
  const double V2   =   P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz];
  const double B2   =   P[Bx]*P[Bx] + P[By]*P[By] + P[Bz]*P[Bz];
  const double Bv   =   P[Bx]*P[vx] + P[By]*P[vy] + P[Bz]*P[vz];
  const double W2   =   1.0 / (1.0 - V2);
  const double W    =   sqrt(W2);
  const double b0   =   W*Bv;
  const double b2   =   (B2 + b0*b0) / W2;
  const double bx   =  (P[Bx] + b0 * W*P[vx]) / W;
  const double by   =  (P[By] + b0 * W*P[vy]) / W;
  const double bz   =  (P[Bz] + b0 * W*P[vz]) / W;
  const double e    =   EOS->Internal(P[rho], EOS->Temperature_p(P[rho], P[pre]))/P[rho];
  const double e_   =   e      + 0.5 * b2 / P[rho];
  const double p_   =   P[pre] + 0.5 * b2;
  const double h_   =   1 + e_ + p_ / P[rho];

  U[ddd] = P[rho] * W;
  U[tau] = P[rho] * h_ * W2 - p_    - b0*b0 - U[ddd];
  U[Sx ] = P[rho] * h_ * W2 * P[vx] - b0*bx;
  U[Sy ] = P[rho] * h_ * W2 * P[vy] - b0*by;
  U[Sz ] = P[rho] * h_ * W2 * P[vz] - b0*bz;
  U[Bx ] = P[Bx ];
  U[By ] = P[By ];
  U[Bz ] = P[Bz ];

  if (V2 >= 1.0) {
    return 1;
  }

  return 0;
}
void Rmhd::FluxAndEigenvalues(const double *U,
                              const double *P, double *F,
                              double *ap, double *am, int dimension) const
{
  const double V2   =   P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz];
  const double B2   =   P[Bx]*P[Bx] + P[By]*P[By] + P[Bz]*P[Bz];
  const double Bv   =   P[Bx]*P[vx] + P[By]*P[vy] + P[Bz]*P[vz];
  const double W2   =   1.0 / (1.0 - V2);
  const double W    =   sqrt(W2);
  const double b0   =   W*Bv;
  const double b2   =  (B2 + b0*b0) / W2;
  const double bx   =  (P[Bx] + b0 * W*P[vx]) / W;
  const double by   =  (P[By] + b0 * W*P[vy]) / W;
  const double bz   =  (P[Bz] + b0 * W*P[vz]) / W;
  const double e    =   EOS->Internal(P[rho], EOS->Temperature_p(P[rho], P[pre]))/P[rho];
  const double p_   =   P[pre] + 0.5 * b2;
  const double h    =   1 + e + P[pre]/P[rho];

  switch (dimension) {
  case 1:
    F[ddd] = U[ddd] * P[vx];
    F[tau] = U[tau] * P[vx] - b0*P[Bx] / W + p_*P[vx];
    F[Sx ] = U[Sx]  * P[vx] - bx*P[Bx] / W + p_;
    F[Sy ] = U[Sy]  * P[vx] - by*P[Bx] / W;
    F[Sz ] = U[Sz]  * P[vx] - bz*P[Bx] / W;
    F[Bx ] = 0.0;
    F[By ] = P[vx]*P[By] - P[vy]*P[Bx];
    F[Bz ] = P[vx]*P[Bz] - P[vz]*P[Bx];
    break;
  case 2:
    F[ddd] = U[ddd] * P[vy];
    F[tau] = U[tau] * P[vy] - b0*P[By] / W + p_*P[vy];
    F[Sx ] = U[Sx]  * P[vy] - bx*P[By] / W;
    F[Sy ] = U[Sy]  * P[vy] - by*P[By] / W + p_;
    F[Sz ] = U[Sz]  * P[vy] - bz*P[By] / W;
    F[Bx ] = P[vy]*P[Bx] - P[vx]*P[By];
    F[By ] = 0.0;
    F[Bz ] = P[vy]*P[Bz] - P[vz]*P[By];
    break;
  case 3:
    F[ddd] = U[ddd] * P[vz];
    F[tau] = U[tau] * P[vz] - b0*P[Bz] / W + p_*P[vz];
    F[Sx ] = U[Sx]  * P[vz] - bx*P[Bz] / W;
    F[Sy ] = U[Sy]  * P[vz] - by*P[Bz] / W;
    F[Sz ] = U[Sz]  * P[vz] - bz*P[Bz] / W + p_;
    F[Bx ] = P[vz]*P[Bx] - P[vx]*P[Bz];
    F[By ] = P[vz]*P[By] - P[vy]*P[Bz];
    F[Bz ] = 0.0;
    break;
  }

  /* Begin eigenvalue calculation
   * ---------------------------------------------------------------------------
   *
   */
  if (ap == NULL || am == NULL) return;  // User may skip eigenvalue calculation

  double vi=0, bi=0;
  switch (dimension) {
  case 1:
    vi = P[vx]; bi = bx;
    break;
  case 2:
    vi = P[vy]; bi = by;
    break;
  case 3:
    vi = P[vz]; bi = bz;
    break;
  }

  const double cs2  =  EOS->SoundSpeed2Sr(P[rho], EOS->Temperature_p(P[rho], P[pre]));
  const double W4   =  W2*W2;
  const double v2   =  vi*vi;
  const double v3   =  vi*v2;
  const double v4   =  vi*v3;

  const double K  =    P[rho]*h * (1.0/cs2-1) * W4;
  const double L  =  -(P[rho]*h +   b2/cs2)   * W2;

  const double A4 =    K    - L          -   b0*b0;
  const double A3 = -4*K*vi + L*vi*2     + 2*b0*bi;
  const double A2 =  6*K*v2 + L*(1.0-v2) +   b0*b0 - bi*bi;
  const double A1 = -4*K*v3 - L*vi*2     - 2*b0*bi;
  const double A0 =    K*v4 + L*v2       +   bi*bi;

  QuarticEquation quartic(A4,A3,A2,A1,A0);

  double roots[4];
  int nr = quartic.Solve(roots);

  *ap = -1.0;
  *am =  1.0;

  for (int i=0; i<nr; ++i) {
    if (roots[i] > *ap) *ap = roots[i];
    if (roots[i] < *am) *am = roots[i];
  }

  if (fabs(*ap)>1.0 || fabs(*am)>1.0) {
    *ap =  1.0;
    *am = -1.0;
  }
}
void Rmhd::ConstrainedTransport2d(double *Fx, double *Fy,
                                  int stride[4]) const
{
  const int N = stride[0]/8;
  double *FxBy = (double*) malloc(N*sizeof(double));
  double *FyBx = (double*) malloc(N*sizeof(double));

  double *F, *G;
  int i;

  const int sx=stride[1],sy=stride[2];
  for (i=sx; i<stride[0]-sx; i+=8) {
    F = &Fx[By+i];
    G = &Fy[Bx+i];

    FxBy[i/8] = (2*F[0]+F[sy]+F[-sy]-G[0]-G[sx]-G[-sy]-G[ sx-sy])*0.125;
    FyBx[i/8] = (2*G[0]+G[sx]+G[-sx]-F[0]-F[sy]-F[-sx]-F[-sx+sy])*0.125;
  }
  for (i=0; i<stride[0]; i+=8) {
    Fx[i+Bx] = 0.0;       Fx[i+By] = FxBy[i/8];
    Fy[i+Bx] = FyBx[i/8]; Fy[i+By] = 0.0;
  }

  free(FxBy);
  free(FyBx);
}
void Rmhd::ConstrainedTransport3d(double *Fx, double *Fy, double *Fz,
                                  int stride[4]) const
{
  const int N = stride[0]/8;
  double *FxBy = (double*) malloc(N*sizeof(double));
  double *FxBz = (double*) malloc(N*sizeof(double));

  double *FyBz = (double*) malloc(N*sizeof(double));
  double *FyBx = (double*) malloc(N*sizeof(double));

  double *FzBx = (double*) malloc(N*sizeof(double));
  double *FzBy = (double*) malloc(N*sizeof(double));

  double *F, *G, *H;
  int i;

  const int sx=stride[1],sy=stride[2],sz=stride[3];
  for (i=sx; i<stride[0]-sx; i+=8) {
    F = &Fx[By+i];
    G = &Fy[Bx+i];

    FxBy[i/8] = (2*F[0]+F[sy]+F[-sy]-G[0]-G[sx]-G[-sy]-G[ sx-sy])*0.125;
    FyBx[i/8] = (2*G[0]+G[sx]+G[-sx]-F[0]-F[sy]-F[-sx]-F[-sx+sy])*0.125;

    G = &Fy[Bz+i];
    H = &Fz[By+i];

    FyBz[i/8] = (2*G[0]+G[sz]+G[-sz]-H[0]-H[sy]-H[-sz]-H[ sy-sz])*0.125;
    FzBy[i/8] = (2*H[0]+H[sy]+H[-sy]-G[0]-G[sz]-G[-sy]-G[-sy+sz])*0.125;

    H = &Fz[Bx+i];
    F = &Fx[Bz+i];

    FzBx[i/8] = (2*H[0]+H[sx]+H[-sx]-F[0]-F[sz]-F[-sx]-F[ sz-sx])*0.125;
    FxBz[i/8] = (2*F[0]+F[sz]+F[-sz]-H[0]-H[sx]-H[-sz]-H[-sz+sx])*0.125;
  }
  for (i=0; i<stride[0]; i+=8) {
    Fx[i+Bx] = 0.0;        Fx[i+By] = FxBy[i/8];  Fx[i+Bz] = FxBz[i/8];
    Fy[i+Bx] = FyBx[i/8];  Fy[i+By] = 0.0;        Fy[i+Bz] = FyBz[i/8];
    Fz[i+Bx] = FzBx[i/8];  Fz[i+By] = FzBy[i/8];  Fz[i+Bz] = 0.0;
  }

  free(FxBy);  free(FyBz);  free(FzBx);
  free(FxBz);  free(FyBx);  free(FzBy);
}
