
/*------------------------------------------------------------------------------
 * FILE: rmhd.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * REFERENCES:
 *
 *------------------------------------------------------------------------------
 */

#include <cstdio>
#include <cmath>
#include <cstring>
#include "srhd.hpp"
#include "matrix.h"
#include "logging.hpp"


typedef AdiabaticIdealSrhd Srhd;
static const double bigW = 1e12;

Srhd::AdiabaticIdealSrhd(double AdiabaticGamma)
  : AdiabaticGamma(AdiabaticGamma) { }

std::vector<std::string> Srhd::GetPrimNames() const
{
  std::string vars[5] = { "rho", "pre",
                          "vx" , "vy", "vz" };
  return std::vector<std::string>(vars, vars+5);
}
std::string Srhd::PrintPrim(const double *P) const
{
  char str[512];
  sprintf(str, "P = { %+6.4e, %+6.4e, %+6.4e, %+6.4e, %+6.4e }",
          P[0], P[1], P[2], P[3], P[4]);
  return std::string(str);
}
std::string Srhd::PrintCons(const double *U) const
{
  char str[512];
  sprintf(str, "U = { %+6.4e, %+6.4e, %+6.4e, %+6.4e, %+6.4e }",
          U[0], U[1], U[2], U[3], U[4]);
  return std::string(str);
}
int Srhd::ConsToPrim(const double *U, double *P) const
{
  if (U[ddd] < 0.0 || U[tau] < 0.0) {
    DebugLog.Warning(__FUNCTION__) << "Got negative D or Tau." << std::endl
                                   << PrintCons(U) << std::endl;
    return 1;
  }

  static const double ERROR_TOLR = 1e-6;
  static const int NEWTON_MAX_ITER = 50;

  const double D     = U[ddd];
  const double Tau   = U[tau];
  const double S2    = U[Sx]*U[Sx] + U[Sy]*U[Sy] + U[Sz]*U[Sz];

  int soln_found     = 0;
  int n_iter         = 0;

  double f,g,W_soln=1.0,p=P[pre];

  while (!soln_found) {

    const double v2  = S2 / pow(Tau + D + p, 2);
    const double W2  = 1.0 / (1.0 - v2);
    const double W   = sqrt(W2);
    const double e   = (Tau + D*(1.0 - W) + p*(1.0 - W2)) / (D*W);
    const double Rho = D / W;
    const double h   = 1.0 + e + p/Rho;
    const double cs2 = AdiabaticGamma * p / (Rho * h);

    f = Rho * e * (AdiabaticGamma - 1.0) - p;
    g = v2*cs2 - 1.0;

    p -= f/g;

    if (fabs(f) < ERROR_TOLR) {
      W_soln = W;
      soln_found = 1;
    }
    if (n_iter++ == NEWTON_MAX_ITER) {
      DebugLog.Warning(__FUNCTION__) << "ConsToPrim ran out of iterations." << std::endl
                                     << PrintPrim(P) << std::endl
                                     << PrintCons(U) << std::endl;
      return 1;
    }
  }

  P[rho] = D/W_soln;
  P[pre] = p;
  P[vx]  = U[Sx] / (Tau + D + p);
  P[vy]  = U[Sy] / (Tau + D + p);
  P[vz]  = U[Sz] / (Tau + D + p);

  if (P[pre] < 0.0) {
    DebugLog.Warning(__FUNCTION__) << "Got negative pressure." << std::endl
                                   << PrintPrim(P) << std::endl
                                   << PrintCons(U) << std::endl;
    return 1;
  }
  if (P[rho] < 0.0) {
    DebugLog.Warning(__FUNCTION__) << "Got negative density." << std::endl
                                   << PrintPrim(P) << std::endl
                                   << PrintCons(U) << std::endl;
    return 1;
  }
  if (std::isnan(W_soln) || W_soln > bigW) {
    DebugLog.Warning(__FUNCTION__) << "Got superluminal velocity." << std::endl
                                   << PrintPrim(P) << std::endl
                                   << PrintCons(U) << std::endl;
    return 1;
  }

  return 0;
}
int Srhd::PrimToCons(const double *P, double *U) const
{
  const double V2   =   P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz];
  if (V2 >= 1.0) {
    DebugLog.Warning(__FUNCTION__) << "Got superluminal velocity, V2 = " << V2 << std::endl
                                   << PrintPrim(P) << std::endl
                                   << PrintCons(U) << std::endl;
    return 1;
  }

  const double W2   =   1.0 / (1.0 - V2);
  const double W    =   sqrt(W2);
  const double e    =   P[pre] / (P[rho] * (AdiabaticGamma - 1.0));
  const double h    =   1.0 + e + P[pre]/P[rho];

  U[ddd] = P[rho]*W;
  U[tau] = P[rho]*h*W2 - P[pre] - U[ddd];
  U[Sx]  = P[rho]*h*W2*P[vx];
  U[Sy]  = P[rho]*h*W2*P[vy];
  U[Sz]  = P[rho]*h*W2*P[vz];

  return 0;
}
void Srhd::FluxAndEigenvalues(const double *U,
                              const double *P, double *F,
                              double *ap, double *am, int dimension) const
{
  switch (dimension) {
  case 1:
    F[ddd] = U[ddd] * P[vx ];
    F[tau] = U[Sx ] - U[ddd] * P[vx ];
    F[Sx ] = U[Sx ] * P[vx ] + P[pre];
    F[Sy ] = U[Sy ] * P[vx ];
    F[Sz ] = U[Sz ] * P[vx ];
    break;
  case 2:
    F[ddd] = U[ddd] * P[vy ];
    F[tau] = U[Sy ] - U[ddd] * P[vy ];
    F[Sx ] = U[Sx ] * P[vy ];
    F[Sy ] = U[Sy ] * P[vy ] + P[pre];
    F[Sz ] = U[Sz ] * P[vy ];
    break;
  case 3:
    F[ddd] = U[ddd] * P[vz ];
    F[tau] = U[Sz ] - U[ddd] * P[vz ];
    F[Sx ] = U[Sx ] * P[vz ];
    F[Sy ] = U[Sy ] * P[vz ];
    F[Sz ] = U[Sz ] * P[vz ] + P[pre];
    break;
  }

  if (ap == 0 || am == 0) return; // User may skip eigenvalue calculation

  const double vx2 = P[vx]*P[vx];
  const double vy2 = P[vy]*P[vy];
  const double vz2 = P[vz]*P[vz];
  const double   e = P[pre] / (P[rho] * (AdiabaticGamma - 1.0));
  const double cs2 = AdiabaticGamma * P[pre] / (P[pre] + P[rho] + P[rho]*e);
  const double v2  = vx2 + vy2 + vz2;

  switch (dimension) {
  case 1:
    *ap = (P[vx]*(1-cs2) + sqrt(cs2*(1-v2)*(1-v2*cs2-vx2*(1-cs2))))/(1-v2*cs2);
    *am = (P[vx]*(1-cs2) - sqrt(cs2*(1-v2)*(1-v2*cs2-vx2*(1-cs2))))/(1-v2*cs2);
    break;
  case 2:
    *ap = (P[vy]*(1-cs2) + sqrt(cs2*(1-v2)*(1-v2*cs2-vy2*(1-cs2))))/(1-v2*cs2);
    *am = (P[vy]*(1-cs2) - sqrt(cs2*(1-v2)*(1-v2*cs2-vy2*(1-cs2))))/(1-v2*cs2);
    break;
  case 3:
    *ap = (P[vz]*(1-cs2) + sqrt(cs2*(1-v2)*(1-v2*cs2-vz2*(1-cs2))))/(1-v2*cs2);
    *am = (P[vz]*(1-cs2) - sqrt(cs2*(1-v2)*(1-v2*cs2-vz2*(1-cs2))))/(1-v2*cs2);
    break;
  }

  if (fabs(*ap)>1.0 || fabs(*am)>1.0) {
    DebugLog.Info(__FUNCTION__)
      << "superluminal eigenvalues: " << *ap << " " << *am
      << " ... resetting to +/- 1.0 " << std::endl;
    *ap =  1.0;
    *am = -1.0;
  }

  return;
}


void Srhd::Eigensystem(const double *U, const double *P_,
                       double *L, double *R, double *lam, int dim) const
// -----------------------------------------------------------------------------
//
// Authors: Jonathan Zrake, Bez Laderman: NYU CCPP
//
// Date: May 7th, 2012
//
// This piece of code implements the left and right eigenvectors of the ideal
// special relativistic hydrodynamics equations. The formulas are an exact
// translation of those given in the literature:
//
// R. Donat, J.A. Font, J.M. Ibanez, & A. Marquina
// JCP, 1998, 146, 58
//
// http://www.sciencedirect.com/science/article/pii/S0021999198959551
//
//
// Having these eigenvectors in a hydrodynamics code is good. They can be used
// for any scheme which requires flux splitting with characteristic
// decomposition, such as high order ENO or WENO schemes.
//
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
      {0, 0, 0,-1, 0},
      {0, 0, 1, 0, 0},
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
      {0, 0, 0, 1, 0},
      {0, 0,-1, 0, 0},
      {0, 0, 0, 0, 1}},
  // Sz
     {{1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 0, 0,-1},
      {0, 0, 0, 1, 0},
      {0, 0, 1, 0, 0}}};

  double P[5];
  matrix_vector_product(T[dim-1][0], P_, P, 5, 5);

  const double D = P[ddd]; // rest mass density
  const double p = P[pre]; // pressure
  const double u = P[vx];  // vx (three velocity)
  const double v = P[vy];  // vy
  const double w = P[vz];  // vz

  const double sie = (p/D) / (AdiabaticGamma - 1); // specific internal energy
  const double h = 1 + sie + p/D;                  // specific enthalpy
  const double cs2 = AdiabaticGamma * p / (D*h);   // sound speed squared
  const double V2 = u*u + v*v + w*w;
  const double W = 1.0 / sqrt(1 - V2);             // Lorentz factor
  const double W2 = W*W;
  const double K = h;                              // for gamma-law only, K = h
  const double hW = h*W;

  // equations (14) and (15)
  const double lp = (u*(1-cs2) + sqrt(cs2*(1-V2)*(1-V2*cs2-u*u*(1-cs2))))/(1-V2*cs2);
  const double lm = (u*(1-cs2) - sqrt(cs2*(1-V2)*(1-V2*cs2-u*u*(1-cs2))))/(1-V2*cs2);

  const double Ap = (1 - u*u) / (1 - u*lp);
  const double Am = (1 - u*u) / (1 - u*lm);

  // Equations (17) through (20)
  // ---------------------------------------------------------------------------
  const double RT[5][5] =
    {{1, hW*Ap - 1, hW*Ap*lp, hW*v, hW*w},
     {1, hW*Am - 1, hW*Am*lm, hW*v, hW*w},
     {K/hW, 1-K/hW, u, v, w},
     {W*v, 2*h*W2*v - W*v, 2*h*W2*u*v, h*(1+2*W2*v*v), 2*h*W2*v*w},
     {W*w, 2*h*W2*w - W*w, 2*h*W2*u*w, 2*h*W2*v*w, h*(1+2*W2*w*w)}};

  double RR[5][5];
  for (int n=0; n<5; ++n) {
    for (int m=0; m<5; ++m) {
      RR[n][m] = RT[m][n];
    }
  }

  // NOTES
  // ---------------------------------------------------------------------------
  // (1) Font describes the columns if the left eigen vector matrix
  // horizontally, which is how they are written below. So we take the
  // transpose at the end of the day.
  //
  // (2) Font's notation uses L_{-/+} for the last left eigenvectors, but that
  // naming is weird, since each L_+ contains lm and Am and vice-versa.
  // ---------------------------------------------------------------------------

  const double Delta = h*h*h*W*(K-1)*(1-u*u)*(Ap*lp - Am*lm); // equation (21)
  const double a = W / (K-1);
  const double b = 1 / (h*(1 - u*u));
  const double c = 1 / (h*(1 - u*u));
  const double d = -h*h / Delta;
  const double e = +h*h / Delta;

  const double LL[5][5] =
    {{d*(hW*Am*(u-lm) - u - W2*(V2 - u*u)*(2*K - 1)*(u - Am*lm) + K*Am*lm),
      d*(-u - W2*(V2 - u*u)*(2*K - 1)*(u - Am*lm) + K*Am*lm),
      d*(1 + W2*(V2 - u*u)*(2*K - 1)*(1 - Am) - K*Am),
      d*(W2*v*(2*K - 1)*Am*(u - lm)),
      d*(W2*w*(2*K - 1)*Am*(u - lm))},
     {e*(hW*Ap*(u-lp) - u - W2*(V2 - u*u)*(2*K - 1)*(u - Ap*lp) + K*Ap*lp),
      e*(-u - W2*(V2 - u*u)*(2*K - 1)*(u - Ap*lp) + K*Ap*lp),
      e*(1 + W2*(V2 - u*u)*(2*K - 1)*(1 - Ap) - K*Ap),
      e*(W2*v*(2*K - 1)*Ap*(u - lp)),
      e*(W2*w*(2*K - 1)*Ap*(u - lp))},
     {a*(h-W), -a*W, a*W*u, a*W*v, a*W*w},
     {-b*v, -b*v, b*u*v, b*(1-u*u), 0},
     {-c*w, -c*w, c*u*w, 0, c*(1-u*u)}};

  matrix_matrix_product(S[dim-1][0], RR[0], R, 5, 5, 5);
  matrix_matrix_product(LL[0], T[dim-1][0], L, 5, 5, 5);

  lam[0] = lm;
  lam[1] = lp;
  lam[2] = u;
  lam[3] = u;
  lam[4] = u;
}
