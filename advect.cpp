
/*------------------------------------------------------------------------------
 * FILE: euler.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * REFERENCES:
 *
 *------------------------------------------------------------------------------
 */

#include <cstdio>
#include <cmath>
#include "advect.hpp"

typedef ScalarAdvection Advect;

Advect::ScalarAdvection(double WaveSpeed) : WaveSpeed(WaveSpeed) { }

std::vector<std::string> Advect::GetPrimNames() const
{
  std::string vars[1] = { "u" };
  return std::vector<std::string>(vars, vars+1);
}
std::string Advect::PrintPrim(const double *P) const
{
  char str[512];
  sprintf(str, "u = { %+6.4e }", P[0]);
  return std::string(str);
}
std::string Advect::PrintCons(const double *U) const
{
  char str[512];
  sprintf(str, "u = { %+6.4e }", U[0]);
  return std::string(str);
}
int Advect::ConsToPrim(const double *U, double *P) const
{
  P[0] = U[0];
  return 0;
}
int Advect::PrimToCons(const double *P, double *U) const
{
  U[0] = P[0];
  return 0;
}
void Advect::FluxAndEigenvalues(const double *U,
                                const double *P, double *F,
                                double *ap, double *am, int dimension) const
{
  *F = WaveSpeed*U[0];

  if (ap) *ap = WaveSpeed > 0.0 ? WaveSpeed : 0.0;
  if (am) *am = WaveSpeed < 0.0 ? WaveSpeed : 0.0;
}

void Advect::Eigensystem(const double *U, const double *P,
                         double *L, double *R, double *lam, int dim) const
{
  L[0] = R[0] = 1.0;
  lam[0] = WaveSpeed;
}
