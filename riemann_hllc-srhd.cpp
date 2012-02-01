
/*------------------------------------------------------------------------------
 * FILE: riemann_hllc-srhd.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * PURPOSE: Obtain the intercell flux between two constant states using a 3-wave
 *   approximation.
 *
 * REFERENCES: Mignone & Bodo (2005)
 *             An HLLC Riemann solver for relativistic flows - I. Hydrodynamics
 *------------------------------------------------------------------------------
 */

#include <cmath>
#include <cstring>
#include "riemann_hllc.hpp"

enum { ddd, tau, Sx, Sy, Sz }; // Conserved
enum { rho, pre, vx, vy, vz }; // Primitive
static const double SMALL_A  = 1e-10;


int HllcSrhdRiemannSolver::IntercellFlux
(const double *pl, const double *pr, double *U, double *F, double s, int dim)
{
  AdiabaticIdealSrhd &fluid = Mara->GetFluid<AdiabaticIdealSrhd>();

  int i;
  double epl, epr, eml, emr;
  double Ul[5], Ur[5];
  double Pl[5], Pr[5];
  double Fl[5], Fr[5];

  int S1=0,S2=0,S3=0;
  int v1=0,v2=0,v3=0;

  switch (dim) {
  case 1:
    S1=Sx; S2=Sy; S3=Sz;
    v1=vx; v2=vy; v3=vz;
    break;

  case 2:
    S1=Sy; S2=Sz; S3=Sx;
    v1=vy; v2=vz; v3=vx;
    break;

  case 3:
    S1=Sz; S2=Sx; S3=Sy;
    v1=vz; v2=vx; v3=vy;
    break;
  }

  std::memcpy(Pl,pl,5*sizeof(double));
  std::memcpy(Pr,pr,5*sizeof(double));

  if (fluid.PrimToCons(Pl,Ul)) return 1;
  if (fluid.PrimToCons(Pr,Ur)) return 1;

  fluid.FluxAndEigenvalues(Ul, Pl, Fl, &epl, &eml, dim);
  fluid.FluxAndEigenvalues(Ur, Pr, Fr, &epr, &emr, dim);

  Ul[tau] += Ul[ddd];  Fl[tau] += Fl[ddd]; // Change in convention of total energy
  Ur[tau] += Ur[ddd];  Fr[tau] += Fr[ddd];

  double ap = (epl>epr) ? epl : epr;
  double am = (eml<emr) ? eml : emr;

  double ml = (fabs(am)<fabs(ap)) ? fabs(ap) : fabs(am);
  if (MaxLambda < ml) MaxLambda = ml;

  double F_hll[5], U_hll[5];
  for (i=0; i<5; ++i) {
    U_hll[i] = (ap*Ur[i] - am*Ul[i] +       (Fl[i] - Fr[i])) / (ap - am);
    F_hll[i] = (ap*Fl[i] - am*Fr[i] + ap*am*(Ur[i] - Ul[i])) / (ap - am);
  }

  double Ul_[5], Ur_[5], lc; // The star states

  const double a =  F_hll[tau];
  const double b = -F_hll[S1 ] - U_hll[tau];
  const double c =  U_hll[S1 ];

  const double v1_ = lc = (fabs(a) < SMALL_A) ? -c/b : (-b - sqrt(b*b - 4*a*c)) / (2*a);
  const double p_  = -F_hll[tau]*v1_ + F_hll[S1];

  Ul_[ddd] = (am - Pl[v1]) / (am - v1_) * Ul[ddd];
  Ur_[ddd] = (ap - Pr[v1]) / (ap - v1_) * Ur[ddd];

  Ul_[tau] = (am*Ul[tau] - Ul[S1] + p_*v1_) / (am - v1_);
  Ur_[tau] = (ap*Ur[tau] - Ur[S1] + p_*v1_) / (ap - v1_);

  Ul_[S1 ] = (Ul_[tau] + p_)*v1_;
  Ur_[S1 ] = (Ur_[tau] + p_)*v1_;

  Ul_[S2 ] = (am - Pl[v1]) / (am - v1_) * Ul[S2];
  Ur_[S2 ] = (ap - Pr[v1]) / (ap - v1_) * Ur[S2];

  Ul_[S3 ] = (am - Pl[v1]) / (am - v1_) * Ul[S3];
  Ur_[S3 ] = (ap - Pr[v1]) / (ap - v1_) * Ur[S3];

  if (U != 0) {
    if      (         s<=am ) for (i=0; i<5; ++i) U[i] = Ul [i];
    else if ( am<s && s<=lc ) for (i=0; i<5; ++i) U[i] = Ul_[i];
    else if ( lc<s && s<=ap ) for (i=0; i<5; ++i) U[i] = Ur_[i];
    else if ( ap<s          ) for (i=0; i<5; ++i) U[i] = Ur [i];
  }
  {
    if      (         s<=am ) for (i=0; i<5; ++i) F[i] = Fl[i];
    else if ( am<s && s<=lc ) for (i=0; i<5; ++i) F[i] = Fl[i] + am*(Ul_[i]-Ul[i]);
    else if ( lc<s && s<=ap ) for (i=0; i<5; ++i) F[i] = Fr[i] + ap*(Ur_[i]-Ur[i]);
    else if ( ap<s          ) for (i=0; i<5; ++i) F[i] = Fr[i];
  }

  if (U) U[tau] -= U[ddd]; // Change in convention of total energy
  /*  */ F[tau] -= F[ddd];

  return 0;
}
