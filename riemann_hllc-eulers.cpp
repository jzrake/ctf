
/*------------------------------------------------------------------------------
 * FILE: riemann_hllc-eulers.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * PURPOSE: Obtain the intercell flux between two constant states using a 3-wave
 *   approximation.
 *
 * REFERENCES: Toro (1999) Riemann Solvers and Numerical Methods for Fluid Dynamics
 *
 *------------------------------------------------------------------------------
 */

#include <cmath>
#include <cstring>
#include "riemann_hllc.hpp"

enum { rho, nrg, px, py, pz }; // Conserved
enum { RHO, pre, vx, vy, vz }; // Primitive

int HllcEulersRiemannSolver::IntercellFlux
(const double *pl, const double *pr, double *U, double *F, double s, int dim)
{
  AdiabaticIdealEulers &fluid = Mara->GetFluid<AdiabaticIdealEulers>();

  int i;
  double epl, epr, eml, emr;
  double Ul[5], Ur[5];
  double Pl[5], Pr[5];
  double Fl[5], Fr[5];

  int p1=0,p2=0,p3=0;
  int v1=0,v2=0,v3=0;

  switch (dim) {
  case 1:
    p1=px; p2=py; p3=pz;
    v1=vx; v2=vy; v3=vz;
    break;

  case 2:
    p1=py; p2=pz; p3=px;
    v1=vy; v2=vz; v3=vx;
    break;

  case 3:
    p1=pz; p2=px; p3=py;
    v1=vz; v2=vx; v3=vy;
    break;
  }

  std::memcpy(Pl,pl,5*sizeof(double));
  std::memcpy(Pr,pr,5*sizeof(double));

  fluid.PrimToCons(Pl,Ul);
  fluid.PrimToCons(Pr,Ur);

  fluid.FluxAndEigenvalues(Ul, Pl, Fl, &epl, &eml, dim);
  fluid.FluxAndEigenvalues(Ur, Pr, Fr, &epr, &emr, dim);

  double ap = (epl>epr) ? epl : epr;
  double am = (eml<emr) ? eml : emr;

  double ml = (fabs(am)<fabs(ap)) ? fabs(ap) : fabs(am);
  if (MaxLambda < ml) MaxLambda = ml;

  double Ul_[5], Ur_[5]; // The star states
  double lc = ((Pr[pre] - Pr[rho]*Pr[v1]*(ap - Pr[v1])) -
	       (Pl[pre] - Pl[rho]*Pl[v1]*(am - Pl[v1]))) /
    (Pl[rho]*(am - Pl[v1]) - Pr[rho]*(ap - Pr[v1])); // eqn 10.58
  {
    double fact = Pl[rho] * (am - Pl[v1]) / (am - lc); // eqn 10.33
    Ul_[rho] = fact;
    Ul_[nrg] = fact * (Ul[nrg]/Pl[rho] + (lc - Pl[v1])*(lc + Pl[pre]/(Pl[rho]*(am-Pl[v1]))));
    Ul_[p1 ] = fact * lc;
    Ul_[p2 ] = fact * Pl[v2];
    Ul_[p3 ] = fact * Pl[v3];
  }
  {
    double fact = Pr[rho] * (ap - Pr[v1]) / (ap - lc); // eqn 10.33
    Ur_[rho] = fact;
    Ur_[nrg] = fact * (Ur[nrg]/Pr[rho] + (lc - Pr[v1])*(lc + Pr[pre]/(Pr[rho]*(ap-Pr[v1]))));
    Ur_[p1 ] = fact * lc;
    Ur_[p2 ] = fact * Pr[v2];
    Ur_[p3 ] = fact * Pr[v3];
  }

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

  return 0;
}
