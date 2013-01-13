
/*------------------------------------------------------------------------------
 * FILE: riemann_hll.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * PURPOSE: Obtain the intercell flux between two constant states using a 2-wave
 *   approximation.
 *
 * REFERENCES:
 *
 *------------------------------------------------------------------------------
 */

#include <cstring>
#include "riemann_hll.hpp"
#define MAXNQ 8 // Used for static array initialization


int HllRiemannSolver::IntercellFlux
(const double *pl, const double *pr, double *U, double *F, double s, int dim)
{
  FluidEquations &fluid = *Mara->fluid;

  int i, NQ = fluid.GetNq();
  double epl=1, epr=1, eml=1, emr=1;
  double Ul[MAXNQ], Ur[MAXNQ];
  double Pl[MAXNQ], Pr[MAXNQ];
  double Fl[MAXNQ], Fr[MAXNQ];

  std::memcpy(Pl,pl,NQ*sizeof(double));
  std::memcpy(Pr,pr,NQ*sizeof(double));

  if (fluid.PrimToCons(Pl,Ul)) return 1;
  if (fluid.PrimToCons(Pr,Ur)) return 1;

  fluid.FluxAndEigenvalues(Ul, Pl, Fl, &epl, &eml, dim);
  fluid.FluxAndEigenvalues(Ur, Pr, Fr, &epr, &emr, dim);

  double ap = (epl>epr) ? epl : epr;
  double am = (eml<emr) ? eml : emr;

  double ml = (fabs(am)<fabs(ap)) ? fabs(ap) : fabs(am);
  if (MaxLambda < ml) MaxLambda = ml;

  double F_hll[MAXNQ], U_hll[MAXNQ];
  for (i=0; i<NQ; ++i) {
    U_hll[i] = (ap*Ur[i] - am*Ul[i] +       (Fl[i] - Fr[i])) / (ap - am);
    F_hll[i] = (ap*Fl[i] - am*Fr[i] + ap*am*(Ur[i] - Ul[i])) / (ap - am);
  }

  if (U) {
    if      (         s<=am ) for (i=0; i<NQ; ++i) U[i] = Ul   [i];
    else if ( am<s && s<=ap ) for (i=0; i<NQ; ++i) U[i] = U_hll[i];
    else if ( ap<s          ) for (i=0; i<NQ; ++i) U[i] = Ur   [i];
  }
  {
    if      (         s<=am ) for (i=0; i<NQ; ++i) F[i] = Fl   [i];
    else if ( am<s && s<=ap ) for (i=0; i<NQ; ++i) F[i] = F_hll[i];
    else if ( ap<s          ) for (i=0; i<NQ; ++i) F[i] = Fr   [i];
  }

  return 0;
}
