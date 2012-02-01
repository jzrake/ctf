
/*------------------------------------------------------------------------------
 * FILE: riemann_exact-eulers.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP
 *
 * DESCRIPTION:
 *
 * REFERENCES: E.F Toro (1999)
 *             Riemann Solvers and Numerical Methods for Fluid Dynamics
 *------------------------------------------------------------------------------
 */

#include <cmath>
#include <cstring>
#include "logging.hpp"
#include "eos.hpp"
#include "riemann_exact-eulers.hpp"



class EulersWavePattern : public EquationSystemBaseClass
{
private:
  const double AdiabaticGamma;
  const double *Pl, *Pr;

  int v1, v2, v3;
  double gm0, gm1, gm2, gm3, gm4, gm5;
  double pL, pR, uL, uR, aL, aR, AL, AR, BL, BR;

  enum { rho, nrg, px, py, pz }; // Conserved
  enum { RHO, pre, vx, vy, vz }; // Primitive

public:
  EulersWavePattern(const double *Pl, const double *Pr, double AdiabaticGamma, int dim);
  void Function(const double *x, double *y) const;
  void Jacobian(const double *x, double *J) const;
  void SampleSolution(double p21, double x2t, double *P) const;
  double EstimateSolution(int attempt);
} ;



int ExactEulersRiemannSolver::IntercellFlux(const double *pl, const double *pr,
					    double *U_out, double *F, double s, int dim)
{
  switch (dim) {
  case 1:
    v1=vx; v2=vy; v3=vz;
    break;
  case 2:
    v1=vy; v2=vz; v3=vx;
    break;
  case 3:
    v1=vz; v2=vx; v3=vy;
    break;
  }

  NewtonRaphesonSolver solver(500, 1e-12);
  EulersWavePattern eqn(pl, pr, Mara->GetEos<AdiabaticEos>().Gamma, dim);

  double U[5], P[5], ap, am, p;
  int Attempt = 0;
  while (Attempt < 3) {
    try {
      p = eqn.EstimateSolution(Attempt);
      solver.Solve(eqn, &p);
      eqn.SampleSolution(p, s, P);
      fluid.PrimToCons(P, U);
      fluid.FluxAndEigenvalues(U, P, F, &ap, &am, dim);

      double ml = (fabs(am)<fabs(ap)) ? fabs(ap) : fabs(am);
      if (MaxLambda < ml) MaxLambda = ml;
      if (U_out) std::memcpy(U_out, U, 5*sizeof(double));
      return 0;
    }
    catch (const std::exception &e) {
      DebugLog.Warning(__FUNCTION__)
	<< "failed to find p* in Riemann solution on attempt "
	<< Attempt << std::endl;
      ++Attempt;
    }
  }
  if (BackupRiemannSolver) {
    DebugLog.Error(__FUNCTION__)
      << "failed to find p* in Riemann solution. "
      << "Reverting to backup riemann solver." << std::endl;
    BackupRiemannSolver->IntercellFlux(pl, pr, U_out, F, s, dim);
  }
  else {
    throw NoPressureSolution();
  }
  return 0;
}

EulersWavePattern::EulersWavePattern
(const double *Pl, const double *Pr, double AdiabaticGamma, int dim)
  : EquationSystemBaseClass(1),
    AdiabaticGamma(AdiabaticGamma),
    Pl(Pl), Pr(Pr), v1(vx), v2(vy), v3(vz)
{
  switch (dim) {
  case 1:
    v1=vx; v2=vy; v3=vz;
    break;
  case 2:
    v1=vy; v2=vz; v3=vx;
    break;
  case 3:
    v1=vz; v2=vx; v3=vy;
    break;
  }

  gm0 = AdiabaticGamma;
  gm1 = (gm0+1) / (2*gm0);
  gm2 = (gm0-1) / (2*gm0);
  gm3 = (gm0-1) / (gm0+1);
  gm4 = 2.0 / (gm0+1);
  gm5 = 2.0 / (gm0-1);

  pL = Pl[pre];
  pR = Pr[pre];
  uL = Pl[v1];
  uR = Pr[v1];

  AL = 2.0/((gm0+1)*Pl[rho]); // eqn 4.8
  AR = 2.0/((gm0+1)*Pr[rho]);
  BL = gm3*pL;
  BR = gm3*pR;

  aL = sqrt(gm0*Pl[pre] / Pl[rho]); // sound speed
  aR = sqrt(gm0*Pr[pre] / Pr[rho]);
}

void EulersWavePattern::Function(const double *x, double *y) const
{
  const double p = *x;
  const double fL = (p > pL)       ?     // eqn 4.6
    (p - pL) * sqrt(AL / (p + BL)) :     // shock wave
    2*aL/(gm0-1) * (pow(p/pL, gm2) - 1); // rarefaction wave

  const double fR = (p > pR)       ?     // eqn 4.7
    (p - pR) * sqrt(AR / (p + BR)) :     // shock wave
    2*aR/(gm0-1) * (pow(p/pR, gm2) - 1); // rarefaction wave

  *y = fL + fR + (Pr[v1] - Pl[v1]); // eqn 4.5
}
void EulersWavePattern::Jacobian(const double *x, double *J) const
{
  const double p = *x;
  const double fpL = (p > pL)                       ? // eqn 4.37
    sqrt(AL/(BL + p)) * (1 - 0.5*(p - pL)/(BL + p)) : // shock wave
    1.0/(Pl[rho] * aL) * pow(p/pL, -gm1);             // rarefaction wave

  const double fpR = (p > pR)                       ? // eqn 4.37
    sqrt(AR/(BR + p)) * (1 - 0.5*(p - pR)/(BR + p)) : // shock wave
    1.0/(Pr[rho] * aR) * pow(p/pR, -gm1);             // rarefaction wave

  *J = fpL + fpR;
}
double EulersWavePattern::EstimateSolution(int attempt)
{
  if (attempt == 0) { // eqn 4.46
    const double gm1 = (gm0-1) / (2*gm0);
    const double gm3 = 1.0/gm1;
    return pow(((aL+aR)+0.5*(gm0-1)*(uL-uR)) / (aL/pow(pL,gm1)+aR/pow(pR,gm1)), gm3);
  }
  else if (attempt == 1) { // eqn 4.47
    double pPV = 0.5*(pL+pR) + 0.125*(uL-uR)*(Pl[rho]+Pr[rho])*(aL+aR);
    if (pPV < 1e-6) pPV = 1e-6;
    return pPV;
  }
  else { // eqn 4.49
    return 0.5*(pL + pR);
  }
}
void EulersWavePattern::SampleSolution(double p_, double S, double *P) const
{
  /* Underscore after variable is Toro's (*), which indicates the Star Region */

  const double aL_ = aL*pow(p_/pL, gm2); // eqn 4.54
  const double aR_ = aR*pow(p_/pR, gm2);

  const double fL = (p_ > pL)        ?    // eqn 4.6
    (p_ - pL) * sqrt(AL / (p_ + BL)) :    // shock wave
    2*aL/(gm0-1) * (pow(p_/pL, gm2) - 1); // rarefaction wave

  const double fR = (p_ > pR)        ?    // eqn 4.7
    (p_ - pR) * sqrt(AR / (p_ + BR)) :    // shock wave
    2*aR/(gm0-1) * (pow(p_/pR, gm2) - 1); // rarefaction wave

  const double u_ = 0.5*(uL + uR) - 0.5*(fL - fR); // eqn 4.88

  const double SL  = uL - aL*sqrt(gm1 * p_/pL + gm2); // eqn 4.52
  const double SHL = uL - aL;
  const double STL = u_ - aL_;
  const double STR = u_ + aR_;
  const double SHR = uR + aR;
  const double SR  = uR + aR*sqrt(gm1 * p_/pR + gm2); // eqn 4.59

  const double FrhoL_ = Pl[rho]*pow(p_/pL, 1.0/gm0); // eqn 4.53 (fan)
  const double FrhoR_ = Pr[rho]*pow(p_/pR, 1.0/gm0); // eqn 4.60
  const double SrhoL_ = Pl[rho]*(p_/pL + gm3) / (gm3*p_/pL + 1); // eqn 4.50 (shock)
  const double SrhoR_ = Pr[rho]*(p_/pR + gm3) / (gm3*p_/pR + 1); // eqn 4.57

  /* --------------------------------------------------------------------------
    The logic below is the sampling procedure outlined by Toro in the flow chart
    Figure 4.14. It consists all-together of 10 cases: 6 regions (L/R, L/R fan,
    and L/R star) and then two sub cases for each of L/R fan and L/R star,
    corresponding to whether the L/R going wave is a fan or a shock.
    --------------------------------------------------------------------------- */
  /*                    |   |   |                                               */
  /*                    |   |   |                                               */
  if (S < u_) {         // sampling left of particle velocity characteristic
    if (p_ > pL) {          // left-moving wave is a shock
      if (S < SL) {             // left region
	std::memcpy(P, Pl, 5*sizeof(double));
      }
      else {                    // left star region
	P[rho] = SrhoL_;
	P[pre] = p_;
	P[v1 ] = u_;
	P[v2 ] = Pl[v2];
	P[v3 ] = Pl[v3];
      }
    }
    else {                  // left-moving wave is a rarefaction
      if (S < SHL) {            // left region
	std::memcpy(P, Pl, 5*sizeof(double));
      }
      else {
        if (S < STL) {          // left fan
	  P[rho] = Pl[rho]*pow(gm4 + gm3*(uL - S)/aL, gm5); // eqn 4.56
	  P[pre] = Pl[pre]*pow(gm4 + gm3*(uL - S)/aL, 1.0/gm2);
	  P[v1 ] = gm4*(aL + uL/gm5 + S);
	  P[v2 ] = Pl[v2];
	  P[v3 ] = Pl[v3];
        }
        else {                  // left star region
	  P[rho] = FrhoL_;
	  P[pre] = p_;
	  P[v1 ] = u_;
	  P[v2 ] = Pl[v2];
	  P[v3 ] = Pl[v3];
        }
      }
    }
  }
  /* -------------------------------------------------------------------------- */
  /*                    |   |   |                                               */
  /*                    |   |   |                                               */
  else {                // sampling right of particle velocity characteristic
    if (p_ > pR) {          // right-moving wave is a shock
      if (S > SR) {             // right region
	std::memcpy(P, Pr, 5*sizeof(double));
      }
      else {                    // right star region
	P[rho] = SrhoR_;
	P[pre] = p_;
	P[v1 ] = u_;
	P[v2 ] = Pr[v2];
	P[v3 ] = Pr[v3];
      }
    }
    else {                  // right-moving wave is a rarefaction
      if (S > SHR) {            // right region
	std::memcpy(P, Pr, 5*sizeof(double));
      }
      else {
        if (S > STR) {          // right fan
	  P[rho] = Pr[rho]*pow(gm4 - gm3*(uR - S)/aR, gm5); // eqn 4.63
	  P[pre] = Pr[pre]*pow(gm4 - gm3*(uR - S)/aR, 1.0/gm2);
	  P[v1 ] = gm4*(-aR + uR/gm5 + S);
	  P[v2 ] = Pr[v2];
	  P[v3 ] = Pr[v3];
        }
        else {                  // right star region
	  P[rho] = FrhoR_;
	  P[pre] = p_;
	  P[v1 ] = u_;
	  P[v2 ] = Pr[v2];
	  P[v3 ] = Pr[v3];
        }
      }
    }
  }
}
