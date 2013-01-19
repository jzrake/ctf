
/*------------------------------------------------------------------------------
 * FILE: riemann_hlld-rmhd.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * PURPOSE: Obtain the intercell flux between two constant states using a 5-wave
 *   approximation.
 *
 * REFERENCES: Mignone, Ugliano, & Bodo (2008) A ﬁve-wave Harten–Lax–van Leer
 *   Riemann solver for relativistic magnetohydrodynamics
 *
 *------------------------------------------------------------------------------
 */

#include <cstdio>
#include "nrsolver.hpp"
#include "secant.hpp"
#include "riemann_hlld-rmhd.hpp"
#include "riemann_hllc.hpp"
#include "config.h"

static double PrimFromHllState[8];
static double ConsFromHllState[8];
static double FluxFromHllState[8];

int HlldRmhdRiemannSolver::NumberOfSuccess = 0;
int HlldRmhdRiemannSolver::NumberOfFailure = 0;


class HlldEquation48 : public EquationSystemBaseClass
{
private:
  enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
  enum { rho, pre, vx, vy, vz };             // Primitive

  struct HlldOuterwaveJumpState
  {
    double E,D;
    double B[3], m[3];
  } ;
  struct HlldInnerwaveState
  {
    double U[8];
    double B[3], v[3];
  } ;

  const double *Pl, *Pr;
  double *Ul, *Ur;
  double *Fl, *Fr;
  const double ap, am;
  HlldOuterwaveJumpState L, R;

public:
  HlldEquation48(const double *Pl, const double *Pr,
                 double *Ul, double *Ur,
                 double *Fl, double *Fr,
                 double ap, double am);
  void Function(const double *x, double *y) const;
  void Jacobian(const double *x, double *J) const;
  void SampleSolution(double p21, double x2t, double *P) const;
  void ReconstructSolution(double p, double *U, double *F, double s) const;
  void PrintFunction(double p0, double p1) const;
  void PrintVariables(double p) const;
  double EstimateSolution(int attempt);

private:
  void FillInnerWaveStates(double p, double *Kl, double *Kr, double *Bc,
                           double &etaL, double &etaR, double &wL, double &wR,
                           HlldInnerwaveState &aL, HlldInnerwaveState &aR) const;

  inline double dot(const double v[3], const double w[3]) const
  {
    return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
  }
  inline double dot(const double v[3]) const
  {
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
  }

public:
  class Failure : public std::exception { };
  class NegativeIntermediatePressure : public Failure
  {
  public:
    virtual const char *what() const throw()
    {
      return "HLLD solver converged to a negative intermediate pressure.";
    }
  } ;
  class SuperLuminalWavespeed : public Failure
  {
  public:
    virtual const char *what() const throw()
    {
      return "HLLD solver found superluminal wave speeds.";
    }
  } ;
  class Equation54Violated : public Failure
  {
  public:
    virtual const char *what() const throw()
    {
      return "HLLD solver found unacceptable solution, see equation (54).";
    }
  } ;
  class JumpInNormalBfield : public Failure
  {
  public:
    virtual const char *what() const throw()
    {
      return "HLLD solver found normal B-field was discontinuous.";
    }
  } ;
  class InfiniteIntermediatePressure : public Failure
  {
  public:
    virtual const char *what() const throw()
    {
      return "Newton-Rapheson solver hit infinite pressure.";
    }
  } ;
  class UnmatchedContactSpeed : public Failure
  {
  public:
    virtual const char *what() const throw()
    {
      return "HLLD found normal jump in velocities across the contact wave.";
    }
  } ;
  class DegenerateWaves : public Failure
  {
  public:
    virtual const char *what() const throw()
    {
      return "HLLD found some waves are degenerate, solution may be unreliable";
    }
  } ;
} ;



int HlldRmhdRiemannSolver::IntercellFlux(const double *pl, const double *pr, double *u,
                                         double *f, double s, int dim)
{
  AdiabaticIdealRmhd &fluid = Mara->GetFluid<AdiabaticIdealRmhd>();

  // Convention here is that user input variables u,f,pl/pr are lower case, and
  // get indexed with B1,B2,B3 etc. Local variables U,F,Pl/Pr are input for HLLD
  // and get indexed with Bx,By,Bz etc. HLLD expects all intercell fluxes to be
  // normal to the x-axis, so we permute the input and output states.
  // ---------------------------------------------------------------------------
  int i;
  double epl, epr, eml, emr;
  double Ul[8], Ur[8];
  double Pl[8], Pr[8];
  double Fl[8], Fr[8];

  int B1=0,B2=0,B3=0;
  int S1=0,S2=0,S3=0;
  int v1=0,v2=0,v3=0;

  switch (dim) {
  case 1:
    B1=Bx; B2=By; B3=Bz;
    S1=Sx; S2=Sy; S3=Sz;
    v1=vx; v2=vy; v3=vz;
    break;

  case 2:
    B1=By; B2=Bz; B3=Bx;
    S1=Sy; S2=Sz; S3=Sx;
    v1=vy; v2=vz; v3=vx;
    break;

  case 3:
    B1=Bz; B2=Bx; B3=By;
    S1=Sz; S2=Sx; S3=Sy;
    v1=vz; v2=vx; v3=vy;
    break;
  }

  Pl[pre] = pl[pre];  Pr[pre] = pr[pre];
  Pl[rho] = pl[rho];  Pr[rho] = pr[rho];
  Pl[vx ] = pl[v1];   Pr[vx ] = pr[v1];
  Pl[vy ] = pl[v2];   Pr[vy ] = pr[v2];
  Pl[vz ] = pl[v3];   Pr[vz ] = pr[v3];
  Pl[Bx ] = pl[B1];   Pr[Bx ] = pr[B1];
  Pl[By ] = pl[B2];   Pr[By ] = pr[B2];
  Pl[Bz ] = pl[B3];   Pr[Bz ] = pr[B3];

  Pl[Bx] = Pr[Bx] = 0.5*(pl[B1] + pr[B1]); // This is to prevent jumps in Bx

  if (fluid.PrimToCons(Pl, Ul)) return 1;
  if (fluid.PrimToCons(Pr, Ur)) return 1;

  fluid.FluxAndEigenvalues(Ul, Pl, Fl, &epl, &eml, 1);
  fluid.FluxAndEigenvalues(Ur, Pr, Fr, &epr, &emr, 1);

  const double am = (eml<emr) ? eml : emr;
  const double ap = (epl>epr) ? epl : epr;

  double ml = (fabs(am)<fabs(ap)) ? fabs(ap) : fabs(am);
  if (MaxLambda < ml) MaxLambda = ml;

  double *P_hll = PrimFromHllState;
  double *U_hll = ConsFromHllState;
  double *F_hll = FluxFromHllState;

  for (i=0; i<8; ++i) {
    P_hll[i] = 0.5*(Pl[i] + Pr[i]);
    U_hll[i] = (ap*Ur[i] - am*Ul[i] +       (Fl[i] - Fr[i])) / (ap - am);
    F_hll[i] = (ap*Fl[i] - am*Fr[i] + ap*am*(Ur[i] - Ul[i])) / (ap - am);
  }
  const int hll_c2p_failed = fluid.ConsToPrim(U_hll, P_hll);

  double U[8], F[8];
  HlldEquation48 eqn(Pl, Pr, Ul, Ur, Fl, Fr, ap, am);

  NewtonRaphesonSolver solver1(15, 1e-12);
  SecantMethodSolver   solver2(15, 1e-12);

  int Attempt=0, HlldSuccess=0;
  double p_star;

  while (!hll_c2p_failed && Attempt < 3 && !HlldSuccess) {
    try {
      if (fabs(Pl[Bx]) < 1e-8 || fabs(Pr[Bx]) < 1e-8) {
        p_star = eqn.EstimateSolution(1);
        eqn.ReconstructSolution(p_star, U, F, s);
      }
      else {
        p_star = eqn.EstimateSolution(Attempt);
        solver1.Solve(eqn, &p_star);
        eqn.ReconstructSolution(p_star, U, F, s);
      }
      HlldSuccess = true;
    }
    catch (const HlldEquation48::Failure &e) { }
    catch (const NewtonRaphesonSolver::Failure &e) { }

    ++Attempt;
  }

  if (HlldSuccess) {
    ++NumberOfSuccess;

    // Flux / Conserved states rotated back
    // ---------------------------------------------------------------------------
    if (u) {
      u[ddd] = U[ddd];
      u[tau] = U[tau];
      u[S1 ] = U[Sx];
      u[S2 ] = U[Sy];
      u[S3 ] = U[Sz];
      u[B1 ] = U[Bx];
      u[B2 ] = U[By];
      u[B3 ] = U[Bz];
    }
    if (f) {
      f[ddd] = F[ddd];
      f[tau] = F[tau];
      f[S1 ] = F[Sx];
      f[S2 ] = F[Sy];
      f[S3 ] = F[Sz];
      f[B1 ] = F[Bx];
      f[B2 ] = F[By];
      f[B3 ] = F[Bz];
    }
    return 0;
  }
  else {
    ++NumberOfFailure;
    HllcRmhdRiemannSolver hllc;
    return hllc.IntercellFlux(pl, pr, u, f, s, dim);
  }
}

HlldEquation48::HlldEquation48(const double *Pl, const double *Pr,
                               double *Ul, double *Ur,
                               double *Fl, double *Fr,
                               double ap, double am)
  : EquationSystemBaseClass(1),
    Pl(Pl), Pr(Pr),
    Ul(Ul), Ur(Ur),
    Fl(Fl), Fr(Fr),
    ap(ap), am(am)
{
  // Change the convention of total energy to be consistent with Mignone
  // ---------------------------------------------------------------------------
  Ul[tau] += Ul[ddd]; // Change in convention of total energy
  Fl[tau] += Fl[ddd];
  Ur[tau] += Ur[ddd];
  Fr[tau] += Fr[ddd];

  L.E = am*Ul[tau] - Fl[tau];
  R.E = ap*Ur[tau] - Fr[tau];
  L.D = am*Ul[ddd] - Fl[ddd];
  R.D = ap*Ur[ddd] - Fr[ddd];

  // Fill in the jump across the fast waves, eqn (12)
  // ---------------------------------------------------------------------------
  for (int k=0; k<3; ++k) {
    L.B[k] = am*Ul[Bx+k] - Fl[Bx+k];
    R.B[k] = ap*Ur[Bx+k] - Fr[Bx+k];
    L.m[k] = am*Ul[Sx+k] - Fl[Sx+k];
    R.m[k] = ap*Ur[Sx+k] - Fr[Sx+k];
  }
}

void HlldEquation48::Jacobian(const double *x, double *J) const
{
  const double dp = 1e-8;
  const double p0 = x[0]-0.5*dp;
  const double p1 = x[0]+0.5*dp;
  double f0, f1;
  Function(&p0, &f0);
  Function(&p1, &f1);
  J[0] = (f1-f0)/dp;
}
void HlldEquation48::FillInnerWaveStates(double p, double *Kl, double *Kr, double *Bc,
                                         double &etaL, double &etaR, double &wL, double &wR,
                                         HlldInnerwaveState &aL, HlldInnerwaveState &aR) const
/* --------------------------------------------------------------------------------------
   Notes:
   - Kl(r) is the same in the aL(R) and cL(R) region.
   - Bc holds everywhere between the Alfven modes, it does not change across the contact.
   - In denominator of (45) Mignone is referring to the fast wavespeeds and L/R Bx value.
   -------------------------------------------------------------------------------------- */
{
  {
    double A = L.m[0] - am*L.E + p*(1 - am*am);
    double G = L.B[1]*L.B[1] + L.B[2]*L.B[2];
    double C = L.B[1]*L.m[1] + L.B[2]*L.m[2];
    double Q = -A - G + (Ul[Bx]*Ul[Bx])*(1 - am*am);
    double X = Ul[Bx]*(A*am*Ul[Bx] + C) - (A + G)*(am*p + L.E);

    if (fabs(X) < 1e-8) X = 1e-8;

    aL.v[0] = (Ul[Bx]*(A*Ul[Bx] + am*C) - (A + G)*(p + L.m[0]))  / X; // eqn (23)
    aL.v[1] = (Q*L.m[1] + L.B[1]*(C + Ul[Bx]*(am*L.m[0] - L.E))) / X; // eqn (24)
    aL.v[2] = (Q*L.m[2] + L.B[2]*(C + Ul[Bx]*(am*L.m[0] - L.E))) / X; // eqn (25)
    aL.B[0] = (L.B[0] - Ul[Bx]*aL.v[0]) / (am - aL.v[0]);             // eqn (21)
    aL.B[1] = (L.B[1] - Ul[Bx]*aL.v[1]) / (am - aL.v[0]);
    aL.B[2] = (L.B[2] - Ul[Bx]*aL.v[2]) / (am - aL.v[0]);
  }
  {
    double A = R.m[0] - ap*R.E + p*(1 - ap*ap);
    double G = R.B[1]*R.B[1] + R.B[2]*R.B[2];
    double C = R.B[1]*R.m[1] + R.B[2]*R.m[2];
    double Q = -A - G + (Ur[Bx]*Ur[Bx])*(1 - ap*ap);
    double X = Ur[Bx]*(A*ap*Ur[Bx] + C) - (A + G)*(ap*p + R.E);

    if (fabs(X) < 1e-8) X = 1e-8;

    aR.v[0] = (Ur[Bx]*(A*Ur[Bx] + ap*C) - (A + G)*(p + R.m[0]))  / X; // eqn (23)
    aR.v[1] = (Q*R.m[1] + R.B[1]*(C + Ur[Bx]*(ap*R.m[0] - R.E))) / X; // eqn (24)
    aR.v[2] = (Q*R.m[2] + R.B[2]*(C + Ur[Bx]*(ap*R.m[0] - R.E))) / X; // eqn (25)
    aR.B[0] = (R.B[0] - Ur[Bx]*aR.v[0]) / (ap - aR.v[0]);             // eqn (21)
    aR.B[1] = (R.B[1] - Ur[Bx]*aR.v[1]) / (ap - aR.v[0]);
    aR.B[2] = (R.B[2] - Ur[Bx]*aR.v[2]) / (ap - aR.v[0]);
  }
  wL = p + (L.E - dot(aL.v,L.m)) / (am - aL.v[0]); // eqn (31)
  wR = p + (R.E - dot(aR.v,R.m)) / (ap - aR.v[0]);

  if (wL<0) wL *= -1;
  if (wR<0) wR *= -1;

  etaL = -((Ul[Bx]>0)-(Ul[Bx]<0)) * sqrt(wL);
  etaR = +((Ur[Bx]>0)-(Ur[Bx]<0)) * sqrt(wR);

  for (int k=0; k<3; ++k) {
    Kl[k] = (L.m[k] + p*(k==0) + L.B[k]*etaL) / (am*p + L.E + Ul[Bx]*etaL); // eqn (43)
    Kr[k] = (R.m[k] + p*(k==0) + R.B[k]*etaR) / (ap*p + R.E + Ur[Bx]*etaR);
    Bc[k] = ((aR.B[k]*(Kr[0] - aR.v[0]) + aR.B[0]*aR.v[k]) -                // eqn (45)
             (aL.B[k]*(Kl[0] - aL.v[0]) + aL.B[0]*aL.v[k])) / (Kr[0] - Kl[0]);

    // The following trick assures that the Bc's are not nan's. In the case of
    // the denominator of (45) being zero, the contact region also has zero
    // width and so it won't matter what Bc is set to.
    // -------------------------------------------------------------------------
    if (fabs(Kr[0] - Kl[0]) < 1e-10) Bc[k] = 0.0;
  }
}

void HlldEquation48::Function(const double *x, double *y) const
{
  const double &p = x[0];

  if (fabs(Pl[Bx] - Pr[Bx]) > 1e-12) {
    throw JumpInNormalBfield();
  }
  if (Mara_isinf_cxx(p)) {
    throw InfiniteIntermediatePressure();
  }

  double Kl[3], Kr[3], Bc[3];
  double etaL, etaR, wL, wR;
  HlldInnerwaveState aL, aR;
  FillInnerWaveStates(p, Kl, Kr, Bc, etaL, etaR, wL, wR, aL, aR);

  const double delK = Kr[0] - Kl[0];
  const double Bch[3] = { Bc[0]*delK, Bc[1]*delK, Bc[2]*delK };
  const double Yl = (1 - dot(Kl)) / (etaL*delK - dot(Kl, Bch));
  const double Yr = (1 - dot(Kr)) / (etaR*delK - dot(Kr, Bch));

  if (fabs(delK) < 1e-10) {
    *y = 0.0;
  }
  else if (fabs(Bc[0]) < 1e-10) {
    *y = delK;
  }
  else {
    *y = delK*(1 - Bc[0]*(Yr - Yl));
  }
}
void HlldEquation48::ReconstructSolution
(double p, double *U, double *F, double s) const
{
  double Kl[3], Kr[3], Bc[3];
  double etaL, etaR, wL, wR;
  HlldInnerwaveState aL, aR, cL, cR;
  FillInnerWaveStates(p, Kl, Kr, Bc, etaL, etaR, wL, wR, aL, aR);

  const double denomL = etaL - dot(Kl,Bc);
  const double denomR = etaR - dot(Kr,Bc);
  for (int k=0; k<3; ++k) {
    cL.v[k] = Kl[k] - Bc[k]*(1 - dot(Kl)) / denomL; // eqn (47)
    cR.v[k] = Kr[k] - Bc[k]*(1 - dot(Kr)) / denomR;
    if (fabs(denomL) < 1e-10) cL.v[k] = aL.v[k];
    if (fabs(denomR) < 1e-10) cR.v[k] = aR.v[k];
  }

  const double lamC = cR.v[0]; // speed of the contact wave (L/R equivalent)
  const double lmaL =   Kl[0]; // speed of the L/R Alfven waves
  const double lmaR =   Kr[0];
  const double TOL  =   1e-14; // positive tolerance ~ weaker rejection criterion

  if (p < 0.0) {
    throw NegativeIntermediatePressure();
  }
  if (lmaL < -1.0 || lmaR > 1.0 || fabs(lamC) > 1.0) {
    throw SuperLuminalWavespeed();
  }
  if (!((wL > p) && (aL.v[0] > am-TOL) && (cL.v[0] > lmaL-TOL))) {
    throw Equation54Violated();
  }
  if (!((wR > p) && (aR.v[0] < ap+TOL) && (cR.v[0] < lmaR+TOL))) {
    throw Equation54Violated();
  }
  if (fabs(cR.v[0] - cL.v[0]) > 1e-8) {
    throw UnmatchedContactSpeed();
  }

  double denom;
  int i;
  // States between the fast and Alfven waves
  // ---------------------------------------------------------------------------
  if (fabs(denom = am - aL.v[0]) < 1e-8) {
    for (i=0; i<8; ++i) aL.U[i] = Ul[i];
  }
  else {
    const double vB = dot(aL.v,aL.B);
    aL.U[ddd] =   L.D                           / denom;
    aL.U[tau] = ( L.E + p*aL.v[0] - vB*aL.B[0]) / denom;
    aL.U[Sx ] = (aL.U[tau] + p)*aL.v[0] - vB*aL.B[0];
    aL.U[Sy ] = (aL.U[tau] + p)*aL.v[1] - vB*aL.B[1];
    aL.U[Sz ] = (aL.U[tau] + p)*aL.v[2] - vB*aL.B[2];
    aL.U[Bx ] =  aL.B[0];
    aL.U[By ] =  aL.B[1];
    aL.U[Bz ] =  aL.B[2];
  }

  if (fabs(denom = ap - aR.v[0]) < 1e-8) {
    for (i=0; i<8; ++i) aR.U[i] = Ur[i];
  }
  else {
    const double vB = dot(aR.v,aR.B);
    aR.U[ddd] =   R.D                           / denom;
    aR.U[tau] = ( R.E + p*aR.v[0] - vB*aR.B[0]) / denom;
    aR.U[Sx ] = (aR.U[tau] + p)*aR.v[0] - vB*aR.B[0];
    aR.U[Sy ] = (aR.U[tau] + p)*aR.v[1] - vB*aR.B[1];
    aR.U[Sz ] = (aR.U[tau] + p)*aR.v[2] - vB*aR.B[2];
    aR.U[Bx ] =  aR.B[0];
    aR.U[By ] =  aR.B[1];
    aR.U[Bz ] =  aR.B[2];
  }

  // States between the Alfven and contact waves
  // ---------------------------------------------------------------------------
  if (fabs(denom = lmaL - lamC) < 1e-8) {
    for (i=0; i<8; ++i) cL.U[i] = aL.U[i];
  }
  else {
    const double vBc = dot(cL.v,Bc);
    cL.U[ddd] =  aL.U[ddd] * (lmaL - aL.v[0])                    / denom;
    cL.U[tau] = (lmaL*aL.U[tau] - aL.U[Sx] + p*lamC - vBc*Bc[0]) / denom;
    cL.U[Sx ] = (cL.U[tau] + p)*cL.v[0] - vBc*Bc[0];
    cL.U[Sy ] = (cL.U[tau] + p)*cL.v[1] - vBc*Bc[1];
    cL.U[Sz ] = (cL.U[tau] + p)*cL.v[2] - vBc*Bc[2];
    cL.U[Bx ] =  Bc[0];
    cL.U[By ] =  Bc[1];
    cL.U[Bz ] =  Bc[2];
  }

  if (fabs(denom = lmaR - lamC) < 1e-8) {
    for (i=0; i<8; ++i) cR.U[i] = aR.U[i];
  }
  else {
    const double vBc = dot(cR.v,Bc);
    cR.U[ddd] = aR.U[ddd] * (lmaR - aR.v[0])                     / denom;
    cR.U[tau] = (lmaR*aR.U[tau] - aR.U[Sx] + p*lamC - vBc*Bc[0]) / denom;
    cR.U[Sx ] = (cR.U[tau] + p)*cR.v[0] - vBc*Bc[0];
    cR.U[Sy ] = (cR.U[tau] + p)*cR.v[1] - vBc*Bc[1];
    cR.U[Sz ] = (cR.U[tau] + p)*cR.v[2] - vBc*Bc[2];
    cR.U[Bx ] =  Bc[0];
    cR.U[By ] =  Bc[1];
    cR.U[Bz ] =  Bc[2];
  }

  // All done, just fill in the flux F and the state U
  // ---------------------------------------------------------------------------
  if      (          s<=am  ) for (i=0; i<8; ++i) F[i] = Fl[i];
  else if (  am<s && s<=lmaL) for (i=0; i<8; ++i) F[i] = Fl[i] + am*(aL.U[i]-Ul[i]);
  else if (lmaL<s && s<=lamC) for (i=0; i<8; ++i) F[i] = Fl[i] + am*(aL.U[i]-Ul[i]) + lmaL*(cL.U[i]-aL.U[i]);
  else if (lamC<s && s<=lmaR) for (i=0; i<8; ++i) F[i] = Fr[i] + ap*(aR.U[i]-Ur[i]) + lmaR*(cR.U[i]-aR.U[i]);
  else if (lmaR<s && s<=ap  ) for (i=0; i<8; ++i) F[i] = Fr[i] + ap*(aR.U[i]-Ur[i]);
  else if (  ap<s           ) for (i=0; i<8; ++i) F[i] = Fr[i];

  if      (          s<=am  ) for (i=0; i<8; ++i) U[i] =   Ul[i];
  else if (  am<s && s<=lmaL) for (i=0; i<8; ++i) U[i] = aL.U[i];
  else if (lmaL<s && s<=lamC) for (i=0; i<8; ++i) U[i] = cL.U[i];
  else if (lamC<s && s<=lmaR) for (i=0; i<8; ++i) U[i] = cR.U[i];
  else if (lmaR<s && s<=ap  ) for (i=0; i<8; ++i) U[i] = aR.U[i];
  else if (  ap<s           ) for (i=0; i<8; ++i) U[i] =   Ur[i];

  F[tau] -= F[ddd];
  U[tau] -= U[ddd];
}
double HlldEquation48::EstimateSolution(int attempt)
{
  // Try the HLL total pressure of the input states first
  // ---------------------------------------------------------------------------
  if (attempt == 0) {
    const double *P   =   PrimFromHllState;
    const double V2   =   P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz];
    const double B2   =   P[Bx]*P[Bx] + P[By]*P[By] + P[Bz]*P[Bz];
    const double Bv   =   P[Bx]*P[vx] + P[By]*P[vy] + P[Bz]*P[vz];
    const double W2   =   1.0 / (1.0 - V2);
    const double W    =   sqrt(W2);
    const double b0   =   W*Bv;
    const double b2   =   (B2 + b0*b0) / W2;

    return P[pre] + 0.5*b2;
  }
  // Next try equation (55)
  // ---------------------------------------------------------------------------
  else if (attempt == 1) {
    double *F_hll = FluxFromHllState;
    double *U_hll = ConsFromHllState;

    F_hll[tau] += F_hll[ddd];
    U_hll[tau] += U_hll[ddd];

    const double a = 1.0;
    const double b = U_hll[tau] - F_hll[Sx];
    const double c = F_hll[tau] * U_hll[Sx] - U_hll[tau] * F_hll[Sx];
    const double p1 = fabs((-b + sqrt(b*b - 4*a*c))/(2*a));
    const double p2 = fabs((-b - sqrt(b*b - 4*a*c))/(2*a));

    F_hll[tau] -= F_hll[ddd];
    U_hll[tau] -= U_hll[ddd];

    return (p1<p2) ? p1 : p2;
  }
  // Try p=0 until they give up
  // ---------------------------------------------------------------------------
  else return 0.0;
}

void HlldEquation48::PrintFunction(double p0, double p1) const
{
  FILE *outf = fopen("eqn48.dat", "w");
  const int NumP = 10000;
  const double dp = (p1-p0) / NumP;

  for (int i=0; i<NumP; ++i) {
    double p = p0 + i*dp;
    double eqn48;
    Function(&p, &eqn48);
    fprintf(outf, "%8.6f %8.6f\n", p, eqn48);
  }

  fclose(outf);
}
