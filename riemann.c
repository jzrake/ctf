
/*------------------------------------------------------------------------------
 * FILE: riemann.c
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * PURPOSE: Obtains the intercell flux between two constant states
 *
 * POLICY:
 *
 * The input states must contain valid (matching) conserved and primitive
 * values. Fluxes and eigenvalues will then be obtained as needed within. The
 * output state is guarenteed to have a conserved value and a flux. The flux in
 * the case of HLL is not the flux derived from those values, but rather from
 * the HLL formula directly, which satisfies the conservative jump
 * conditions. After the call is executed the user is free to obtain the
 * primitives by calling c2p on the output state, although primitives might
 * already have been calculated, as in the case of the exact solver.
 *
 * ------------------------------------------------------------------------------
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#define FLUIDS_PRIVATE_DEFS
#define FLUIDS_INDEX_VARS
#include "fluids.h"

static const long FLUIDS_FLUX[3] = {FLUIDS_FLUX0, FLUIDS_FLUX1, FLUIDS_FLUX2};
static const long FLUIDS_EVAL[3] = {FLUIDS_EVAL0, FLUIDS_EVAL1, FLUIDS_EVAL2};

static double _soln_f(fluid_riemann *R, double p);
static double _soln_g(fluid_riemann *R, double p);
static double _soln_estimate(fluid_riemann *R, int attempt);
static int _soln_solve(fluid_riemann *R, double *p);
static int _hll_exec(fluid_riemann *R);
static int _hll_sample(fluid_riemann *R, fluid_state *S, double s);
static int _nrhyd_exact_exec(fluid_riemann *R);
static int _nrhyd_exact_sample(fluid_riemann *R, fluid_state *S, double s);

struct fluid_riemann
{
  fluid_state *SL;
  fluid_state *SR;
  double p_solution;
  double gm0, gm1, gm2, gm3, gm4, gm5;
  double pL, pR, uL, uR, aL, aR, AL, AR, BL, BR;
  double ap, am;
  double *U_hll;
  double *F_hll;
  int v1, v2, v3;
  int dim;
  int solver;
} ;

fluid_riemann *fluids_riemann_new()
{
  fluid_riemann *R = (fluid_riemann*) malloc(sizeof(fluid_riemann));
  fluid_riemann riem = {
    .SL = NULL,
    .SR = NULL,
    .gm0 = 0.0,
    .gm1 = 0.0,
    .gm2 = 0.0,
    .gm3 = 0.0,
    .gm4 = 0.0,
    .gm5 = 0.0,
    .pL = 0.0,
    .pR = 0.0,
    .uL = 0.0,
    .uR = 0.0,
    .aL = 0.0,
    .aR = 0.0,
    .AL = 0.0,
    .AR = 0.0,
    .BL = 0.0,
    .BR = 0.0,
    .U_hll = NULL,
    .F_hll = NULL,
    .v1 = 0,
    .v2 = 0,
    .v3 = 0,
    .dim = 0,
    .solver = FLUIDS_RIEMANN_HLL,
  } ;
  *R = riem;
  return R;
}
int fluids_riemann_del(fluid_riemann *R)
{
  free(R->U_hll);
  free(R->F_hll);
  free(R);
  return 0;
}
int fluids_riemann_setdim(fluid_riemann *R, int dim)
{
  R->dim = dim;
  return 0;
}
int fluids_riemann_setstateL(fluid_riemann *R, fluid_state *S)
{
  R->SL = S;
  return 0;
}
int fluids_riemann_setstateR(fluid_riemann *R, fluid_state *S)
{
  R->SR = S;
  return 0;
}
int fluids_riemann_setsolver(fluid_riemann *R, int solver)
{
  R->solver = solver;
  return 0;
}

int fluids_riemann_execute(fluid_riemann *R)
{
  fluids_update(R->SL, FLUIDS_EVAL[R->dim] | FLUIDS_FLUX[R->dim]);
  fluids_update(R->SR, FLUIDS_EVAL[R->dim] | FLUIDS_FLUX[R->dim]);

  switch (R->dim) {
  case 0: R->v1=vx; R->v2=vy; R->v3=vz; break;
  case 1: R->v1=vy; R->v2=vz; R->v3=vx; break;
  case 2: R->v1=vz; R->v2=vx; R->v3=vy; break;
  default: return FLUIDS_ERROR_BADARG;
  }

  switch (R->solver) {
  case FLUIDS_RIEMANN_HLL:
    return _hll_exec(R);
  case FLUIDS_RIEMANN_HLLC:
    return FLUIDS_ERROR_BADREQUEST;
  case FLUIDS_RIEMANN_EXACT:
    return _nrhyd_exact_exec(R);
  default:
    return FLUIDS_ERROR_BADREQUEST;
  }
}

int fluids_riemann_sample(fluid_riemann *R, fluid_state *S, double s)
{
  switch (R->solver) {
  case FLUIDS_RIEMANN_HLL:
    return _hll_sample(R, S, s);
  case FLUIDS_RIEMANN_HLLC:
    return FLUIDS_ERROR_BADREQUEST;
  case FLUIDS_RIEMANN_EXACT:
    _nrhyd_exact_sample(R, S, s);
  default:
    return FLUIDS_ERROR_BADREQUEST;
  }
  return 0;
}

int _hll_exec(fluid_riemann *R)
{
  int nw = R->SL->nwaves;
  double epl = R->SL->eigenvalues[R->dim][nw - 1];
  double epr = R->SR->eigenvalues[R->dim][nw - 1];
  double eml = R->SL->eigenvalues[R->dim][0];
  double emr = R->SR->eigenvalues[R->dim][0];
  double ap = R->ap = (epl>epr) ? epl : epr;
  double am = R->am = (eml<emr) ? eml : emr;
  double *Ul = R->SL->conserved;
  double *Ur = R->SR->conserved;
  double *Fl = R->SL->flux[R->dim];
  double *Fr = R->SR->flux[R->dim];
  R->U_hll = (double*) realloc(R->U_hll, nw * sizeof(double));
  R->F_hll = (double*) realloc(R->F_hll, nw * sizeof(double));
  for (int i=0; i<nw; ++i) {
    R->U_hll[i] = (ap*Ur[i] - am*Ul[i] +       (Fl[i] - Fr[i])) / (ap - am);
    R->F_hll[i] = (ap*Fl[i] - am*Fr[i] + ap*am*(Ur[i] - Ul[i])) / (ap - am);
  }
  return 0;
}

int _hll_sample(fluid_riemann *R, fluid_state *S, double s)
{
  double ap = R->ap;
  double am = R->am;
  double *Ul = R->SL->conserved;
  double *Ur = R->SR->conserved;
  double *Fl = R->SL->flux[R->dim];
  double *Fr = R->SR->flux[R->dim];
  double *U = S->conserved;
  double *F = S->flux[R->dim];
  int i;
  int nw = R->SL->nwaves;
  if      (         s<=am ) for (i=0; i<nw; ++i) U[i] = Ul[i];
  else if ( am<s && s<=ap ) for (i=0; i<nw; ++i) U[i] = R->U_hll[i];
  else if ( ap<s          ) for (i=0; i<nw; ++i) U[i] = Ur[i];  
  if      (         s<=am ) for (i=0; i<nw; ++i) F[i] = Fl[i];
  else if ( am<s && s<=ap ) for (i=0; i<nw; ++i) F[i] = R->F_hll[i];
  else if ( ap<s          ) for (i=0; i<nw; ++i) F[i] = Fr[i];
  return 0;
}

int _nrhyd_exact_exec(fluid_riemann *R)
{
  double *Pl = R->SL->primitive;
  double *Pr = R->SR->primitive;

  R->gm0 = R->SR->gammalawindex;
  R->gm1 = (R->gm0+1) / (2*R->gm0);
  R->gm2 = (R->gm0-1) / (2*R->gm0);
  R->gm3 = (R->gm0-1) / (R->gm0+1);
  R->gm4 = 2.0 / (R->gm0+1);
  R->gm5 = 2.0 / (R->gm0-1);

  R->pL = Pl[pre];
  R->pR = Pr[pre];
  R->uL = Pl[R->v1];
  R->uR = Pr[R->v1];

  R->AL = 2.0/((R->gm0+1)*Pl[rho]); // eqn 4.8
  R->AR = 2.0/((R->gm0+1)*Pr[rho]);
  R->BL = R->gm3*R->pL;
  R->BR = R->gm3*R->pR;

  R->aL = sqrt(R->gm0*Pl[pre] / Pl[rho]); // sound speed
  R->aR = sqrt(R->gm0*Pr[pre] / Pr[rho]);

  double p;
  int err = 1;

  if (err) {
    p = _soln_estimate(R, 0);
    err = _soln_solve(R, &p);
  }
  if (err) {
    p = _soln_estimate(R, 1);
    err = _soln_solve(R, &p);
  }
  if (err) {
    p = _soln_estimate(R, 2);
    err = _soln_solve(R, &p);
  }
  if (err == 0) {
    R->p_solution = p;
  }
  return err;
}

int _soln_solve(fluid_riemann *R, double *p)
{
  double f, g;
  int niter = 0;
  do {
    if (++niter > 25) {
      return FLUIDS_ERROR_RIEMANN;
    }
    f = _soln_f(R, *p);
    g = _soln_g(R, *p);
    *p -= f/g;
  } while (fabs(f) > 1e-12);
  return 0;
}

/* function whose root is needed */
double _soln_f(fluid_riemann *R, double p)
{
  double *Pl = R->SL->primitive;
  double *Pr = R->SR->primitive;
  double gm0 = R->gm0;
  double gm2 = R->gm2;
  double aL = R->aL;
  double aR = R->aR;
  double AL = R->AL;
  double AR = R->AR;
  double BL = R->BL;
  double BR = R->BR;
  double pL = R->pL;
  double pR = R->pR;

  double fL = (p > pL)             ?     // eqn 4.6
    (p - pL) * sqrt(AL / (p + BL)) :     // shock wave
    2*aL/(gm0-1) * (pow(p/pL, gm2) - 1); // rarefaction wave

  double fR = (p > pR)             ?     // eqn 4.7
    (p - pR) * sqrt(AR / (p + BR)) :     // shock wave
    2*aR/(gm0-1) * (pow(p/pR, gm2) - 1); // rarefaction wave

  return fL + fR + (Pr[R->v1] - Pl[R->v1]); // eqn 4.5
}

/* df/fx */
double _soln_g(fluid_riemann *R, double p)
{
  double *Pl = R->SL->primitive;
  double *Pr = R->SR->primitive;
  double gm1 = R->gm1;
  double aL = R->aL;
  double aR = R->aR;
  double AL = R->AL;
  double AR = R->AR;
  double BL = R->BL;
  double BR = R->BR;
  double pL = R->pL;
  double pR = R->pR;

  double fpL = (p > pL)                             ? // eqn 4.37
    sqrt(AL/(BL + p)) * (1 - 0.5*(p - pL)/(BL + p)) : // shock wave
    1.0/(Pl[rho] * aL) * pow(p/pL, -gm1);             // rarefaction wave

  double fpR = (p > pR)                             ? // eqn 4.37
    sqrt(AR/(BR + p)) * (1 - 0.5*(p - pR)/(BR + p)) : // shock wave
    1.0/(Pr[rho] * aR) * pow(p/pR, -gm1);             // rarefaction wave

  return fpL + fpR;
}

double _soln_estimate(fluid_riemann *R, int attempt)
{
  double *Pl = R->SL->primitive;
  double *Pr = R->SR->primitive;
  double gm0 = R->gm0;
  double aL = R->aL;
  double aR = R->aR;
  double pL = R->pL;
  double pR = R->pR;
  double uL = R->uL;
  double uR = R->uR;

  if (attempt == 0) { // eqn 4.46
    double gm1 = (gm0-1) / (2*gm0);
    double gm3 = 1.0/gm1;
    return pow(((aL+aR)+0.5*(gm0-1)*(uL-uR)) / (aL/pow(pL,gm1)+aR/pow(pR,gm1)),
	       gm3);
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

int _nrhyd_exact_sample(fluid_riemann *R, fluid_state *S, double s)
{
  /* Underscore after variable is Toro's (*), which indicates the Star Region */
  double p_ = R->p_solution;
  double *Pl = R->SL->primitive;
  double *Pr = R->SR->primitive;
  int v1 = R->v1;
  int v2 = R->v2;
  int v3 = R->v3;
  double gm0 = R->gm0;
  double gm1 = R->gm1;
  double gm2 = R->gm2;
  double gm3 = R->gm3;
  double gm4 = R->gm4;
  double gm5 = R->gm5;
  double aL = R->aL;
  double aR = R->aR;
  double pL = R->pL;
  double pR = R->pR;
  double uL = R->uL;
  double uR = R->uR;
  double AL = R->AL;
  double AR = R->AR;
  double BL = R->BL;
  double BR = R->BR;
  double aL_ = aL*pow(p_/pL, gm2); // eqn 4.54
  double aR_ = aR*pow(p_/pR, gm2);

  double fL = (p_ > pL)              ?    // eqn 4.6
    (p_ - pL) * sqrt(AL / (p_ + BL)) :    // shock wave
    2*aL/(gm0-1) * (pow(p_/pL, gm2) - 1); // rarefaction wave

  double fR = (p_ > pR)              ?    // eqn 4.7
    (p_ - pR) * sqrt(AR / (p_ + BR)) :    // shock wave
    2*aR/(gm0-1) * (pow(p_/pR, gm2) - 1); // rarefaction wave

  double u_ = 0.5*(uL + uR) - 0.5*(fL - fR); // eqn 4.88
  double SL  = uL - aL*sqrt(gm1 * p_/pL + gm2); // eqn 4.52
  double SHL = uL - aL;
  double STL = u_ - aL_;
  double STR = u_ + aR_;
  double SHR = uR + aR;
  double SR  = uR + aR*sqrt(gm1 * p_/pR + gm2); // eqn 4.59

  double FrhoL_ = Pl[rho]*pow(p_/pL, 1.0/gm0); // eqn 4.53 (fan)
  double FrhoR_ = Pr[rho]*pow(p_/pR, 1.0/gm0); // eqn 4.60
  double SrhoL_ = Pl[rho]*(p_/pL + gm3) / (gm3*p_/pL + 1); // eqn 4.50 (shock)
  double SrhoR_ = Pr[rho]*(p_/pR + gm3) / (gm3*p_/pR + 1); // eqn 4.57
  double *P = S->primitive;

  /* ---------------------------------------------------------------------------
     The logic below is the sampling procedure outlined by Toro in the flow
     chart Figure 4.14. It consists all-together of 10 cases: 6 regions (L/R,
     L/R fan, and L/R star) and then two sub cases for each of L/R fan and L/R
     star, corresponding to whether the L/R going wave is a fan or a shock.
     ------------------------------------------------------------------------ */
  /*                    |   |   |                                             */
  /*                    |   |   |                                             */
  if (s < u_) {         // sampling left of particle velocity characteristic
    if (p_ > pL) {          // left-moving wave is a shock
      if (s < SL) {             // left region
        memcpy(P, Pl, 5*sizeof(double));
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
      if (s < SHL) {            // left region
        memcpy(P, Pl, 5*sizeof(double));
      }
      else {
        if (s < STL) {          // left fan
          P[rho] = Pl[rho]*pow(gm4 + gm3*(uL - s)/aL, gm5); // eqn 4.56
          P[pre] = Pl[pre]*pow(gm4 + gm3*(uL - s)/aL, 1.0/gm2);
          P[v1 ] = gm4*(aL + uL/gm5 + s);
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
  /* ------------------------------------------------------------------------ */
  /*                    |   |   |                                             */
  /*                    |   |   |                                             */
  else {                // sampling right of particle velocity characteristic
    if (p_ > pR) {          // right-moving wave is a shock
      if (s > SR) {             // right region
        memcpy(P, Pr, 5*sizeof(double));
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
      if (s > SHR) {            // right region
        memcpy(P, Pr, 5*sizeof(double));
      }
      else {
        if (s > STR) {          // right fan
          P[rho] = Pr[rho]*pow(gm4 - gm3*(uR - s)/aR, gm5); // eqn 4.63
          P[pre] = Pr[pre]*pow(gm4 - gm3*(uR - s)/aR, 1.0/gm2);
          P[v1 ] = gm4*(-aR + uR/gm5 + s);
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
  fluids_p2c(S);
  fluids_update(S, FLUIDS_FLUX[R->dim]);
  return 0;
}
