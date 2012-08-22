
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define FISH_PRIVATE_DEFS
#include "fish.h"
#include "weno.h"

#define MAXQ 8 // for small statically-declared arrays

static void _plm(fluid_state **src, fluid_state *L, fluid_state *R, int Q);
static void _weno5(fluid_state **src, fluid_state *L, fluid_state *R, int Q);
static const long FLUIDS_FLUX[3] = {FLUIDS_FLUX0, FLUIDS_FLUX1, FLUIDS_FLUX2};
static const long FLUIDS_EVAL[3] = {FLUIDS_EVAL0, FLUIDS_EVAL1, FLUIDS_EVAL2};

fish_state *fish_new(void)
{
  fish_state *S = (fish_state*) malloc(sizeof(fish_state));
  fish_state state = {
    .fluid = FLUIDS_NRHYD,
    .riemannsolver = FLUIDS_RIEMANN_HLL,
    .reconstruction = FISH_PLM,
    .plmtheta = 2.0,
  } ;
  *S = state;
  return S;
}

int fish_del(fish_state *S)
{
  free(S);
  return 0;
}

int fish_setfluid(fish_state *S, int fluid)
{
  S->fluid = fluid;
  return 0;
}

int fish_setriemannsolver(fish_state *S, int riemannsolver)
{
  S->riemannsolver = riemannsolver;
  return 0;
}

int fish_setreconstruction(fish_state *S, int reconstruction)
{
  S->reconstruction = reconstruction;
  return 0;
}

int fish_setplmtheta(fish_state *S, double plmtheta)
{
  S->plmtheta = plmtheta;
  return 0;
}

int fish_intercellflux(fish_state *S, fluid_state **fluid, double *F, int N,
		       int dim)
{
  int Q = fluids_getnwaves(S->fluid);
  fluid_state *S_ = fluids_new();
  fluids_setfluid(S_, S->fluid);
  fluids_alloc(S_, FLUIDS_CONSERVED | FLUIDS_PRIMITIVE | FLUIDS_FLUX[dim]);

  fluid_riemann *R = fluids_riemann_new();
  fluids_riemann_setsolver(R, S->riemannsolver);
  fluids_riemann_setdim(R, 0);

  fluid_state *SL = fluids_new();
  fluid_state *SR = fluids_new();
  fluids_setfluid(SL, S->fluid);
  fluids_setfluid(SR, S->fluid);
  fluids_alloc(SL, FLUIDS_FLAGSALL);
  fluids_alloc(SR, FLUIDS_FLAGSALL);

  switch (S->reconstruction) {
  case FISH_NONE:
    for (int n=0; n<N-1; ++n) {
      fluids_riemann_setstateL(R, fluid[n]);
      fluids_riemann_setstateR(R, fluid[n+1]);
      fluids_riemann_execute(R);
      fluids_riemann_sample(R, S_, 0.0);
      fluids_getattrib(S_, &F[Q*n], FLUIDS_FLUX[dim]);
    }
    break;
  case FISH_PLM:
    reconstruct_set_plm_theta(S->plmtheta);
    for (int n=1; n<N-2; ++n) {
      _plm(&fluid[n], SL, SR, Q);
      fluids_riemann_setstateL(R, SL);
      fluids_riemann_setstateR(R, SR);
      fluids_riemann_execute(R);
      fluids_riemann_sample(R, S_, 0.0);
      fluids_getattrib(S_, &F[Q*n], FLUIDS_FLUX[dim]);
    }
    break;
  case FISH_WENO5:
    for (int n=2; n<N-3; ++n) {
      _weno5(&fluid[n], SL, SR, Q);
      fluids_riemann_setstateL(R, SL);
      fluids_riemann_setstateR(R, SR);
      fluids_riemann_execute(R);
      fluids_riemann_sample(R, S_, 0.0);
      fluids_getattrib(S_, &F[Q*n], FLUIDS_FLUX[dim]);
    }
  default:
    break;
  }

  fluids_del(SL);
  fluids_del(SR);
  fluids_del(S_);
  fluids_riemann_del(R);
  return 0;
}

void _plm(fluid_state **src, fluid_state *L, fluid_state *R, int Q)
{
  double P0[MAXQ], P1[MAXQ], P2[MAXQ], P3[MAXQ];
  double Pl[MAXQ], Pr[MAXQ];

  fluids_getattrib(src[-1], P0, FLUIDS_PRIMITIVE);
  fluids_getattrib(src[ 0], P1, FLUIDS_PRIMITIVE);
  fluids_getattrib(src[ 1], P2, FLUIDS_PRIMITIVE);
  fluids_getattrib(src[ 2], P3, FLUIDS_PRIMITIVE);

  for (int q=0; q<Q; ++q) {
    double v[4] = {P0[q], P1[q], P2[q], P3[q]};
    Pl[q] = reconstruct(&v[1], PLM_C2R);
    Pr[q] = reconstruct(&v[2], PLM_C2L);
  }

  fluids_setattrib(L, Pl, FLUIDS_PRIMITIVE);
  fluids_setattrib(R, Pr, FLUIDS_PRIMITIVE);
}


void _weno5(fluid_state **src, fluid_state *L, fluid_state *R, int Q)
{
  double P[6][MAXQ];
  double Pl[MAXQ], Pr[MAXQ];

  for (int j=-2; j<4; ++j) {
    fluids_getattrib(src[j], P[j+2], FLUIDS_PRIMITIVE);
  }

  for (int q=0; q<Q; ++q) {
    double v[6];
    for (int j=0; j<6; ++j) {
      v[j] = P[j][q];
    }
    Pl[q] = reconstruct(&v[2], WENO5_FD_C2R);
    Pr[q] = reconstruct(&v[3], WENO5_FD_C2L);
  }

  fluids_setattrib(L, Pl, FLUIDS_PRIMITIVE);
  fluids_setattrib(R, Pr, FLUIDS_PRIMITIVE);
}
