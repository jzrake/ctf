
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define FISH_PRIVATE_DEFS
#include "fish.h"
#include "weno.h"

#define MAXQ 8 // for small statically-declared arrays

static void _pcm(fluids_state **src, fluids_state *L, fluids_state *R, long flag);
static void _plm(fluids_state **src, fluids_state *L, fluids_state *R, long flag);
static void _weno5(fluids_state **src, fluids_state *L, fluids_state *R, int Q);
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
int fish_getfluid(fish_state *S, int *fluid)
{
  *fluid = S->fluid;
  return 0;
}
int fish_getriemannsolver(fish_state *S, int *riemannsolver)
{
  *riemannsolver = S->riemannsolver;
  return 0;
}
int fish_getreconstruction(fish_state *S, int *reconstruction)
{
  *reconstruction = S->reconstruction;
  return 0;
}
int fish_getplmtheta(fish_state *S, double *plmtheta)
{
  *plmtheta = S->plmtheta;
  return 0;
}

int fish_intercellflux(fish_state *S, fluids_state **fluid, double *F, int N,
                       int dim)
{
  fluids_descr *D;
  fluids_state_getdescr(fluid[0], &D);

  int Q = fluids_descr_getncomp(D, FLUIDS_PRIMITIVE);
  fluids_state *S_ = fluids_state_new();
  fluids_state *SL = fluids_state_new();
  fluids_state *SR = fluids_state_new();

  /* Assumes all states have the same descriptor, after all what sense does this
     make otherwise? */
  fluids_state_setdescr(S_, D);
  fluids_state_setdescr(SL, D);
  fluids_state_setdescr(SR, D);
  fluids_state_cache(S_, FLUIDS_CACHE_CREATE);
  fluids_state_cache(SL, FLUIDS_CACHE_CREATE);
  fluids_state_cache(SR, FLUIDS_CACHE_CREATE);

  fluids_riemn *R = fluids_riemn_new();
  fluids_riemn_setsolver(R, S->riemannsolver);
  fluids_riemn_setdim(R, dim);
  fluids_riemn_setstateL(R, SL);
  fluids_riemn_setstateR(R, SR);

  switch (S->reconstruction) {
  case FISH_NONE:
    for (int n=0; n<N-1; ++n) {
      _pcm(&fluid[n], SL, SR, FLUIDS_PRIMITIVE);
      _pcm(&fluid[n], SL, SR, FLUIDS_GRAVITY);
      fluids_riemn_execute(R);
      fluids_riemn_sample(R, S_, 0.0);
      fluids_state_derive(S_, &F[Q*n], FLUIDS_FLUX[dim]);
    }
    break;
  case FISH_PLM:
    reconstruct_set_plm_theta(S->plmtheta);
    for (int n=1; n<N-2; ++n) {
      _plm(&fluid[n], SL, SR, FLUIDS_PRIMITIVE);
      _plm(&fluid[n], SL, SR, FLUIDS_GRAVITY);
      fluids_riemn_execute(R);
      fluids_riemn_sample(R, S_, 0.0);
      fluids_state_derive(S_, &F[Q*n], FLUIDS_FLUX[dim]);
    }
    break;
  case FISH_WENO5:
    for (int n=2; n<N-3; ++n) {
      /*
      _weno5(&fluid[n], SL, SR, Q);
      fluids_riemn_execute(R);
      fluids_riemn_sample(R, S_, 0.0);
      fluids_state_derive(S_, &F[Q*n], FLUIDS_FLUX[dim]);
      */
    }
  default:
    break;
  }

  fluids_state_del(SL);
  fluids_state_del(SR);
  fluids_state_del(S_);
  fluids_riemn_del(R);
  return 0;
}

void _pcm(fluids_state **src, fluids_state *L, fluids_state *R, long flag)
{
  double Pl[MAXQ], Pr[MAXQ];
  fluids_descr *D;
  int fluid;
  fluids_state_getdescr(src[0], &D);
  fluids_descr_getfluid(D, &fluid);
  fluids_state_getattr(src[0], Pl, flag);
  fluids_state_getattr(src[1], Pr, flag);
  fluids_state_setattr(L, Pl, flag);
  fluids_state_setattr(R, Pr, flag);
}

void _plm(fluids_state **src, fluids_state *L, fluids_state *R, long flag)
{
  double P[4][MAXQ];
  double Pl[MAXQ], Pr[MAXQ];
  fluids_descr *D;
  int fluid;
  fluids_state_getdescr(src[0], &D);
  fluids_descr_getfluid(D, &fluid);
  int Q = fluids_descr_getncomp(D, flag);
  for (int j=-1; j<3; ++j) {
    fluids_state_getattr(src[j], P[j+1], flag);
  }
  for (int q=0; q<Q; ++q) {
    double v[4];
    for (int j=0; j<4; ++j) {
      v[j] = P[j][q];
    }
    Pl[q] = reconstruct(&v[1], PLM_C2R);
    Pr[q] = reconstruct(&v[2], PLM_C2L);
  }
  fluids_state_setattr(L, Pl, flag);
  fluids_state_setattr(R, Pr, flag);
}

void _weno5(fluids_state **src, fluids_state *L, fluids_state *R, int Q)
{
  double P[6][MAXQ];
  double Pl[MAXQ], Pr[MAXQ];
  double Gl[MAXQ], Gr[MAXQ];
  fluids_descr *D;
  int fluid;
  fluids_state_getdescr(src[0], &D);
  fluids_descr_getfluid(D, &fluid);
  for (int j=-2; j<4; ++j) {
    fluids_state_getattr(src[j], P[j+2], FLUIDS_PRIMITIVE);
  }
  for (int q=0; q<Q; ++q) {
    double v[6];
    for (int j=0; j<6; ++j) {
      v[j] = P[j][q];
    }
    Pl[q] = reconstruct(&v[2], WENO5_FD_C2R);
    Pr[q] = reconstruct(&v[3], WENO5_FD_C2L);
  }
  fluids_state_setattr(L, Pl, FLUIDS_PRIMITIVE);
  fluids_state_setattr(R, Pr, FLUIDS_PRIMITIVE);
  if (fluid == FLUIDS_GRAVS ||
      fluid == FLUIDS_GRAVP ||
      fluid == FLUIDS_GRAVE) {
    fluids_state_getattr(src[2], Gl, FLUIDS_GRAVITY);
    fluids_state_getattr(src[3], Gr, FLUIDS_GRAVITY);
    fluids_state_setattr(L, Gl, FLUIDS_GRAVITY);
    fluids_state_setattr(R, Gr, FLUIDS_GRAVITY);
  }
}
