
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define FISH_PRIVATE_DEFS
#include "fish.h"
#include "weno.h"

static void _reconstr(fluid_state **src, fluid_state *L, fluid_state *R);

fish_state *fish_new(void)
{
  fish_state *S = (fish_state*) malloc(sizeof(fish_state));
  fish_state state = {

  } ;
  *S = state;
  return S;
}

int fish_del(fish_state *S)
{
  free(S);
  return 0;
}

int fish_intercellflux(fish_state *S, fluid_state **fluid, double *F, int N)
{
  fluid_state *S_ = fluids_new();
  fluids_setfluid(S_, FLUIDS_NRHYD);
  fluids_alloc(S_, FLUIDS_CONSERVED | FLUIDS_PRIMITIVE | FLUIDS_FLUX0);

  fluid_riemann *R = fluids_riemann_new();
  fluids_riemann_setsolver(R, FLUIDS_RIEMANN_HLL);
  fluids_riemann_setdim(R, 0);

  if (0) {
    for (int n=0; n<N-1; ++n) {
      fluids_riemann_setstateL(R, fluid[n]);
      fluids_riemann_setstateR(R, fluid[n+1]);
      fluids_riemann_execute(R);
      fluids_riemann_sample(R, S_, 0.0);
      fluids_getattrib(S_, &F[5*n], FLUIDS_FLUX0);
    }
  }
  else {
    fluid_state *SL = fluids_new();
    fluid_state *SR = fluids_new();
    fluids_setfluid(SL, FLUIDS_NRHYD);
    fluids_setfluid(SR, FLUIDS_NRHYD);
    fluids_alloc(SL, FLUIDS_FLAGSALL);
    fluids_alloc(SR, FLUIDS_FLAGSALL);
    for (int n=1; n<N-2; ++n) {
      _reconstr(&fluid[n], SL, SR);
      fluids_riemann_setstateL(R, SL);
      fluids_riemann_setstateR(R, SR);
      fluids_riemann_execute(R);
      fluids_riemann_sample(R, S_, 0.0);
      fluids_getattrib(S_, &F[5*n], FLUIDS_FLUX0);
    }
    fluids_del(SL);
    fluids_del(SR);
  }
  fluids_del(S_);
  fluids_riemann_del(R);
  return 0;
}

void _reconstr(fluid_state **src, fluid_state *L, fluid_state *R)
{
  double P0[5], P1[5], P2[5], P3[5];
  double Pl[5], Pr[5];

  fluids_getattrib(src[-1], P0, FLUIDS_PRIMITIVE);
  fluids_getattrib(src[ 0], P1, FLUIDS_PRIMITIVE);
  fluids_getattrib(src[ 1], P2, FLUIDS_PRIMITIVE);
  fluids_getattrib(src[ 2], P3, FLUIDS_PRIMITIVE);

  for (int q=0; q<5; ++q) {
    double v[4] = {P0[q], P1[q], P2[q], P3[q]};
    Pl[q] = reconstruct(&v[1], PLM_C2R);
    Pr[q] = reconstruct(&v[2], PLM_C2L);
  }

  fluids_setattrib(L, Pl, FLUIDS_PRIMITIVE);
  fluids_setattrib(R, Pr, FLUIDS_PRIMITIVE);
}
