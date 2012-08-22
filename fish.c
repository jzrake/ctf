
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define FISH_PRIVATE_DEFS
#include "fish.h"


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

  for (int n=0; n<N-1; ++n) {
    fluids_riemann_setstateL(R, fluid[n]);
    fluids_riemann_setstateR(R, fluid[n+1]);
    fluids_riemann_execute(R);
    fluids_riemann_sample(R, S_, 0.0);
    fluids_getattrib(S_, &F[5*n], FLUIDS_FLUX0);
  }

  fluids_del(S_);
  fluids_riemann_del(R);
  return 0;
}
