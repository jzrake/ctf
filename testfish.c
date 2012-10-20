

#include <stdio.h>
#include <assert.h>
#include <math.h>
#define FISH_PRIVATE_DEFS
#include "fish.h"


#define asserteq(x,y) assert(fabs(x-y) < 1e-12)

// Passes when a trivial 1d problem produces the correct intercell fluxes
// -----------------------------------------------------------------------------
int test1()
{
  fish_state *S = fish_new();
  fish_setreconstruction(S, FISH_PLM);
  fish_setriemannsolver(S, FLUIDS_RIEMANN_EXACT);

  fluids_descr *D = fluids_descr_new();
  fluids_descr_setfluid(D, FLUIDS_GRAVP);
  fluids_descr_setgamma(D, 1.4);
  fluids_descr_seteos(D, FLUIDS_EOS_GAMMALAW);

  double P[5] = {1, 1, 1, 1, 1};
  double G[4] = {0, 0, 0, 0};
  fluids_state *fluid[100];

  double Fiph[500];
  for (int n=0; n<100; ++n) {
    fluid[n] = fluids_state_new();
    fluids_state_setdescr(fluid[n], D);
    fluids_state_setattr(fluid[n], P, FLUIDS_PRIMITIVE);
    fluids_state_setattr(fluid[n], G, FLUIDS_GRAVITY);
  }

  fish_intercellflux(S, fluid, Fiph, 100, 0);
  for (int n=1; n<100-2; ++n) {
    asserteq(Fiph[5*n], 1.0);
  }
  for (int n=0; n<100; ++n) {
    fluids_state_del(fluid[n]);
  }

  fluids_descr_del(D);
  fish_del(S);
  printf("TEST 1 PASSED\n");
  return 0;
}

int main()
{
  test1();
  return 0;
}
