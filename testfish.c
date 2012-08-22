

#include <stdio.h>
#include <assert.h>
#include <math.h>
#define FISH_PRIVATE_DEFS
#include "fish.h"


#define asserteq(x,y) assert(fabs(x-y) < 1e-12)

static long fields()
{
  long modes = 0;
  modes |= FLUIDS_CONSERVED;
  modes |= FLUIDS_PRIMITIVE;
  modes |= FLUIDS_FLUXALL;
  modes |= FLUIDS_EVALSALL;
  modes |= FLUIDS_LEVECSALL;
  modes |= FLUIDS_REVECSALL;
  modes |= FLUIDS_JACOBIANALL;
  return modes;
}

// Passes when a trivial 1d problem produces the correct intercell fluxes
// -----------------------------------------------------------------------------
int test1()
{
  fish_state *S = fish_new();
  double P[5] = {1, 1, 1, 1, 1};
  double gam = 1.4;
  fluid_state *fluid[100];
  double Fiph[500];
  for (int n=0; n<100; ++n) {
    fluid[n] = fluids_new();
    fluids_setfluid(fluid[n], FLUIDS_NRHYD);
    fluids_alloc(fluid[n], fields());
    fluids_setattrib(fluid[n], &gam, FLUIDS_GAMMALAWINDEX);
    fluids_setattrib(fluid[n], P, FLUIDS_PRIMITIVE);
  }
  fish_intercellflux(S, fluid, Fiph, 100);
  for (int n=0; n<100-1; ++n) {
    asserteq(Fiph[5*n], 1.0);
  }
  for (int n=0; n<100; ++n) {
    fluids_del(fluid[n]);
  }
  fish_del(S);
  printf("TEST 1 PASSED\n");
  return 0;
}

int main()
{
  test1();
  return 0;
}
