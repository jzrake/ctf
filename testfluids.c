
#include <stdio.h>
#include <assert.h>
#include <math.h>
#define FLUIDS_PRIVATE_DEFS
#include "fluids.h"
#include "matrix.h"

#define asserteq(x,y) assert(fabs(x-y) < 1e-12)

// Passes when get/set attributes work correctly
// -----------------------------------------------------------------------------
int test1()
{
  double x[5] = {1, 1, 1, 1, 1};
  double y[5];

  fluids_descr *D = fluids_descr_new();
  fluids_state *S = fluids_state_new();

  fluids_descr_setfluid(D, FLUIDS_NRHYD);
  fluids_descr_setgamma(D, 1.4);
  fluids_descr_seteos(D, FLUIDS_EOS_GAMMALAW);

  fluids_state_setdescr(S, D);

  fluids_state_setattr(S, x, FLUIDS_PRIMITIVE);
  fluids_state_getattr(S, y, FLUIDS_PRIMITIVE);

  fluids_descr_del(D);
  fluids_state_del(S);

  for (int n=0; n<5; ++n) {
    asserteq(y[n], 1.0);
  }
  printf("TEST 1 PASSED\n");
  return 0;
}

int main()
{
  test1();
  return 0;
}
