
#include <stdio.h>
#include <assert.h>
#include <math.h>
#define FLUIDS_PRIVATE_DEFS
#include "fluids.h"
#include "matrix.h"

#define asserteq(x,y) assert(fabs(x-y) < 1e-12)

// Passes when get/set attributes on primitive struct work correctly
// -----------------------------------------------------------------------------
int test1()
{
  double x[5] = {1, 1, 1, 1, 1};
  double y[5];
  double gam;

  fluids_descr *D = fluids_descr_new();
  fluids_state *S = fluids_state_new();

  fluids_descr_setfluid(D, FLUIDS_NRHYD);
  fluids_descr_setgamma(D, 1.4);
  fluids_descr_seteos(D, FLUIDS_EOS_GAMMALAW);

  fluids_state_setdescr(S, D);
  fluids_state_setattr(S, x, FLUIDS_PRIMITIVE);
  fluids_state_getattr(S, y, FLUIDS_PRIMITIVE);
  fluids_descr_getgamma(D, &gam);

  fluids_state_del(S);
  fluids_descr_del(D);

  asserteq(1.4, gam);
  for (int n=0; n<5; ++n) {
    asserteq(y[n], 1.0);
  }
  printf("TEST 1 PASSED\n");
  return 0;
}

// Passes when derived quantities work correctly
// -----------------------------------------------------------------------------
int test2()
{
  double x[5] = {1, 1, 1, 1, 1};
  double U[5];

  fluids_descr *D = fluids_descr_new();
  fluids_state *S = fluids_state_new();

  fluids_descr_setfluid(D, FLUIDS_NRHYD);
  fluids_descr_setgamma(D, 1.4);
  fluids_descr_seteos(D, FLUIDS_EOS_GAMMALAW);

  fluids_state_setdescr(S, D);
  fluids_state_setattr(S, x, FLUIDS_PRIMITIVE);
  fluids_state_derive(S, U, FLUIDS_CONSERVED);

  asserteq(U[0], 1.0);
  asserteq(U[1], 4.0);
  asserteq(U[2], 1.0);
  asserteq(U[3], 1.0);
  asserteq(U[4], 1.0);

  fluids_state_del(S);
  fluids_descr_del(D);

  printf("TEST 2 PASSED\n");
  return 0;
}

int main()
{
  test1();
  test2();
  return 0;
}
