
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
  double P[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
  double Q[5];
  double gam;

  fluids_descr *D = fluids_descr_new();
  fluids_state *S = fluids_state_new();

  fluids_descr_setfluid(D, FLUIDS_NRHYD);
  fluids_descr_setgamma(D, 1.4);
  fluids_descr_seteos(D, FLUIDS_EOS_GAMMALAW);

  fluids_state_setdescr(S, D);
  fluids_state_setattr(S, P, FLUIDS_PRIMITIVE);
  fluids_state_getattr(S, Q, FLUIDS_PRIMITIVE);
  fluids_descr_getgamma(D, &gam);

  fluids_state_del(S);
  fluids_descr_del(D);

  asserteq(1.4, gam);
  for (int n=0; n<5; ++n) {
    asserteq(Q[n], 1.0);
  }
  printf("TEST 1 PASSED\n");
  return 0;
}

// Passes when derived quantities work correctly
// -----------------------------------------------------------------------------
int test2()
{
  double P[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
  double U[5], F[5];

  fluids_descr *D = fluids_descr_new();
  fluids_state *S = fluids_state_new();

  fluids_descr_setfluid(D, FLUIDS_NRHYD);
  fluids_descr_setgamma(D, 1.4);
  fluids_descr_seteos(D, FLUIDS_EOS_GAMMALAW);

  fluids_state_setdescr(S, D);
  fluids_state_setattr(S, P, FLUIDS_PRIMITIVE);
  fluids_state_derive(S, F, FLUIDS_FLUX0);
  fluids_state_derive(S, U, FLUIDS_CONSERVED);

  asserteq(U[0], 1.0);
  asserteq(U[1], 4.0);
  asserteq(U[2], 1.0);
  asserteq(U[3], 1.0);
  asserteq(U[4], 1.0);

  asserteq(F[0], 1.0);
  asserteq(F[1], 5.0);
  asserteq(F[2], 2.0);
  asserteq(F[3], 1.0);
  asserteq(F[4], 1.0);

  fluids_state_del(S);
  fluids_descr_del(D);

  printf("TEST 2 PASSED\n");
  return 0;
}

// Passes when the fromcons function works correctly
// -----------------------------------------------------------------------------
int test3()
{
  double U[5] = { 1.0, 4.0, 1.0, 1.0, 1.0 };
  double F[5], P[5];

  fluids_descr *D = fluids_descr_new();
  fluids_state *S = fluids_state_new();

  fluids_descr_setfluid(D, FLUIDS_NRHYD);
  fluids_descr_setgamma(D, 1.4);
  fluids_descr_seteos(D, FLUIDS_EOS_GAMMALAW);

  fluids_state_setdescr(S, D);
  fluids_state_fromcons(S, U);
  fluids_state_derive(S, F, FLUIDS_FLUX0);
  fluids_state_getattr(S, P, FLUIDS_PRIMITIVE);

  asserteq(P[0], 1.0);
  asserteq(P[1], 1.0);
  asserteq(P[2], 1.0);
  asserteq(P[3], 1.0);
  asserteq(P[4], 1.0);

  asserteq(F[0], 1.0);
  asserteq(F[1], 5.0);
  asserteq(F[2], 2.0);
  asserteq(F[3], 1.0);
  asserteq(F[4], 1.0);

  fluids_state_del(S);
  fluids_descr_del(D);

  printf("TEST 3 PASSED\n");
  return 0;
}

int main()
{
  test1();
  test2();
  test3();
  return 0;
}
