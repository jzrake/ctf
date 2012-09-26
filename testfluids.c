
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
  fluids_state_fromcons(S, U, FLUIDS_CACHE_ERASE);
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

// Passes when
// (1) L.R = I
// (2) L.A.R = diag{ lam0, lam1, lam2, lam3, lam4 }
// -----------------------------------------------------------------------------
/*
int test4()
{
  double x[5] = {1, 1, 1, 1, 1};
  double gam = 1.4;
  double V[5];
  double L[25];
  double R[25];
  double I[25];
  double A[25];
  double LA[25];
  double LAR[25];
  fluid_state *S = fluids_new();
  fluids_setfluid(S, FLUIDS_NRHYD);
  fluids_alloc(S, fields());
  fluids_setattrib(S, &gam, FLUIDS_GAMMALAWINDEX);
  fluids_setattrib(S, x, FLUIDS_PRIMITIVE);
  fluids_getattrib(S, V, FLUIDS_EVAL0);
  fluids_getattrib(S, L, FLUIDS_LEVECS0);
  fluids_getattrib(S, R, FLUIDS_REVECS0);
  fluids_getattrib(S, A, FLUIDS_JACOBIAN0);
  fluids_del(S);
  matrix_matrix_product(L, R, I, 5, 5, 5);
  matrix_matrix_product(L, A, LA, 5, 5, 5);
  matrix_matrix_product(LA, R, LAR, 5, 5, 5);
  for (int m=0; m<5; ++m) {
    for (int n=0; n<5; ++n) {
      asserteq(I[m*5 + n], (m==n));
      asserteq(LAR[m*5 + n], (m==n) * V[m]);
    }
  }
  printf("TEST 4 PASSED\n");
  return 0;
}
*/
int main()
{
  printf("sizeof(fluid_state) = %ld\n", sizeof(fluids_state));
  printf("sizeof(fluid_cache) = %ld\n", sizeof(fluids_cache));
  test1();
  test2();
  test3();
  return 0;
}
