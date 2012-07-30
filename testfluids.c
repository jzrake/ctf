
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
  double gam = 1.4;
  double Gam;
  fluid_state *S = fluids_new();
  fluids_setfluid(S, FLUIDS_NRHYD);
  fluids_setattrib(S, &gam, FLUIDS_GAMMALAWINDEX);
  fluids_getattrib(S, &Gam, FLUIDS_GAMMALAWINDEX);
  fluids_setattrib(S, x, FLUIDS_PRIMITIVE);
  fluids_getattrib(S, y, FLUIDS_PRIMITIVE);
  fluids_del(S);
  asserteq(Gam, gam);
  for (int n=0; n<5; ++n) {
    asserteq(y[n], 1.0);
  }
  printf("TEST 1 PASSED\n");
  return 0;
}

// Passes when fluxes and sound speeds are computed correctly
// -----------------------------------------------------------------------------
int test2()
{
  double x[5] = {1, 1, 1, 1, 1};
  double F[5];
  double G[5];
  double gam = 1.4;
  double cs2;
  fluid_state *S = fluids_new();
  fluids_setfluid(S, FLUIDS_NRHYD);
  fluids_setattrib(S, &gam, FLUIDS_GAMMALAWINDEX);
  fluids_setattrib(S, x, FLUIDS_PRIMITIVE);
  fluids_p2c(S);
  fluids_update(S, FLUIDS_FLUX0 | FLUIDS_FLUX1 | FLUIDS_SOUNDSPEEDSQUARED);
  fluids_getattrib(S, F, FLUIDS_FLUX0);
  fluids_getattrib(S, G, FLUIDS_FLUX1);
  fluids_getattrib(S, &cs2, FLUIDS_SOUNDSPEEDSQUARED);
  fluids_del(S);
  asserteq(F[0], 1.0);
  asserteq(F[1], 5.0);
  asserteq(F[2], 2.0);
  asserteq(G[0], 1.0);
  asserteq(G[1], 5.0);
  asserteq(G[2], 1.0);
  asserteq(cs2, 1.4);
  printf("TEST 2 PASSED\n");
  return 0;
}

// Passes when sound speed calculation is skipped appropriately
// -----------------------------------------------------------------------------
int test3()
{
  double x[5] = {1, 1, 1, 1, 1};
  double gam = 1.4;
  double lam1[5], cs2;
  fluid_state *S = fluids_new();
  fluids_setfluid(S, FLUIDS_NRHYD);
  fluids_setattrib(S, &gam, FLUIDS_GAMMALAWINDEX);
  fluids_setattrib(S, x, FLUIDS_PRIMITIVE);
  fluids_p2c(S);
  fluids_update(S, FLUIDS_FLUX0 | FLUIDS_EVAL1);
  fluids_getattrib(S, lam1, FLUIDS_EVAL1);
  fluids_getattrib(S, &cs2, FLUIDS_SOUNDSPEEDSQUARED);
  fluids_del(S);
  asserteq(lam1[1], 1.0);
  asserteq(cs2, 0.0);
  printf("TEST 3 PASSED\n");
  return 0;
}

// Passes when
// (1) L.R = I
// (2) L.A.R = diag{ lam0, lam1, lam2, lam3, lam4 }
// -----------------------------------------------------------------------------
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
  fluids_setattrib(S, &gam, FLUIDS_GAMMALAWINDEX);
  fluids_setattrib(S, x, FLUIDS_PRIMITIVE);
  fluids_p2c(S);
  fluids_update(S,
                FLUIDS_EVAL0 |
                FLUIDS_LEVECS0 |
                FLUIDS_REVECS0 |
                FLUIDS_JACOBIAN0);
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

// Passes when the exact riemann solver creates a sane response to the trivial
// Riemann problem. Also writes a file containing thee solution to the Brio-Wu
// shocktube.
// -----------------------------------------------------------------------------
int test5()
{
  int solvers[2] = {FLUIDS_RIEMANN_HLL, FLUIDS_RIEMANN_EXACT};
  for (int solver=0; solver<2; ++solver) {
    double Pl[5] = {1, 1, 1, 1, 1};
    double Pr[5] = {1, 1, 1, 1, 1};
    double P_[5];
    double gam = 1.4;
    fluid_state *SL = fluids_new();
    fluid_state *SR = fluids_new();
    fluid_state *S_ = fluids_new();

    fluids_setfluid(SL, FLUIDS_NRHYD);
    fluids_setfluid(SR, FLUIDS_NRHYD);
    fluids_setfluid(S_, FLUIDS_NRHYD);

    fluids_setattrib(SL, &gam, FLUIDS_GAMMALAWINDEX);
    fluids_setattrib(SR, &gam, FLUIDS_GAMMALAWINDEX);
    fluids_setattrib(S_, &gam, FLUIDS_GAMMALAWINDEX);
    fluids_setattrib(SL, Pl, FLUIDS_PRIMITIVE);
    fluids_setattrib(SR, Pr, FLUIDS_PRIMITIVE);
    fluids_p2c(SL);
    fluids_p2c(SR);

    fluid_riemann *R = fluids_riemann_new();
    fluids_riemann_setsolver(R, solvers[solver]);
    fluids_riemann_setdim(R, 0);
    fluids_riemann_setstateL(R, SL);
    fluids_riemann_setstateR(R, SR);
    fluids_riemann_execute(R);
    fluids_riemann_sample(R, S_, 0.2);
    fluids_c2p(S_);
    fluids_getattrib(S_, P_, FLUIDS_PRIMITIVE);

    for (int n=0; n<5; ++n) {
      asserteq(P_[n], Pl[n]);
    }
    printf("TEST 5.%d PASSED\n", solver);

    fluids_riemann_del(R);
    fluids_del(SL);
    fluids_del(SR);
    fluids_del(S_);
  }
  return 0;
}

int main()
{
  test1();
  test2();
  test3();
  test4();
  test5();
  return 0;
}
