
#include <stdio.h>
#include <assert.h>
#include <math.h>
#define FLUIDS_PRIVATE_DEFS
#include "fluids.h"
#include "matrix.h"

#define asserteq(x,y) assert(fabs(x-y) < 1e-12)
static long fluidtype = FLUIDS_NRHYD;

// Passes when get/set attributes on primitive struct work correctly
// -----------------------------------------------------------------------------
int test1()
{
  double P[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
  double G[4] = { 1.0, 1.0, 1.0, 1.0 };
  double Q[5];
  double R[5];
  double gam;

  fluids_descr *D = fluids_descr_new();
  fluids_state *S = fluids_state_new();

  fluids_descr_setfluid(D, fluidtype);
  fluids_descr_setgamma(D, 1.4);
  fluids_descr_seteos(D, FLUIDS_EOS_GAMMALAW);

  fluids_state_setdescr(S, D);
  fluids_state_setattr(S, P, FLUIDS_PRIMITIVE);
  fluids_state_getattr(S, Q, FLUIDS_PRIMITIVE);
  if (fluidtype == FLUIDS_GRAVP) {
    fluids_state_setattr(S, G, FLUIDS_GRAVITY);    
  }
  fluids_descr_getgamma(D, &gam);

  asserteq(1.4, gam);
  for (int n=0; n<5; ++n) {
    asserteq(Q[n], 1.0);
  }

  fluids_state_mapbuffer(S, R, FLUIDS_PRIMITIVE);
  R[0] = 2.0;
  R[1] = 2.0;
  R[2] = 2.0;
  R[3] = 2.0;
  R[4] = 2.0;
  fluids_state_getattr(S, Q, FLUIDS_PRIMITIVE);
  for (int n=0; n<5; ++n) {
    asserteq(Q[n], 2.0);
  }

  fluids_state_del(S);
  fluids_descr_del(D);

  printf("TEST 1 PASSED\n");
  return 0;
}

// Passes when derived quantities work correctly
// -----------------------------------------------------------------------------
int test2()
{
  double P[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
  double G[4] = { 1.0, 1.0, 1.0, 1.0 };
  double U[5], F[5];

  fluids_descr *D = fluids_descr_new();
  fluids_state *S = fluids_state_new();

  fluids_descr_setfluid(D, fluidtype);
  fluids_descr_setgamma(D, 1.4);
  fluids_descr_seteos(D, FLUIDS_EOS_GAMMALAW);

  fluids_state_setdescr(S, D);
  fluids_state_setattr(S, P, FLUIDS_PRIMITIVE);
  if (fluidtype == FLUIDS_GRAVP) {
    fluids_state_setattr(S, G, FLUIDS_GRAVITY);    
  }
  fluids_state_derive(S, F, FLUIDS_FLUX0);
  fluids_state_derive(S, U, FLUIDS_CONSERVED);

  asserteq(U[0], 1.0);
  asserteq(U[1], 4.0);
  asserteq(U[2], 1.0);
  asserteq(U[3], 1.0);
  asserteq(U[4], 1.0);

  asserteq(F[0], 1.0);
  asserteq(F[1], 5.0);
  if (fluidtype != FLUIDS_GRAVP) {
    asserteq(F[2], 2.0);
    asserteq(F[3], 1.0);
    asserteq(F[4], 1.0);
  }
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

  fluids_descr_setfluid(D, fluidtype);
  fluids_descr_setgamma(D, 1.4);
  fluids_descr_seteos(D, FLUIDS_EOS_GAMMALAW);

  fluids_state_setdescr(S, D);
  fluids_state_fromcons(S, U, FLUIDS_CACHE_DEFAULT);
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
  asserteq(fluids_descr_getncomp(D, FLUIDS_PRIMITIVE), 5);

  fluids_state_del(S);
  fluids_descr_del(D);

  printf("TEST 3 PASSED\n");
  return 0;
}

// Passes when the cache logic over many states works correctly
// -----------------------------------------------------------------------------
int test4()
{
  double U[5] = { 1.0, 4.0, 1.0, 1.0, 1.0 };
  double P[5];
  fluids_descr *D = fluids_descr_new();
  fluids_state *S[3];

  fluids_descr_setfluid(D, fluidtype);
  fluids_descr_setgamma(D, 1.4);
  fluids_descr_seteos(D, FLUIDS_EOS_GAMMALAW);

  for (int n=0; n<3; ++n) {
    U[4] = n;
    S[n] = fluids_state_new();
    fluids_state_setdescr(S[n], D);
    fluids_state_cache(S[n], FLUIDS_CACHE_CREATE);
    fluids_state_fromcons(S[n], U, FLUIDS_CACHE_DEFAULT);
  }
  for (int n=0; n<3; ++n) {
    fluids_state_getattr(S[n], P, FLUIDS_PRIMITIVE);
    fluids_state_derive(S[n], U, FLUIDS_CONSERVED);
    asserteq(P[4], n);
    asserteq(U[4], n);
  }

  for (int n=0; n<3; ++n) {
    fluids_state_cache(S[n], FLUIDS_CACHE_ERASE);
  }
  for (int n=0; n<3; ++n) {
    fluids_state_getattr(S[n], P, FLUIDS_PRIMITIVE);
    fluids_state_derive(S[n], U, FLUIDS_CONSERVED);
    asserteq(P[4], n);
    asserteq(U[4], n);
  }

  for (int n=0; n<3; ++n) {
    fluids_state_del(S[n]);
  }
  fluids_descr_del(D);

  printf("TEST 4 PASSED\n");
  return 0;
}

// Passes when
// (1) L.R = I
// (2) L.A.R = diag{ lam0, lam1, lam2, lam3, lam4 }
// -----------------------------------------------------------------------------
int test5()
{
  double P[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
  double V[5], V_[5];
  double L[25];
  double R[25];
  double A[25];
  double I[25];
  double LA[25];
  double LAR[25];

  fluids_descr *D = fluids_descr_new();
  fluids_state *S = fluids_state_new();

  fluids_descr_setfluid(D, fluidtype);
  fluids_descr_setgamma(D, 1.4);
  fluids_descr_seteos(D, FLUIDS_EOS_GAMMALAW);

  fluids_state_setdescr(S, D);
  fluids_state_setattr(S, P, FLUIDS_PRIMITIVE);
  fluids_state_derive(S, V, FLUIDS_EVAL0);
  fluids_state_derive(S, L, FLUIDS_LEVECS0);
  fluids_state_derive(S, R, FLUIDS_REVECS0);
  fluids_state_derive(S, A, FLUIDS_JACOBIAN0);

  fluids_state_del(S);
  fluids_descr_del(D);

  matrix_matrix_product(L, R, I, 5, 5, 5);
  matrix_matrix_product(L, A, LA, 5, 5, 5);
  matrix_matrix_product(LA, R, LAR, 5, 5, 5);

  V_[0] = V[0]; // permute eigenvalues to agree with Jacobian
  V_[1] = V[4];
  V_[2] = V[1];
  V_[3] = V[2];
  V_[4] = V[3];

  for (int m=0; m<5; ++m) {
    for (int n=0; n<5; ++n) {
      asserteq(I[m*5 + n], (m==n));
      asserteq(LAR[m*5 + n], (m==n) * V_[m]);
    }
  }

  printf("TEST 5 PASSED\n");
  return 0;
}

// Passes when the exact riemann solver creates a sane response to the trivial
// Riemann problem.
// -----------------------------------------------------------------------------
int test6()
{
  int solvers[3] = {FLUIDS_RIEMANN_HLL,
                    FLUIDS_RIEMANN_HLLC,
                    FLUIDS_RIEMANN_EXACT};
  for (int solver=0; solver<3; ++solver) {
    double Pl[5] = {1, 1, 1, 1, 1};
    double Pr[5] = {1, 1, 1, 1, 1};
    double Ul[5];
    double P_[5];
    double U_[5];
    fluids_descr *D = fluids_descr_new();
    fluids_state *SL = fluids_state_new();
    fluids_state *SR = fluids_state_new();
    fluids_state *S_ = fluids_state_new();

    fluids_descr_setfluid(D, fluidtype);
    fluids_descr_setgamma(D, 1.4);
    fluids_descr_seteos(D, FLUIDS_EOS_GAMMALAW);

    fluids_state_setdescr(SL, D);
    fluids_state_setdescr(SR, D);
    fluids_state_setdescr(S_, D);

    fluids_state_setattr(SL, Pl, FLUIDS_PRIMITIVE);
    fluids_state_setattr(SR, Pr, FLUIDS_PRIMITIVE);

    fluids_riemn *R = fluids_riemn_new();
    fluids_riemn_setsolver(R, solvers[solver]);

    fluids_riemn_setdim(R, 0);
    fluids_riemn_setstateL(R, SL);
    fluids_riemn_setstateR(R, SR);
    fluids_riemn_execute(R);
    fluids_riemn_sample(R, S_, 0.2);
    fluids_state_getattr(S_, P_, FLUIDS_PRIMITIVE);
    fluids_state_derive(S_, U_, FLUIDS_CONSERVED);
    fluids_state_derive(SL, Ul, FLUIDS_CONSERVED);

    for (int n=0; n<5; ++n) {
      asserteq(P_[n], Pl[n]);
      asserteq(U_[n], Ul[n]);
    }
    printf("TEST 6.%d PASSED\n", solver);

    fluids_riemn_del(R);
    fluids_state_del(SL);
    fluids_state_del(SR);
    fluids_state_del(S_);
    fluids_descr_del(D);
  }
  return 0;
}

int main()
{
  printf("sizeof(fluid_state) = %ld\n", sizeof(fluids_state));
  printf("sizeof(fluid_cache) = %ld\n", sizeof(fluids_cache));
  test1();
  test2();
  if (fluidtype == FLUIDS_NRHYD) {
    test3();
    test4();
    test5();
    test6();
  }
  return 0;
}
