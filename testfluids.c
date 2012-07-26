
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "fluids.h"
#include "matrix.h"

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
  assert(Gam == gam);
  for (int n=0; n<5; ++n) {
    assert(y[n] == 1.0);
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
  assert(F[0] == 1.0);
  assert(F[1] == 5.0);
  assert(F[2] == 2.0);
  assert(G[0] == 1.0);
  assert(G[1] == 5.0);
  assert(G[2] == 1.0);
  assert(cs2 == 1.4);
  printf("F[0] = %f G[0] = %f\n", F[0], G[0]);
  printf("F[1] = %f G[1] = %f\n", F[1], G[1]);
  printf("F[2] = %f G[2] = %f\n", F[2], G[2]);
  printf("cs2 = %f\n", cs2);
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
  fluids_update(S, FLUIDS_FLUX0 | FLUIDS_EVALS1);
  fluids_getattrib(S, lam1, FLUIDS_EVALS1);
  fluids_getattrib(S, &cs2, FLUIDS_SOUNDSPEEDSQUARED);
  fluids_del(S);
  assert(lam1[1] == 1.0);
  assert(cs2 == 0.0);
  printf("TEST 3 PASSED\n");
  return 0;
}

// Passes when the eigenvector calculation goes correctly
// -----------------------------------------------------------------------------
int test4()
{
  double x[5] = {1, 1, 1, 1, 1};
  double gam = 1.4;
  double L[25];
  double R[25];
  double I[25];
  fluid_state *S = fluids_new();
  fluids_setfluid(S, FLUIDS_NRHYD);
  fluids_setattrib(S, &gam, FLUIDS_GAMMALAWINDEX);
  fluids_setattrib(S, x, FLUIDS_PRIMITIVE);
  fluids_p2c(S);
  fluids_update(S, FLUIDS_LEVECS0 | FLUIDS_REVECS0);
  fluids_getattrib(S, L, FLUIDS_LEVECS0);
  fluids_getattrib(S, R, FLUIDS_REVECS0);
  fluids_del(S);
  matrix_matrix_product(L, R, I, 5, 5, 5);
  for (int m=0; m<5; ++m) {
    for (int n=0; n<5; ++n) {
      assert(fabs(I[m*5 + n] - (double)(m==n)) < 1e-12);
      printf("(L.R)[%d][%d] = %f\n", m, n, I[m*5 + n]);      
    }
  }
  printf("TEST 4 PASSED\n");
  return 0;
}

int main()
{
  test1();
  test2();
  test3();
  test4();
  return 0;
}
