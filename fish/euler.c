

#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "fish.h"

static struct fluids_state *fluid[100];
static fluids_descr *descr;
static double dx = 1.0 / 100.0;
static double dt = 0.001;

void init()
{
  descr = fluids_descr_new();
  fluids_descr_setfluid(descr, FLUIDS_NRHYD);
  fluids_descr_setgamma(descr, 1.4);
  fluids_descr_seteos(descr, FLUIDS_EOS_GAMMALAW);

  for (int n=0; n<100; ++n) {
    double P[5] = {0, 0, 0, 0, 0};
    //    double G[4] = {0, 0, 0, 0};
    if (n<50) {
      P[0] = 1.0;
      P[1] = 1.0;
    }
    else {
      P[0] = 0.125;
      P[1] = 0.1;
    }
    fluid[n] = fluids_state_new();
    fluids_state_setdescr(fluid[n], descr);
    fluids_state_setattr(fluid[n], P, FLUIDS_PRIMITIVE);
    //    fluids_state_setattr(fluid[n], G, FLUIDS_GRAVITY);
  }
}

void finish()
{
  for (int n=0; n<100; ++n) {
    fluids_state_del(fluid[n]);
  }
  fluids_descr_del(descr);
}

int timederiv(double *L)
{
  fish_state *S = fish_new();
  fish_setparami(S, FISH_WENO5, FISH_RECONSTRUCTION);
  fish_setparami(S, FLUIDS_RIEMANN_EXACT, FISH_RIEMANN_SOLVER);
  fish_setparami(S, FISH_SPECTRAL, FISH_SOLVER_TYPE);
  fish_setparamd(S, 2.0, FISH_PLM_THETA);

  for (int m=0; m < 5 * 100; ++m) {
    L[m] = 0.0;
  }
  int nzone = 100;
  fish_timederivative(S, fluid, 1, &nzone, &dx, L);

  /*
  double Fiph[500];
  fish_intercellflux(S, fluid, Fiph, 100, 0);

  for (int n=3; n<100-3; ++n) {
    for (int q=0; q<5; ++q) {
      L[5*n + q] = -(Fiph[5*n + q] - Fiph[5*(n-1) + q]) / dx;
    }
  }
  */

  for (int q=0; q<5; ++q) {
    L[5* 0 + q] = 0.0;
    L[5* 1 + q] = 0.0;
    L[5* 2 + q] = 0.0;
    L[5*97 + q] = 0.0;
    L[5*98 + q] = 0.0;
    L[5*99 + q] = 0.0;
  }
  fish_del(S);
  return 0;
}

int advance()
{
  double L[500], U0[500];

  for (int m=0; m < 5 * 100; ++m) { U0[m] = 0.0; }
  for (int n=0; n<100; ++n) {
    fluids_state_derive(fluid[n], &U0[5*n], FLUIDS_CONSERVED);
  }

  timederiv(L);
  for (int n=0; n<100; ++n) {
    double U[5];
    fluids_state_derive(fluid[n], U, FLUIDS_CONSERVED);
    for (int q=0; q<5; ++q) {
      U[q] += L[5*n + q] * dt * 0.5;
    }
    fluids_state_fromcons(fluid[n], U, FLUIDS_CACHE_DEFAULT);
  }

  timederiv(L);
  for (int n=0; n<100; ++n) {
    double U[5];
    fluids_state_derive(fluid[n], U, FLUIDS_CONSERVED);
    for (int q=0; q<5; ++q) {
      U[q] = U0[5*n + q] + L[5*n + q] * dt;
    }
    fluids_state_fromcons(fluid[n], U, FLUIDS_CACHE_DEFAULT);
  }
  return 0;
}

int fish_run_euler(double *result)
{
  double t = 0.0;
  init();
  for (int n=0; n<80; ++n) {
    clock_t start = clock();
    advance();
    t += dt;
    clock_t del = clock() - start;
    if (n % 20 == 0) {
      printf("n=%d t=%3.2f %f kz/s\n", n, t, 100.0 / (1e3*del / CLOCKS_PER_SEC));
    }
  }
  for (int n=0; n<100; ++n) {
    double P[5];
    fluids_state_getattr(fluid[n], P, FLUIDS_PRIMITIVE);
    memcpy(&result[5*n], P, 5 * sizeof(double));
  }
  finish();
  return 0;
}
