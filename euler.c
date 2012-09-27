

#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
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
  fish_setfluid(S, FLUIDS_NRHYD);
  fish_setriemannsolver(S, FLUIDS_RIEMANN_HLLC);
  fish_setreconstruction(S, FISH_PLM);
  fish_setplmtheta(S, 2.0);

  double Fiph[500];
  fish_intercellflux(S, fluid, Fiph, 100, 0);

  for (int n=3; n<100-3; ++n) {
    for (int q=0; q<5; ++q) {
      L[5*n + q] = -(Fiph[5*n + q] - Fiph[5*(n-1) + q]) / dx;
    }
  }

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
    fluids_state_fromcons(fluid[n], U, FLUIDS_CACHE_RESET);
  }

  timederiv(L);
  for (int n=0; n<100; ++n) {
    double U[5];
    fluids_state_derive(fluid[n], U, FLUIDS_CONSERVED);
    for (int q=0; q<5; ++q) {
      U[q] = U0[5*n + q] + L[5*n + q] * dt;
    }
    fluids_state_fromcons(fluid[n], U, FLUIDS_CACHE_RESET);
  }
  return 0;
}

int main()
{
  init();
  for (int n=0; n<200; ++n) {
    clock_t start = clock();
    advance();
    clock_t del = clock() - start;
    printf("running at %f kz/s\n", 100.0 / (1e3*del / CLOCKS_PER_SEC));
  }
  FILE *outf = fopen("euler.dat", "w");
  for (int n=0; n<100; ++n) {
    double P[5];
    fluids_state_getattr(fluid[n], P, FLUIDS_PRIMITIVE);
    fprintf(outf, "%d %f %f %f\n", n, P[0], P[1], P[2]);
  }
  fclose(outf);
  finish();
  return 0;
}
