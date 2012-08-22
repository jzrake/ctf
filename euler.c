

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "fish.h"

static struct fluid_state *fluid[100];
static double AdiabaticGamma = 1.4;
static double dx = 1.0 / 100.0;
static double dt = 0.001;

void init()
{
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
    fluid[n] = fluids_new();
    fluids_setfluid(fluid[n], FLUIDS_NRHYD);
    fluids_alloc(fluid[n], FLUIDS_FLAGSALL);
    fluids_setattrib(fluid[n], &AdiabaticGamma, FLUIDS_GAMMALAWINDEX);
    fluids_setattrib(fluid[n], P, FLUIDS_PRIMITIVE);
  }
}

void finish()
{
  for (int n=0; n<100; ++n) {
    fluids_del(fluid[n]);
  }
}

int timederiv(double *L)
{
  fish_state *S = fish_new();
  fish_setfluid(S, FLUIDS_NRHYD);
  fish_setriemannsolver(S, FLUIDS_RIEMANN_EXACT);
  fish_setreconstruction(S, FISH_PLM);
  fish_setplmtheta(S, 2.0);

  double Fiph[500];
  fish_intercellflux(S, fluid, Fiph, 100, 0);

  for (int n=2; n<100-2; ++n) {
    for (int q=0; q<5; ++q) {
      L[5*n + q] = -(Fiph[5*n + q] - Fiph[5*(n-1) + q]) / dx;
    }
  }

  for (int q=0; q<5; ++q) {
    L[5* 0 + q] = 0.0;
    L[5* 1 + q] = 0.0;
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
    fluids_getattrib(fluid[n], &U0[5*n], FLUIDS_CONSERVED);
  }

  timederiv(L);
  for (int n=0; n<100; ++n) {
    double U[5];
    fluids_getattrib(fluid[n], U, FLUIDS_CONSERVED);
    for (int q=0; q<5; ++q) {
      U[q] += L[5*n + q] * dt * 0.5;
    }
    fluids_setattrib(fluid[n], U, FLUIDS_CONSERVED);
  }

  timederiv(L);
  for (int n=0; n<100; ++n) {
    double U[5];
    fluids_getattrib(fluid[n], U, FLUIDS_CONSERVED);
    for (int q=0; q<5; ++q) {
      U[q] = U0[5*n + q] + L[5*n + q] * dt;
    }
    fluids_setattrib(fluid[n], U, FLUIDS_CONSERVED);
  }
  return 0;
}

int main()
{
  init();
  for (int n=0; n<200; ++n) {
    advance();
  }
  FILE *outf = fopen("euler.dat", "w");
  for (int n=0; n<100; ++n) {
    double P[5];
    fluids_getattrib(fluid[n], P, FLUIDS_PRIMITIVE);
    fprintf(outf, "%d %f %f %f\n", n, P[0], P[1], P[2]);
  }
  fclose(outf);
  finish();
  return 0;
}
