

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

int advance()
{
  fish_state *S = fish_new();
  double L[500], Fiph[500];
  fish_intercellflux(S, fluid, Fiph, 100);

  for (int n=1; n<100-1; ++n) {
    for (int q=0; q<5; ++q) {
      L[5*n + q] = -(Fiph[5*n + q] - Fiph[5*(n-1) + q]) / dx;
    }
  }

  for (int q=0; q<5; ++q) {
    L[5* 0 + q] = 0.0;
    L[5*99 + q] = 0.0;
  }

  for (int n=0; n<100; ++n) {
    double U[5];
    fluids_getattrib(fluid[n], U, FLUIDS_CONSERVED);
    for (int q=0; q<5; ++q) {
      U[q] += L[5*n + q] * dt;
    }
    fluids_setattrib(fluid[n], U, FLUIDS_CONSERVED);
  }

  fish_del(S);
  return 0;
}

int main()
{
  init();
  for (int n=0; n<100; ++n) {
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
