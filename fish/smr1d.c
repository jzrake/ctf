

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define FISH_PRIVATE_DEFS
#include "fish.h"

double fish_block_maxwavespeed(fish_block *block)
{
  int TotalZones = fish_block_totalstates(block);
  fluids_state **fluid = fish_block_getfluid(block);
  double a = 0.0;
  for (int n=0; n<TotalZones; ++n) {
    double A[5];
    fluids_state_derive(fluid[n], A, FLUIDS_EVAL0);
    for (int q=0; q<5; ++q) {
      a = fabs(A[q]) > a ? fabs(A[q]) : a;
    }
  }
  return a;
}

int fish_block_timederivative(fish_block *block, fish_state *scheme)
{
  int TotalZones = fish_block_totalstates(block);
  double *L = block->time_derivative;
  fluids_state **fluid = block->fluid;
  double dx;
  for (int m=0; m<5*TotalZones; ++m) L[m] = 0.0;
  fish_block_gridspacing(block, 0, &dx);
  fish_timederivative(scheme, fluid, 1, &TotalZones, &dx, L);
  return 0;
}

int fish_block_evolve(fish_block *block, double dt)
{
  int TotalZones = fish_block_totalstates(block);
  double *L = block->time_derivative;
  double U1[5];

  fluids_state **fluid = fish_block_getfluid(block);

  for (int n=0; n<TotalZones; ++n) {
    fluids_state_derive(fluid[n], U1, FLUIDS_CONSERVED);
    for (int q=0; q<5; ++q) {
      U1[q] += L[5*n + q] * dt;
    }
    fluids_state_fromcons(fluid[n], U1, FLUIDS_CACHE_DEFAULT);
  }

  return 0;
}

int fish_block_fillguard(fish_block *block)
{
  int Ng = fish_block_getguard(block);
  fish_block *BL, *BR;

  fish_block_getneighbor(block, 0, FISH_LEFT, &BL);
  fish_block_getneighbor(block, 0, FISH_RIGHT, &BR);

  fluids_state **fluid0 = fish_block_getfluid(block);
  fluids_state **fluidL = fish_block_getfluid(BL);
  fluids_state **fluidR = fish_block_getfluid(BR);

  int Nx0 = fish_block_getsize(block, 0);
  int NxL = fish_block_getsize(BL, 0);

  for (int n=0; n<Ng; ++n) {
    fluids_state_copy(fluid0[n           ], fluidL[NxL + n]);
    fluids_state_copy(fluid0[Nx0 + Ng + n], fluidR[Ng  + n]);
  }
  return 0;
}
