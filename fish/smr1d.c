

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define FISH_PRIVATE_DEFS
#include "fish.h"

static void timederiv(fish_state *scheme, fish_block *block, double *L);
static void set_bc(fish_block *block);

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

int fish_block_advance(fish_block *block, fish_state *scheme, double dt)
{
  int TotalZones = fish_block_totalstates(block);
  fluids_state **fluid = fish_block_getfluid(block);
  double *L  = (double*) malloc(TotalZones * 5 * sizeof(double));
  double *U0 = (double*) malloc(TotalZones * 5 * sizeof(double));
  double U1[5];
  int method;
  fish_getparami(scheme, &method, FISH_TIME_UPDATE);

  for (int m=0; m < 5 * TotalZones; ++m) { U0[m] = 0.0; }
  for (int n=0; n<TotalZones; ++n) {
    fluids_state_derive(fluid[n], &U0[5*n], FLUIDS_CONSERVED);
  }

  switch (method) {
  case FISH_MIDPOINT:
    timederiv(scheme, block, L);
    for (int n=0; n<TotalZones; ++n) {
      fluids_state_derive(fluid[n], U1, FLUIDS_CONSERVED);
      for (int q=0; q<5; ++q) {
	U1[q] = U0[5*n + q] + L[5*n + q] * dt * 0.5;
      }
      fluids_state_fromcons(fluid[n], U1, FLUIDS_CACHE_DEFAULT);
    }
    set_bc(block);

    timederiv(scheme, block, L);
    for (int n=0; n<TotalZones; ++n) {
      fluids_state_derive(fluid[n], U1, FLUIDS_CONSERVED);
      for (int q=0; q<5; ++q) {
        U1[q] = U0[5*n + q] + L[5*n + q] * dt;
      }
      fluids_state_fromcons(fluid[n], U1, FLUIDS_CACHE_DEFAULT);
    }
    set_bc(block);
  break;

  case FISH_SHUOSHER_RK3:
    /* ******************************* Step 1 ******************************* */
    timederiv(scheme, block, L);
    for (int n=0; n<TotalZones; ++n) {
      fluids_state_derive(fluid[n], U1, FLUIDS_CONSERVED);
      for (int q=0; q<5; ++q) {
	int m = 5*n + q;
	U1[q] = U0[m] + L[m] * dt;
      }
      fluids_state_fromcons(fluid[n], U1, FLUIDS_CACHE_DEFAULT);
    }
    set_bc(block);

    /* ******************************* Step 2 ******************************* */
    timederiv(scheme, block, L);
    for (int n=0; n<TotalZones; ++n) {
      fluids_state_derive(fluid[n], U1, FLUIDS_CONSERVED);
      for (int q=0; q<5; ++q) {
	int m = 5*n + q;
	U1[q] = 3./4*U0[m] + 1./4*U1[q] + 1./4 * dt * L[m];
      }
      fluids_state_fromcons(fluid[n], U1, FLUIDS_CACHE_DEFAULT);
    }
    set_bc(block);

    /* ******************************* Step 3 ******************************* */
    timederiv(scheme, block, L);
    for (int n=0; n<TotalZones; ++n) {
      fluids_state_derive(fluid[n], U1, FLUIDS_CONSERVED);
      for (int q=0; q<5; ++q) {
	int m = 5*n + q;
	U1[q] = 1./3*U0[m] + 2./3*U1[q] + 2./3 * dt * L[m];
      }
      fluids_state_fromcons(fluid[n], U1, FLUIDS_CACHE_DEFAULT);
    }
    set_bc(block);
    break;
  }

  free(L);
  free(U0);
  return 0;
}


void timederiv(fish_state *scheme, fish_block *block, double *L)
{
  int TotalZones = fish_block_totalstates(block);
  double dx;
  fluids_state **fluid = fish_block_getfluid(block);
  for (int m=0; m < 5 * TotalZones; ++m) {
    L[m] = 0.0;
  }
  fish_block_gridspacing(block, 0, &dx);
  fish_timederivative(scheme, fluid, 1, &TotalZones, &dx, L);
}


void set_bc(fish_block *block)
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
}
