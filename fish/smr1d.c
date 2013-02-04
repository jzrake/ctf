

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "fish.h"

static void timederiv(fish_state *scheme, fish_block *block, double *L);
static void set_bc(fish_state *scheme, fish_block *block, double *U, int nq);

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

    timederiv(scheme, block, L);
    for (int n=0; n<TotalZones; ++n) {
      fluids_state_derive(fluid[n], U1, FLUIDS_CONSERVED);
      for (int q=0; q<5; ++q) {
        U1[q] = U0[5*n + q] + L[5*n + q] * dt;
      }
      fluids_state_fromcons(fluid[n], U1, FLUIDS_CACHE_DEFAULT);
    }
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
  set_bc(scheme, block, L, 5);
}


void set_bc(fish_state *scheme, fish_block *block, double *U, int nq)
{
  int N = fish_block_totalstates(block);
  int Ng = fish_block_getguard(block);
  int BC;

  fish_getparami(scheme, &BC, FISH_BOUNDARY_CONDITIONS);

  switch (BC) {
  case FISH_PERIODIC:
    for (int i=0; i<Ng; ++i) {
      for (int q=0; q<nq; ++q) {
	U[(    i    )*nq + q] = U[( N - 2*Ng + i)*nq + q];
	U[(N - i - 1)*nq + q] = U[(-1 + 2*Ng - i)*nq + q];
      }
    }
    break;
  case FISH_OUTFLOW:
    for (int i=0; i<Ng; ++i) {
      for (int q=0; q<nq; ++q) {
	U[(    i    )*nq + q] = U[(    Ng    )*nq + q];
	U[(N - i - 1)*nq + q] = U[(N - Ng - 1)*nq + q];
      }
    }
    break;
  default:
    break; // warning!
  }
}
