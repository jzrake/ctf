

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "fish.h"

static fish_state *scheme;
static fish_block *block;
static fluids_descr *descr;
static double *primitive;
static int TotalZones;
static int Nx = 128;
static int Ng = 2;
static void timederiv(double *L);
static void set_bc(double *U, int nq);

int fish_smr1d()
{
  TotalZones = Nx + 2*Ng;

  descr = fluids_descr_new();
  scheme = fish_state_new();
  block = fish_block_new();

  fluids_descr_setfluid(descr, FLUIDS_NRHYD);

  fish_setparami(scheme, FLUIDS_RIEMANN_HLLC, FISH_RIEMANN_SOLVER);
  fish_setparami(scheme, FISH_PLM, FISH_RECONSTRUCTION);
  fish_setparami(scheme, FISH_GODUNOV, FISH_SOLVER_TYPE);
  fish_setparami(scheme, FISH_PERIODIC, FISH_BOUNDARY_CONDITIONS);
  fish_setparami(scheme, FISH_SHUOSHER_RK3, FISH_TIME_UPDATE);
  fish_setparamd(scheme, 2.0, FISH_PLM_THETA);

  fish_block_setdescr(block, descr);
  fish_block_setrank(block, 1);
  fish_block_setsize(block, 0, Nx);
  fish_block_setrange(block, 0, 0.0, 1.0);
  fish_block_setguard(block, Ng);
  fish_block_allocate(block);

  int ntot = fish_block_totalstates(block);
  primitive = (double*) malloc(ntot * 5 * sizeof(double));
  fish_block_mapbuffer(block, primitive, FLUIDS_PRIMITIVE);

  double dx;
  fish_block_gridspacing(block, 0, &dx);

  for (int n=0; n<Nx + 2*Ng; ++n) {
    double x = (n - Ng) * dx;
    primitive[5*n + 0] = 1.0 + 0.1 * sin(4 * M_PI * x);
    primitive[5*n + 1] = 1.0;
    primitive[5*n + 2] = 0.0;
    primitive[5*n + 3] = 0.0;
    primitive[5*n + 4] = 0.0;
  }

  fish_block_advance(block, 0.01);

  fluids_descr_del(descr);
  fish_state_del(scheme);
  fish_block_del(block);
  free(primitive);
  return 0;
}


int fish_block_advance(fish_block *block, double dt)
{
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
    timederiv(L);
    for (int n=0; n<TotalZones; ++n) {
      fluids_state_derive(fluid[n], U1, FLUIDS_CONSERVED);
      for (int q=0; q<5; ++q) {
	U1[q] = U0[5*n + q] + L[5*n + q] * dt * 0.5;
      }
      fluids_state_fromcons(fluid[n], U1, FLUIDS_CACHE_DEFAULT);
    }

    timederiv(L);
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
    timederiv(L);
    for (int n=0; n<TotalZones; ++n) {
      fluids_state_derive(fluid[n], U1, FLUIDS_CONSERVED);
      for (int q=0; q<5; ++q) {
	int m = 5*n + q;
	U1[q] = U0[m] + L[m] * dt;
      }
      fluids_state_fromcons(fluid[n], U1, FLUIDS_CACHE_DEFAULT);
    }

    /* ******************************* Step 2 ******************************* */
    timederiv(L);
    for (int n=0; n<TotalZones; ++n) {
      fluids_state_derive(fluid[n], U1, FLUIDS_CONSERVED);
      for (int q=0; q<5; ++q) {
	int m = 5*n + q;
	U1[q] = 3./4*U0[m] + 1./4*U1[q] + 1./4 * dt * L[m];
      }
      fluids_state_fromcons(fluid[n], U1, FLUIDS_CACHE_DEFAULT);
    }

    /* ******************************* Step 3 ******************************* */
    timederiv(L);
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


void timederiv(double *L)
{
  double dx;
  fluids_state **fluid = fish_block_getfluid(block);
  for (int m=0; m < 5 * TotalZones; ++m) {
    L[m] = 0.0;
  }
  fish_block_gridspacing(block, 0, &dx);
  fish_timederivative(scheme, fluid, 1, &TotalZones, &dx, L);
  set_bc(L, 5);
}


void set_bc(double *U, int nq)
{
  int N = TotalZones;
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
