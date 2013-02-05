

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define FISH_PRIVATE_DEFS
#include "fish.h"


fish_block *fish_block_new()
{
  fish_block *B = (fish_block*) malloc(sizeof(fish_block));
  fish_block block = {
    .rank = 1,
    .guard = 0,
    .size = { 1, 1, 1 },
    .x0 = { 0.0, 0.0, 0.0 },
    .x1 = { 1.0, 1.0, 1.0 },
    .fluid = NULL,
    .descr = NULL,
    .error = NULL,
  } ;
  *B = block;
  return B;
}

int fish_block_del(fish_block *B)
{
  fish_block_deallocate(B);
  free(B);
  return 0;
}

char *fish_block_geterror(fish_block *B)
{
  return B->error;
}

int fish_block_getsize(fish_block *B, int dim)
{
  if (dim < B->rank) {
    return B->size[dim];
  }
  else {
    B->error = "argument 'dim' must be smaller than the rank of the block";
    return FISH_ERROR;
  }
}

int fish_block_setsize(fish_block *B, int dim, int size)
{
  if (dim < B->rank) {
    B->size[dim] = size;
    return 0;
  }
  else {
    B->error = "argument 'dim' must be smaller than the rank of the block";
    return FISH_ERROR;
  }
}

int fish_block_getrange(fish_block *B, int dim, double *x0, double *x1)
{
  if (dim < B->rank) {
    *x0 = B->x0[dim];
    *x1 = B->x1[dim];
    return 0;
  }
  else {
    B->error = "argument 'dim' must be smaller than the rank of the block";
    return FISH_ERROR;
  }
}

int fish_block_setrange(fish_block *B, int dim, double x0, double x1)
{
  if (dim < B->rank) {
    B->x0[dim] = x0;
    B->x1[dim] = x1;
    return 0;
  }
  else {
    B->error = "argument 'dim' must be smaller than the rank of the block";
    return FISH_ERROR;
  }
}

int fish_block_getrank(fish_block *B)
{
  return B->rank;
}

int fish_block_setrank(fish_block *B, int rank)
{
  if (rank >= 1 && rank <= 3) {
    B->rank = rank;
    return 0;
  }
  else {
    B->error = "rank must be 1, 2, or 3";
    return FISH_ERROR;
  }
}

int fish_block_getguard(fish_block *B)
{
  return B->guard;
}

int fish_block_setguard(fish_block *B, int guard)
{
  B->guard = guard;
  return 0;
}

int fish_block_getdescr(fish_block *B, fluids_descr **D)
{
  *D = B->descr;
  return 0;
}

int fish_block_setdescr(fish_block *B, fluids_descr *D)
{
  B->descr = D;
  return 0;
}

int fish_block_totalstates(fish_block *B)
{
  int ng = B->guard;
  switch (B->rank) {
  case 1: return (B->size[0]+2*ng);
  case 2: return (B->size[0]+2*ng) * (B->size[1]+2*ng);
  case 3: return (B->size[0]+2*ng) * (B->size[1]+2*ng) * (B->size[2]+2*ng);
  default: return 0;
  }
}

int fish_block_allocate(fish_block *B)
{
  if (B->descr == NULL) {
    B->error = "block's fluid descriptor must be set before allocating";
    return FISH_ERROR;
  }
  if (B->fluid != NULL) {
    fish_block_deallocate(B);
  }
  int ntot = fish_block_totalstates(B);
  int nprm = fluids_descr_getncomp(B->descr, FLUIDS_PRIMITIVE);

  B->time_derivative = (double*) malloc(ntot * nprm * sizeof(double));
  B->temp_conserved  = (double*) malloc(ntot * nprm * sizeof(double));
  B->fluid = (fluids_state**) malloc(ntot * sizeof(fluids_state*));

  for (int n=0; n<ntot; ++n) {
    B->fluid[n] = fluids_state_new();
    fluids_state_setdescr(B->fluid[n], B->descr);
  }

  return 0;
}

int fish_block_deallocate(fish_block *B)
{
  int ntot = fish_block_totalstates(B);
  if (B->fluid) {
    for (int n=0; n<ntot; ++n) {
      fluids_state_del(B->fluid[n]);
    }
    free(B->fluid);
    free(B->time_derivative);
    free(B->temp_conserved);
    B->fluid = NULL;
    B->time_derivative = NULL;
    B->temp_conserved = NULL;
  }
  return 0;
}

int fish_block_mapbuffer(fish_block *B, double *x, long flag)
{
  int nz = fish_block_totalstates(B);
  int nq = fluids_descr_getncomp(B->descr, flag);
  for (int n=0; n<nz; ++n) {
    fluids_state_mapbuffer(B->fluid[n], &x[nq*n], flag);
  }
  return 0;
}

fluids_state **fish_block_getfluid(fish_block *B)
{
  return B->fluid;
}

int fish_block_getneighbor(fish_block *B, int dim, int LR, fish_block **B1)
{
  if (dim < B->rank) {
    switch (LR) {
    case FISH_LEFT : *B1 = B->neighborL[dim]; return 0;
    case FISH_RIGHT: *B1 = B->neighborR[dim]; return 0;
    default:
      B->error = "argument 'LR' must be FISH_LEFT or FISH_RIGHT";
      return FISH_ERROR;
    }
  }
  else {
    B->error = "argument 'dim' must be smaller than the rank of the block";
    return FISH_ERROR;
  }
}

int fish_block_setneighbor(fish_block *B, int dim, int LR, fish_block *B1)
{
  if (dim < B->rank) {
    switch (LR) {
    case FISH_RIGHT: B->neighborR[dim] = B1; return 0;
    case FISH_LEFT : B->neighborL[dim] = B1; return 0;
    default:
      B->error = "argument 'LR' must be FISH_LEFT or FISH_RIGHT";
      return FISH_ERROR;
    }
  }
  else {
    B->error = "argument 'dim' must be smaller than the rank of the block";
    return FISH_ERROR;
  }
}

double fish_block_gridspacing(fish_block *B, int dim)
{
  if (dim < B->rank) {
    return (B->x1[dim] - B->x0[dim]) / B->size[dim];
  }
  else {
    B->error = "argument 'dim' must be smaller than the rank of the block";
    return 0.0;
  }
}

double fish_block_positionatindex(fish_block *B, int dim, int index)
// -----------------------------------------------------------------------------
// Returns the physical coordinates of the center of zone `index` along
// dimension `dim`. The index includes padding, so that i=0 refers to ng zones
// to the left of the block boundary, where ng is the number of guard zones.
// -----------------------------------------------------------------------------
{
  if (dim < B->rank) {
    double dx = fish_block_gridspacing(B, dim);
    return B->x0[dim] + dx * (index - B->guard + 0.5);
  }
  else {
    B->error = "argument 'dim' must be smaller than the rank of the block";
    return 0.0;
  }
}

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
  double dx = fish_block_gridspacing(block, 0);
  for (int m=0; m<5*TotalZones; ++m) L[m] = 0.0;
  fish_timederivative(scheme, fluid, 1, &TotalZones, &dx, L);
  return 0;
}

int fish_block_evolve(fish_block *block, double *W, double dt)
{
  int TotalZones = fish_block_totalstates(block);
  double *U = block->temp_conserved;
  double *L = block->time_derivative;
  double U1[5];

  fluids_state **fluid = fish_block_getfluid(block);

  for (int n=0; n<TotalZones; ++n) {
    fluids_state_derive(fluid[n], U1, FLUIDS_CONSERVED);
    for (int q=0; q<5; ++q) {
      double u1 = U1[q];
      double u0 = U[5*n + q];
      double du = L[5*n + q] * dt;
      U1[q] = W[0]*u0 + W[1]*u1 + W[2]*du;
    }
    fluids_state_fromcons(fluid[n], U1, FLUIDS_CACHE_DEFAULT);
  }

  return 0;
}

int fish_block_fillconserved(fish_block *block)
{
  int TotalZones = fish_block_totalstates(block);
  int Nq = fluids_descr_getncomp(block->descr, FLUIDS_PRIMITIVE);
  double *U = block->temp_conserved;
  fluids_state **fluid = block->fluid;
  for (int n=0; n<TotalZones; ++n) {
    fluids_state_derive(fluid[n], &U[Nq*n], FLUIDS_CONSERVED);
  }
  return 0;
}

int fish_block_fillguard(fish_block *block)
{
  int Ng = fish_block_getguard(block);
  fish_block *BL=NULL, *BR=NULL;

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
