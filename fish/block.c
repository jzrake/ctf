

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

int fish_block_gridspacing(fish_block *B, int dim, double *dx)
{
  if (dim < B->rank) {
    *dx = (B->x1[dim] - B->x0[dim]) / B->size[dim];
    return 0;
  }
  else {
    B->error = "argument 'dim' must be smaller than the rank of the block";
    return FISH_ERROR;
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
