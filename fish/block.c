

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define FISH_PRIVATE_DEFS
#include "fish.h"

#define CHECK(c,m) do{if(!(c)){B->error=m;return -1;}B->error=NULL;}while(0)

fish_block *fish_block_new()
{
  fish_block *B = (fish_block*) malloc(sizeof(fish_block));
  fish_block block = {
    .allocated = 0,
    .rank = 1,
    .guard = 0,
    .size = { 1, 1, 1 },
    .x0 = { 0.0, 0.0, 0.0 },
    .x1 = { 1.0, 1.0, 1.0 },
    .neighborL = { NULL, NULL, NULL },
    .neighborR = { NULL, NULL, NULL },
    .children = { NULL, NULL, NULL, NULL,
		  NULL, NULL, NULL, NULL },
    .parent = NULL,
    .fluid = NULL,
    .descr = NULL,
    .error = NULL,
    .pid = 0,
    .pstart = { 0, 0, 0 },
  } ;
  *B = block;
  return B;
}

int fish_block_del(fish_block *B)
{
  CHECK(1, NULL);
  if (B->allocated) {
    fish_block_deallocate(B);
  }
  free(B);
  return 0;
}

char *fish_block_geterror(fish_block *B)
{
  return B->error;
}

int fish_block_getsize(fish_block *B, int dim)
{
  CHECK(dim < B->rank,
	"argument 'dim' must be smaller than the rank of the block");
  return B->size[dim];
}

int fish_block_setsize(fish_block *B, int dim, int size)
{
  CHECK(dim < B->rank,
	"argument 'dim' must be smaller than the rank of the block");
  B->size[dim] = size;
  return 0;
}

int fish_block_getrange(fish_block *B, int dim, double *x0, double *x1)
{
  CHECK(dim < B->rank,
	"argument 'dim' must be smaller than the rank of the block");
  *x0 = B->x0[dim];
  *x1 = B->x1[dim];
  return 0;
}

int fish_block_setrange(fish_block *B, int dim, double x0, double x1)
{
  CHECK(dim < B->rank,
	"argument 'dim' must be smaller than the rank of the block");
  B->x0[dim] = x0;
  B->x1[dim] = x1;
  return 0;
}

int fish_block_getrank(fish_block *B)
{
  return B->rank;
}

int fish_block_setrank(fish_block *B, int rank)
{
  CHECK(!B->allocated, "cannot set rank when allocated");
  CHECK(rank >= 1 && rank <= 3, "rank must be 1, 2, or 3");
  B->rank = rank;
  return 0;
}

int fish_block_getguard(fish_block *B)
{
  CHECK(1, NULL);
  return B->guard;
}

int fish_block_setguard(fish_block *B, int guard)
{
  CHECK(!B->allocated, "cannot set guard zones when allocated");
  B->guard = guard;
  return 0;
}

int fish_block_getdescr(fish_block *B, fluids_descr **D)
{
  CHECK(1, NULL);
  *D = B->descr;
  return 0;
}

int fish_block_setdescr(fish_block *B, fluids_descr *D)
{
  CHECK(!B->allocated, "cannot set fluid descriptor when allocated");
  B->descr = D;
  return 0;
}

int fish_block_totalstates(fish_block *B)
{
  CHECK(1, NULL);
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
  CHECK(!B->allocated, NULL); // will return if allocated, leave no message
  CHECK(B->descr, "block needs a fluid descriptor");

  int ntot = fish_block_totalstates(B);
  int nprm = fluids_descr_getncomp(B->descr, FLUIDS_PRIMITIVE);

  B->time_derivative = (double*) malloc(ntot * nprm * sizeof(double));
  B->temp_conserved  = (double*) malloc(ntot * nprm * sizeof(double));
  B->fluid = (fluids_state**) malloc(ntot * sizeof(fluids_state*));

  for (int n=0; n<ntot; ++n) {
    B->fluid[n] = fluids_state_new();
    fluids_state_setdescr(B->fluid[n], B->descr);
  }

  B->allocated = 1;
  return 0;
}

int fish_block_deallocate(fish_block *B)
{
  CHECK(B->allocated, NULL); // will return if not allocated, leave no message

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
  B->allocated = 0;
  return 0;
}

int fish_block_mapbuffer(fish_block *B, double *x, long flag)
{
  CHECK(B->allocated, "block must already be allocated");
  CHECK(B->descr, "block needs a fluid descriptor");

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
  CHECK(dim < B->rank,
	"argument 'dim' must be smaller than the rank of the block");
  CHECK(LR == FISH_LEFT || LR == FISH_RIGHT,
	"argument 'LR' must be FISH_LEFT or FISH_RIGHT");
  switch (LR) {
  case FISH_LEFT : *B1 = B->neighborL[dim]; return 0;
  case FISH_RIGHT: *B1 = B->neighborR[dim]; return 0;
  }
  return -1;
}

int fish_block_setneighbor(fish_block *B, int dim, int LR, fish_block *B1)
{
  CHECK(dim < B->rank,
	"argument 'dim' must be smaller than the rank of the block");
  CHECK(LR == FISH_LEFT || LR == FISH_RIGHT,
	"argument 'LR' must be FISH_LEFT or FISH_RIGHT");
  switch (LR) {
    case FISH_RIGHT: B->neighborR[dim] = B1; return 0;
    case FISH_LEFT : B->neighborL[dim] = B1; return 0;
  }
  return -1;
}

int fish_block_getchild(fish_block *B, int id, fish_block **B1)
{
  CHECK(id < 8, "argument 'id' must be smaller than 8");
  *B1 = B->children[id];
  return 0;
}

int fish_block_setchild(fish_block *B, int id, fish_block *B1)
// -----------------------------------------------------------------------------
// Establish the block `B1` as the child block of `B` at location `id`, which
// labels the (up to) three dimensional location within the parent block. id=0
// means x-left/y-left/z-left, id=1 means x-right,y-left,z-left and so on. This
// function infers the physical location of the child block, over-riding
// whichever was previously set by the user.
// -----------------------------------------------------------------------------
{
  CHECK(id < 8, "argument 'id' must be smaller than 8");
  B->children[id] = B1;
  B1->parent = B;
  B1->pid = id;
  for (int n=0; n<B->rank; ++n) {
    if (id & 1 << n) { // this block is on the right of dimension n
      B1->pstart[n] = B->size[n] / 2;
      B1->x0[n] = (B->x0[n] + B->x1[n]) * 0.5;
      B1->x1[n] =  B->x1[n];
    }
    else { // this block is on the left of dimension n
      B1->pstart[n] = 0;
      B1->x0[n] =  B->x0[n];
      B1->x1[n] = (B->x0[n] + B->x1[n]) * 0.5;
    }
  }
  return 0;
}

double fish_block_gridspacing(fish_block *B, int dim)
{
  CHECK(dim < B->rank,
	"argument 'dim' must be smaller than the rank of the block");
  return (B->x1[dim] - B->x0[dim]) / B->size[dim];
}

double fish_block_positionatindex(fish_block *B, int dim, int index)
// -----------------------------------------------------------------------------
// Return the physical coordinates of the center of zone `index` along dimension
// `dim`. The index includes padding, so that i=0 refers to ng zones to the left
// of the block boundary, where ng is the number of guard zones.
// -----------------------------------------------------------------------------
{
  CHECK(dim < B->rank,
	"argument 'dim' must be smaller than the rank of the block");
  double dx = fish_block_gridspacing(B, dim);
  return B->x0[dim] + dx * (index - B->guard + 0.5);
}

double fish_block_maxwavespeed(fish_block *B)
{
  CHECK(B->allocated, "block needs to be allocated");
  CHECK(B->descr, "block needs a fluid descriptor");

  int TotalZones = fish_block_totalstates(B);
  fluids_state **fluid = fish_block_getfluid(B);
  double a = 0.0;
  for (int n=0; n<TotalZones; ++n) {
    double A[5];
    fluids_state_derive(fluid[n], A, FLUIDS_EVAL0);
    for (int q=0; q<5; ++q) {
      a = fabs(A[q]) > a ? fabs(A[q]) : a;
    }
  }
  CHECK(1, NULL); // clear error message
  return a;
}

int fish_block_timederivative(fish_block *B, fish_state *scheme)
{
  CHECK(B->allocated, "block needs to be allocated");
  CHECK(B->descr, "block needs a fluid descriptor");

  int TotalZones = fish_block_totalstates(B);
  double *L = B->time_derivative;
  fluids_state **fluid = B->fluid;
  double dx = fish_block_gridspacing(B, 0);
  for (int m=0; m<5*TotalZones; ++m) L[m] = 0.0;
  fish_timederivative(scheme, fluid, 1, &TotalZones, &dx, L);
  return 0;
}

int fish_block_evolve(fish_block *B, double *W, double dt)
{
  CHECK(B->allocated, "block needs to be allocated");
  CHECK(B->descr, "block needs a fluid descriptor");

  int TotalZones = fish_block_totalstates(B);
  double *U = B->temp_conserved;
  double *L = B->time_derivative;
  double U1[5];
  fluids_state **fluid = fish_block_getfluid(B);

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

int fish_block_fillconserved(fish_block *B)
{
  CHECK(B->allocated, "block needs to be allocated");
  CHECK(B->descr, "block needs a fluid descriptor");

  int TotalZones = fish_block_totalstates(B);
  int Nq = fluids_descr_getncomp(B->descr, FLUIDS_PRIMITIVE);
  double *U = B->temp_conserved;
  fluids_state **fluid = B->fluid;
  for (int n=0; n<TotalZones; ++n) {
    fluids_state_derive(fluid[n], &U[Nq*n], FLUIDS_CONSERVED);
  }
  return 0;
}

int fish_block_fillguard(fish_block *B)
/* -----------------------------------------------------------------------------
 *
 * Fill the guard zones of the block `B`. For each (rank-1) dimensional edge, if
 * if that edge borders another block of the same refinement level, then fill
 * guard zones from that block's interior. Otherwise, fill guard zones from the
 * overlying portion of the parent grid's interior.
 *
 *
 *
 * |-------------------------------------------------------|
 * |     |# # #|+ + + + + + + +|# # #|                     |
 * |                     |# # #|+ + + + + + + +|# # #|     |      <- child grids
 * |-------------------------------------------------------|
 * | #   #   # | +   +   +   +   +   +   +   + | #   #   # |      <- parent grid
 * |-------------------------------------------------------|
 *
 *                             ^ sibling grids exchange over fine-fine boundary
 *
 *             ^ parent grid fills child grid over jumps in refinement
 *
 *
 * KEY: (diagram assumes 3 guard zones)
 *
 *  +: interior zone
 *  #: Guard zone
 *  |: Interior or exterior block edge
 *
 * -----------------------------------------------------------------------------
 */
{
  CHECK(B->allocated, "block must be allocated");
  CHECK(B->descr, "block needs a fluid descriptor");

  int Ng = B->guard;
  int Nx = B->size[0];

  fish_block *B0 = B->parent;
  fish_block *BL = B->neighborL[0];
  fish_block *BR = B->neighborR[0];
  double Pl[5], Pr[5], P[5];

  if (BL != NULL) { // fill from sibling to the left
    for (int ic=0; ic<Ng; ++ic) {
      fluids_state_copy(B->fluid[ic], BL->fluid[Nx+ic]);
    }
  }
  else { // fill from parent
    int i0 = B->pstart[0];

    for (int ic=0; ic<Ng; ++ic) { // ic labels which guard zone we are filling

      // decide on the two indices into the parent block which surround the
      // current zone

      int n0 = i0 + (ic + Ng - 1) / 2;
      int n1 = n0 + 1;

      fluids_state_getattr(B0->fluid[n0], Pl, FLUIDS_PRIMITIVE);
      fluids_state_getattr(B0->fluid[n1], Pr, FLUIDS_PRIMITIVE);

      for (int q=0; q<5; ++q) {

	if ((Ng + ic) % 2 == 0) P[q] = Pr[q] - 0.25 * (Pr[q] - Pl[q]);
	else                    P[q] = Pl[q] + 0.25 * (Pr[q] - Pl[q]);

      }
      fluids_state_setattr(B->fluid[ic], P, FLUIDS_PRIMITIVE);
    }
  }
  if (BR != NULL) { // fill from sibling to the right
    for (int ic=0; ic<Ng; ++ic) {
      fluids_state_copy(B->fluid[Nx+Ng+ic], BR->fluid[Ng+ic]);
    }
  }
  else { // fill from parent
    int i0 = B->pstart[0] + Nx / 2 + Ng;

    for (int ic=0; ic<Ng; ++ic) { // ic labels which guard zone we are filling

      int n0 = i0 - 1 + (ic + 1) / 2;
      int n1 = n0 + 1;

      fluids_state_getattr(B0->fluid[n0], Pl, FLUIDS_PRIMITIVE);
      fluids_state_getattr(B0->fluid[n1], Pr, FLUIDS_PRIMITIVE);

      for (int q=0; q<5; ++q) {

	if (ic % 2 == 0) P[q] = Pr[q] - 0.25 * (Pr[q] - Pl[q]);
	else             P[q] = Pl[q] + 0.25 * (Pr[q] - Pl[q]);

      }
      fluids_state_setattr(B->fluid[Nx+Ng+ic], P, FLUIDS_PRIMITIVE);
    }
  }
  return 0;
}
