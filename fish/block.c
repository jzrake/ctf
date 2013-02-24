

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define FISH_PRIVATE_DEFS
#include "fish.h"

#define CHECK(c,m) do{if(!(c)){B->error=m;return -1;}B->error=NULL;}while(0)
#define CHECK2(c,m,e) do{if(!(c)){B->error=m;return e;}B->error=NULL;}while(0)

#define ITERATE_INTERIOR			\
  int Nd = B->rank;				\
  int Ng = B->guard;				\
  int Nx = B->size[0];				\
  int Ny = B->size[1];				\
  int Nz = B->size[2];				\
						\
  int i0 = Nd >= 1 ? Ng : 0;			\
  int j0 = Nd >= 2 ? Ng : 0;			\
  int k0 = Nd >= 3 ? Ng : 0;			\
  int i1 = i0 + Nx;				\
  int j1 = j0 + Ny;				\
  int k1 = k0 + Nz;				\
  int sx=1,sy=1,sz=1;				\
						\
  switch (Nd) {					\
  case 1:					\
    sx = 1;					\
    break;					\
  case 2:					\
    sx = (Ny + 2*Ng);				\
    sy = 1;					\
    break;					\
  case 3:					\
    sx = (Ny + 2*Ng) * (Nz + 2*Ng);		\
    sy = (Nz + 2*Ng);				\
    sz = 1;					\
    break;					\
  }						\
						\
  for (int i=i0; i<i1; ++i) 			\
    for (int j=j0; j<j1; ++j) 			\
      for (int k=k0; k<k1; ++k) 		\


fish_block *fish_block_new()
{
  fish_block *B = (fish_block*) malloc(sizeof(fish_block));
  fish_block block = {
    .allocated = 0,
    .rank = 1,
    .guard = 0,
    .bcL = { FISH_NONE, FISH_NONE, FISH_NONE },
    .bcR = { FISH_NONE, FISH_NONE, FISH_NONE },
    .size = { 1, 1, 1 },
    .x0 = { 0.0, 0.0, 0.0 },
    .x1 = { 1.0, 1.0, 1.0 },
    .parent = NULL,
    .boundaryL = { NULL, NULL, NULL },
    .boundaryR = { NULL, NULL, NULL },
    .children = { NULL, NULL, NULL, NULL,
                  NULL, NULL, NULL, NULL },
    .fluid = NULL,
    .descr = NULL,

    .primitive = NULL,
    .gravity = NULL,
    .passive = NULL,
    .magnetic = NULL,
    .location = NULL,

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

int fish_block_getboundaryflag(fish_block *B, int dim, int LR, int *flag)
{
  CHECK(dim < B->rank,
        "argument 'dim' must be smaller than the rank of the block");
  CHECK(LR == FISH_LEFT || LR == FISH_RIGHT,
        "argument 'LR' must be FISH_LEFT or FISH_RIGHT");
  switch (LR) {
  case FISH_LEFT : *flag = B->bcL[dim]; break;
  case FISH_RIGHT: *flag = B->bcR[dim]; break;
  }
  return 0;
}

int fish_block_setboundaryflag(fish_block *B, int dim, int LR, int flag)
{
  CHECK(dim < B->rank,
        "argument 'dim' must be smaller than the rank of the block");
  CHECK(LR == FISH_LEFT || LR == FISH_RIGHT,
        "argument 'LR' must be FISH_LEFT or FISH_RIGHT");
  switch (LR) {
  case FISH_LEFT : B->bcL[dim] = flag; break;
  case FISH_RIGHT: B->bcR[dim] = flag; break;
  }
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

int fish_block_totalstates(fish_block *B, int mode)
{
  CHECK(1, NULL);
  int ng = mode == FISH_INCLUDING_GUARD ? B->guard : 0;
  switch (B->rank) {
  case 1: return (B->size[0]+2*ng);
  case 2: return (B->size[0]+2*ng) * (B->size[1]+2*ng);
  case 3: return (B->size[0]+2*ng) * (B->size[1]+2*ng) * (B->size[2]+2*ng);
  default: return 0;
  }
}

int fish_block_allocate(fish_block *B)
{
#define MB fluids_state_mapbuffer
  CHECK(!B->allocated, NULL); // will return if allocated, leave no message
  CHECK(B->descr, "block needs a fluid descriptor");

  int nz = fish_block_totalstates(B, FISH_INCLUDING_GUARD);
  int np = fluids_descr_getncomp(B->descr, FLUIDS_PRIMITIVE);
  int ng = fluids_descr_getncomp(B->descr, FLUIDS_GRAVITY);
  int ns = fluids_descr_getncomp(B->descr, FLUIDS_PASSIVE);
  int nm = fluids_descr_getncomp(B->descr, FLUIDS_MAGNETIC);
  int nl = fluids_descr_getncomp(B->descr, FLUIDS_LOCATION);

  CHECK(!np || B->primitive, "primitive not mapped to a buffer");
  CHECK(!ng || B->gravity, "gravity not mapped to a buffer");
  CHECK(!ns || B->passive, "passive not mapped to a buffer");
  CHECK(!nm || B->magnetic, "magnetic not mapped to a buffer");
  CHECK(!nl || B->location, "location not mapped to a buffer");

  B->time_derivative = (double*) malloc(nz * np * sizeof(double));
  B->temp_conserved  = (double*) malloc(nz * np * sizeof(double));
  B->fluid = (fluids_state**) malloc(nz * sizeof(fluids_state*));

  for (int n=0; n<nz; ++n) {
    B->fluid[n] = fluids_state_new();
    fluids_state_setdescr(B->fluid[n], B->descr);
  }

  for (int n=0; n<nz; ++n) {
    if (np) MB(B->fluid[n], &B->primitive[np*n], FLUIDS_PRIMITIVE);
    if (ng) MB(B->fluid[n], &B->gravity[ng*n], FLUIDS_GRAVITY);
    if (ns) MB(B->fluid[n], &B->passive[ns*n], FLUIDS_PASSIVE);
    if (nm) MB(B->fluid[n], &B->magnetic[nm*n], FLUIDS_MAGNETIC);
    if (nl) MB(B->fluid[n], &B->location[nl*n], FLUIDS_LOCATION);
  }

  B->allocated = 1;
  return 0;
#undef MB
}

int fish_block_deallocate(fish_block *B)
{
  CHECK(B->allocated, NULL); // will return if not allocated, leave no message

  int ntot = fish_block_totalstates(B, FISH_INCLUDING_GUARD);
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
  //  CHECK(B->allocated, "block must already be allocated");
  CHECK(B->descr, "block needs a fluid descriptor");

  switch (flag) {
  case FLUIDS_PRIMITIVE: B->primitive = x; break;
  case FLUIDS_GRAVITY  : B->gravity = x; break;
  case FLUIDS_PASSIVE  : B->passive = x; break;
  case FLUIDS_MAGNETIC : B->magnetic = x; break;
  case FLUIDS_LOCATION : B->location = x; break;
  default: CHECK(0, "bad buffer class");
  }

  return 0;
}

fluids_state **fish_block_getfluid(fish_block *B)
{
  return B->fluid;
}

int fish_block_neighbor(fish_block *B, int dim, int LR, fish_block **B1)
{
  *B1 = NULL;
  CHECK(dim < B->rank,
        "argument 'dim' must be smaller than the rank of the block");
  CHECK(LR == FISH_LEFT || LR == FISH_RIGHT,
        "argument 'LR' must be FISH_LEFT or FISH_RIGHT");

  if (B->parent == NULL) {
    *B1 = NULL;
  }
  else {
    switch (LR) {
    case FISH_LEFT:
      if (B->pid == 1) {
        *B1 = B->parent->children[0];
      }
      else if (B->pid == 0) {
        fish_block *PL;
        fish_block_neighbor(B->parent, dim, FISH_LEFT, &PL);
        if (PL) {
          *B1 = PL->children[1];
        }
        else {
          *B1 = NULL;
        }
      }
      break;
    case FISH_RIGHT:
      if (B->pid == 0) {
        *B1 = B->parent->children[1];
      }
      else if (B->pid == 1) {
        fish_block *PR;
        fish_block_neighbor(B->parent, dim, FISH_RIGHT, &PR);
        if (PR) {
          *B1 = PR->children[0];
        }
        else {
          *B1 = NULL;
        }
      }
      break;
    }
  }
  return 0;
}

int fish_block_getchild(fish_block *B, int id, fish_block **B1)
{
  CHECK(id < 8, "argument 'id' must be smaller than 8");
  *B1 = B->children[id];
  return 0;
}

int fish_block_getboundaryblock(fish_block *B, int dim, int LR, fish_block **B1)
{
  CHECK(dim < B->rank,
        "argument 'dim' must be smaller than the rank of the block");
  CHECK(LR == FISH_LEFT || LR == FISH_RIGHT,
        "argument 'LR' must be FISH_LEFT or FISH_RIGHT");
  switch (LR) {
  case FISH_LEFT : *B1 = B->boundaryL[dim]; break;
  case FISH_RIGHT: *B1 = B->boundaryR[dim]; break;
  }
  return 0;
}

int fish_block_setboundaryblock(fish_block *B, int dim, int LR, fish_block *B1)
{
  CHECK(dim < B->rank,
        "argument 'dim' must be smaller than the rank of the block");
  CHECK(LR == FISH_LEFT || LR == FISH_RIGHT,
        "argument 'LR' must be FISH_LEFT or FISH_RIGHT");
  CHECK(fish_block_level(B) == fish_block_level(B1),
        "the boundary block must be at the same level");
  switch (LR) {
  case FISH_LEFT : B->boundaryL[dim] = B1; break;
  case FISH_RIGHT: B->boundaryR[dim] = B1; break;
  }
  return 0;
}

int fish_block_level(fish_block *B)
{
  CHECK(1, NULL);
  return B->parent ? fish_block_level(B->parent) + 1 : 0;
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

  double a = 0.0;

  ITERATE_INTERIOR {
    int m = i*sx+j*sy+k*sz;
    double A[5];
    fluids_state_derive(B->fluid[m], A, FLUIDS_EVAL0);
    for (int q=0; q<5; ++q) {
      a = fabs(A[q]) > a ? fabs(A[q]) : a;
    }
  }
  return a;
}

int fish_block_timederivative(fish_block *B, fish_state *scheme)
{
  CHECK(B->allocated, "block needs to be allocated");
  CHECK(B->descr, "block needs a fluid descriptor");

  int TotalZones = fish_block_totalstates(B, FISH_INCLUDING_GUARD);
  double *L = B->time_derivative;
  fluids_state **fluid = B->fluid;
  double dx[3];
  int Ntot[3];

  for (int n=0; n<B->rank; ++n) {
    dx[n] = fish_block_gridspacing(B, n);
    Ntot[n] = B->size[n] + 2 * B->guard;
  }

  for (int m=0; m<5*TotalZones; ++m) L[m] = 0.0;

  fish_timederivative(scheme, fluid, B->rank, Ntot, dx, L);

  return 0;
}

int fish_block_sourceterms(fish_block *B)
/* -----------------------------------------------------------------------------
 *
 * Add the source terms onto the existing time-derivative data
 *
 * -----------------------------------------------------------------------------
 */
{
  CHECK(B->allocated, "block needs to be allocated");
  CHECK(B->descr, "block needs a fluid descriptor");

  double *L = B->time_derivative;
  double S[5];
  //  double dx = fish_block_gridspacing(B, 0);

  ITERATE_INTERIOR {
    int m = i*sx+j*sy+k*sz;
    fluids_state_derive(B->fluid[m], S, FLUIDS_SOURCETERMS);
    for (int q=0; q<5; ++q) {
      L[5*m + q] += S[q];
    }
  }

  return 0;
}

int fish_block_evolve(fish_block *B, double *W, double dt)
/* -----------------------------------------------------------------------------
 *
 * Evolve the fluids states forward a time dt with a low-storage Runge-Kutta
 * scheme. Use a weighted average of U^n (the conserved quantities at time level
 * n) and U^{n+1}. U^{n+1} is computed in-place based on the stored time
 * derivative L. The weight factors W[0] and W[1] must add to one.
 *
 * -----------------------------------------------------------------------------
 */
{
  CHECK(B->allocated, "block needs to be allocated");
  CHECK(B->descr, "block needs a fluid descriptor");

  double *U = B->temp_conserved; // U^n
  double *L = B->time_derivative;
  double U1[5];

  ITERATE_INTERIOR {
    int m = i*sx+j*sy+k*sz;
    fluids_state_derive(B->fluid[m], U1, FLUIDS_CONSERVED);

    for (int q=0; q<5; ++q) {
      double u0 = U[5*m + q]; // beginning of the time-step
      double u1 = U1[q] + L[5*m + q] * dt;
      U1[q] = W[0]*u0 + W[1]*u1;
    }
    int err = fluids_state_fromcons(B->fluid[m], U1, FLUIDS_CACHE_DEFAULT);
    CHECK2(!err, "conserved to primitive conversion failed", err);
  }
  return 0;
}

int fish_block_fillconserved(fish_block *B)
{
  CHECK(B->allocated, "block needs to be allocated");
  CHECK(B->descr, "block needs a fluid descriptor");

  double *U = B->temp_conserved;
  fluids_state **fluid = B->fluid;

  ITERATE_INTERIOR {
    int m = i*sx+j*sy+k*sz;
    int err = fluids_state_derive(fluid[m], &U[5*m], FLUIDS_CONSERVED);
    CHECK(!err, "found unphysical state");
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
  fish_block *BL = B->boundaryL[0];
  fish_block *BR = B->boundaryR[0];
  double Pl[5], Pr[5], P[5];

  int has_allocated_parent = B->parent ? B->parent->allocated : 0;

  if (BL == NULL) fish_block_neighbor(B, 0, FISH_LEFT, &BL);
  if (BR == NULL) fish_block_neighbor(B, 0, FISH_RIGHT, &BR);

  if (BL != NULL) { // fill from sibling (or boundary block) to the left
    CHECK(BL->allocated, "block needs an allocated sibling");
    for (int ic=0; ic<Ng; ++ic) {
      fluids_state_copy(B->fluid[ic], BL->fluid[Nx+ic]);
    }
  }

  else if (B->bcL[0] != FISH_NONE) { // fill with physical BC
    switch (B->bcL[0]) {
    case FISH_OUTFLOW:

      for (int i=0; i<Ng; ++i) {
	fluids_state_copy(B->fluid[i], B->fluid[Ng]);
      }

      break;
    default:
      CHECK(0, "only outflow BC currently supported");
      break;
    }
  }

  else {
    CHECK(has_allocated_parent,
          "block needs an allocated parent, boundary block, or flag");
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

  if (BR != NULL) { // fill from sibling (or boundary block) to the right
    CHECK(BR->allocated, "block needs an allocated sibling");
    for (int ic=0; ic<Ng; ++ic) {
      fluids_state_copy(B->fluid[Nx+Ng+ic], BR->fluid[Ng+ic]);
    }
  }

  else if (B->bcR[0] != FISH_NONE) { // fill with physical BC
    switch (B->bcR[0]) {
    case FISH_OUTFLOW:

      for (int i=0; i<Ng; ++i) {
	fluids_state_copy(B->fluid[Nx+Ng+i], B->fluid[Nx+Ng-1]);
      }

      break;
    default:
      CHECK(0, "only outflow BC currently supported");
      break;
    }
  }

  else {
    CHECK(has_allocated_parent,
          "block needs an allocated parent, boundary block, or flag");
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

int fish_block_project(fish_block *B)
/* -----------------------------------------------------------------------------
 *
 * Project data from the finer child grid to the corresponding sub-volume of the
 * parent grid.
 *
 * -----------------------------------------------------------------------------
 */
{
  CHECK(B->allocated, "block must be allocated");
  CHECK(B->descr, "block needs a fluid descriptor");
  CHECK(B->parent, "block needs a parent to project onto");
  CHECK(B->parent->allocated, "block needs an allocated parent");

  int Ng = B->guard;
  int Nx = B->size[0];

  fish_block *B0 = B->parent;
  double U0[5], U1[5], Up[5];
  int i0 = Ng + B->pstart[0];

  for (int n=Ng; n<Nx+Ng; n+=2) {

    fluids_state_derive(B->fluid[n+0], U0, FLUIDS_CONSERVED);
    fluids_state_derive(B->fluid[n+1], U1, FLUIDS_CONSERVED);

    for (int q=0; q<5; ++q) {
      Up[q] = 0.5 * (U0[q] + U1[q]);
    }

    fluids_state_fromcons(B0->fluid[i0+(n-Ng)/2], Up, FLUIDS_CACHE_DEFAULT);
  }

  return 0;
}

int fish_block_allocated(fish_block *B)
{
  CHECK(1, NULL);
  return B->allocated;
}

  
