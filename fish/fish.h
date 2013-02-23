
#ifndef FISH_HEADER_INCLUDED
#define FISH_HEADER_INCLUDED

enum {
  FISH_PCM, // piecewise constant reconstruction
  FISH_PLM, // piecewise linear reconstruction
  FISH_WENO5, // weno-5 reconstruction
  FISH_GODUNOV, // conservative finite volume Riemann-solver intercell fluxes
  FISH_SPECTRAL, // conservative finite differencing of characteristic fields
  FISH_DIFFUSION, // apply Lax-Friedrichs diffusion (for when things get dicey)

  // -------------------
  // boundary conditions
  // -------------------
  FISH_NONE,
  FISH_PERIODIC,
  FISH_OUTFLOW,

  // -------------------
  // time update methods
  // -------------------
  FISH_SINGLE,
  FISH_MIDPOINT,
  FISH_TVD_RK2,
  FISH_SHUOSHER_RK3,

  // ---------------------------------------------------------------------------
  // smoothness indicators for WENO reconstruction
  // ---------------------------------------------------------------------------
  FISH_ISK_JIANGSHU96, // original smoothness indicator of Jiang & Shu (1996)
  FISH_ISK_BORGES08, // improved by Borges (2008) NOTE: might be 4th order
  FISH_ISK_SHENZHA10, // improved by Shen & Zha (2010)

  // ---------------------------------------------------------------------------
  // names of parameters for solver description
  // ---------------------------------------------------------------------------

  // ------------------
  // integer parameters
  // ------------------
  FISH_SOLVER_TYPE,
  FISH_RIEMANN_SOLVER,
  FISH_RECONSTRUCTION,
  FISH_SMOOTHNESS_INDICATOR,
  FISH_BOUNDARY_CONDITIONS,
  FISH_TIME_UPDATE,

  // -----------------
  // double parameters
  // -----------------
  FISH_PLM_THETA, // [1 -> 2 (most aggressive)]
  FISH_SHENZHA10_PARAM, // [0 -> ~100 (most aggressive)]

  FISH_LEFT,
  FISH_RIGHT,
  FISH_INTERIOR,
  FISH_INCLUDING_GUARD,

  // -------------------
  // error codes
  // -------------------
  FISH_ERROR = -42,
  FISH_ERROR_BADARG = -43,
} ;

#include "fluids.h"

typedef struct fish_state fish_state;
typedef struct fish_block fish_block;

fish_state *fish_state_new(void);
int fish_state_del(fish_state *S);
int fish_intercellflux(fish_state *S, fluids_state **fluid, double *F, int N, int dim);
int fish_timederivative(fish_state *S, fluids_state **fluid, int ndim, int *shape, double *dx, double *dUdt);
int fish_getparami(fish_state *S, int *param, long flag);
int fish_setparami(fish_state *S, int param, long flag);
int fish_getparamd(fish_state *S, double *param, long flag);
int fish_setparamd(fish_state *S, double param, long flag);

fish_block *fish_block_new(void);
int fish_block_del(fish_block *B);
int fish_block_getrank(fish_block *B);
int fish_block_setrank(fish_block *B, int ndim);
int fish_block_getsize(fish_block *B, int dim);
int fish_block_setsize(fish_block *B, int dim, int size);
int fish_block_getrange(fish_block *B, int dim, double *x0, double *x1);
int fish_block_setrange(fish_block *B, int dim, double x0, double x1);
int fish_block_getguard(fish_block *B);
int fish_block_setguard(fish_block *B, int guard);
int fish_block_getdescr(fish_block *B, fluids_descr **D);
int fish_block_setdescr(fish_block *B, fluids_descr *D);
int fish_block_getchild(fish_block *B, int id, fish_block **B1);
int fish_block_setchild(fish_block *B, int id, fish_block *B1);
int fish_block_getboundaryblock(fish_block *B, int dim, int LR, fish_block **B1);
int fish_block_setboundaryblock(fish_block *B, int dim, int LR, fish_block *B1);
int fish_block_getboundaryflag(fish_block *B, int dim, int LR, int *flag);
int fish_block_setboundaryflag(fish_block *B, int dim, int LR, int flag);
int fish_block_totalstates(fish_block *B, int mode);
int fish_block_allocate(fish_block *B);
int fish_block_deallocate(fish_block *B);
int fish_block_mapbuffer(fish_block *B, double *x, long flag);
int fish_block_timederivative(fish_block *B, fish_state *scheme);
int fish_block_sourceterms(fish_block *B);
int fish_block_evolve(fish_block *B, double *W, double dt);
int fish_block_fillconserved(fish_block *B);
int fish_block_fillguard(fish_block *B);
int fish_block_project(fish_block *B);
int fish_block_allocated(fish_block *B);
int fish_block_neighbor(fish_block *B, int dim, int LR, fish_block **B1);
int fish_block_level(fish_block *B);
int fish_block_solvepoisson(fish_block *B);
double fish_block_gridspacing(fish_block *B, int dim);
double fish_block_positionatindex(fish_block *B, int dim, int index);
double fish_block_maxwavespeed(fish_block *B);
char *fish_block_geterror(fish_block *B);
fluids_state **fish_block_getfluid(fish_block *B);

#ifdef LUA_BLOCK // Functions only called from Lua code
void fish_block_map();
#endif

#ifdef FISH_PRIVATE_DEFS

enum { PCM_C2L, PCM_C2R,
       PLM_C2L, PLM_C2R,
       WENO5_FD_C2R, WENO5_FD_C2L,
       WENO5_FV_C2R, WENO5_FV_C2L,
       WENO5_FV_C2A, WENO5_FV_A2C };
double _reconstruct(fish_state *S, double *v, int type);

struct fish_state {
  int solver_type;
  int riemann_solver;
  int reconstruction;
  int smoothness_indicator;
  int boundary_conditions;
  int time_update;
  double plm_theta;
  double shenzha10_param;
} ;

struct fish_block {
  int allocated; // whether the block holds data
  int rank;      // dimensionality: 1, 2, or 3
  int guard;     // number of guard zones
  int bcL[3];    // boundary conditions flag (left)
  int bcR[3];    // boundary conditions flag (right)
  int size[3];   // number of zones (interior)
  double x0[3];  // lower domain bounds
  double x1[3];  // upper domain bounds
  struct fish_block *parent;
  struct fish_block *children[8];
  struct fish_block *boundaryL[3];
  struct fish_block *boundaryR[3];
  struct fluids_state **fluid;
  struct fluids_descr *descr;

  double *primitive;
  double *gravity;
  double *passive;
  double *magnetic;
  double *location;

  double *temp_conserved;
  double *time_derivative;
  char *error;
  int pid; // id in parent block array
  int pstart[3]; // starting indices into parent grid
} ;

#endif // FISH_PRIVATE_DEFS
#endif // FISH_HEADER_INCLUDED
