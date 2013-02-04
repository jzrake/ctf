
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
  FISH_PERIODIC,
  FISH_OUTFLOW,

  // -------------------
  // time update methods
  // -------------------
  FISH_MIDPOINT,
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
int fish_block_getneighbor(fish_block *B, int dim, int LR, fish_block **B1);
int fish_block_setneighbor(fish_block *B, int dim, int LR, fish_block *B1);
int fish_block_totalstates(fish_block *B);
int fish_block_gridspacing(fish_block *B, int dim, double *dx);
int fish_block_allocate(fish_block *B);
int fish_block_deallocate(fish_block *B);
int fish_block_mapbuffer(fish_block *B, double *x, long flag);
int fish_block_advance(fish_block *B, fish_state *scheme, double dt);
double fish_block_maxwavespeed(fish_block *B);
char *fish_block_geterror(fish_block *B);
fluids_state **fish_block_getfluid(fish_block *B);


void fish_grav1d_init(fluids_descr *descr, int N);
void fish_grav1d_finalize();
void fish_grav1d_setscheme(fish_state *S);
void fish_grav1d_advance(double dt);
void fish_grav1d_getprim(double *prim, double *grav);
void fish_grav1d_setprim(double *prim);
void fish_grav1d_mapbuffer(double *x, long flag);
double fish_grav1d_maxwavespeed();


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
  int rank;
  int guard;
  int size[3];
  double x0[3];
  double x1[3];
  struct fish_block *neighborL[3];
  struct fish_block *neighborR[3];
  struct fluids_state **fluid;
  struct fluids_descr *descr;
  double *time_derivative;
  char *error;
} ;

#endif // FISH_PRIVATE_DEFS
#endif // FISH_HEADER_INCLUDED
