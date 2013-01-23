
#ifndef FISH_HEADER_INCLUDED
#define FISH_HEADER_INCLUDED

enum {
  FISH_PCM, // piecewise constant reconstruction
  FISH_PLM, // piecewise linear reconstruction
  FISH_WENO5, // weno-5 reconstruction
  FISH_GODUNOV, // conservative finite volume Riemann-solver intercell fluxes
  FISH_SPECTRAL, // conservative finite differencing of characteristic fields
  FISH_DIFFUSION, // apply Lax-Friedrichs diffusion (for when things get dicey)

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

  // -----------------
  // double parameters
  // -----------------
  FISH_PLM_THETA, // [1 -> 2 (most aggressive)]
  FISH_SHENZHA10_PARAM, // [0 -> ~100 (most aggressive)]

  FISH_ERROR_BADARG,
} ;

#include "fluids.h"

typedef struct fish_state fish_state;

fish_state *fish_new(void);
int fish_del(fish_state *S);
int fish_intercellflux(fish_state *S, fluids_state **fluid, double *F, int N, int dim);
int fish_timederivative(fish_state *S, fluids_state **fluid, int ndim, int *shape, double *dx, double *dUdt);
int fish_getparami(fish_state *S, int *param, long flag);
int fish_setparami(fish_state *S, int param, long flag);
int fish_getparamd(fish_state *S, double *param, long flag);
int fish_setparamd(fish_state *S, double param, long flag);

void fish_grav1d_init();
void fish_grav1d_finalize();
void fish_grav1d_advance();
void fish_grav1d_getprim(double *prim, double *grav);
void fish_grav1d_setprim(double *prim);
void fish_grav1d_mapbuffer(double *x, long flag);

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
  double plm_theta;
  double shenzha10_param;
} ;

#endif // FISH_PRIVATE_DEFS
#endif // FISH_HEADER_INCLUDED
