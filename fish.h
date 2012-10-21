
#ifndef FISH_HEADER_INCLUDED
#define FISH_HEADER_INCLUDED

enum {
  FISH_PCM, // piecewise constant reconstruction
  FISH_PLM, // piecewise linear reconstruction
  FISH_WENO5, // weno-5 reconstruction
  FISH_GODUNOV, // conservative finite volume Riemann-solver intercell fluxes
  FISH_SPECTRAL, // conservative finite differencing of characteristic fields

  FISH_SCHEME,
  FISH_RIEMANN_SOLVER,
  FISH_RECONSTRUCTION,
  FISH_PLM_THETA,

  FISH_ERROR_BADARG,
} ;

#include "fluids.h"

typedef struct fish_state fish_state;

fish_state *fish_new(void);
int fish_del(fish_state *S);
int fish_evolve(fish_state *S, double dt);
int fish_intercellflux(fish_state *S, fluids_state **fluid, double *F, int N,
                       int dim);
int fish_getparami(fish_state *S, int *param, long flag);
int fish_setparami(fish_state *S, int param, long flag);
int fish_getparamd(fish_state *S, double *param, long flag);
int fish_setparamd(fish_state *S, double param, long flag);

#ifdef FISH_PRIVATE_DEFS
struct fish_state {
  int scheme;
  int riemann_solver;
  int reconstruction;
  double plm_theta;
} ;

#endif // FISH_PRIVATE_DEFS
#endif // FISH_HEADER_INCLUDED
