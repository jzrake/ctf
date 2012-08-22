
#ifndef FISH_HEADER_INCLUDED
#define FISH_HEADER_INCLUDED

#define FISH_FLAGSALL          ((1<<30) - 1)
#define FISH_NONE              -42
#define FISH_PLM               -43
#define FISH_WENO5             -44

#include "cow.h"
#include "fluids.h"

struct fish_state;
typedef struct fish_state fish_state;

fish_state *fish_new(void);
int fish_del(fish_state *S);

int fish_evolve(fish_state *S, double dt);
int fish_intercellflux(fish_state *S, fluid_state **fluid, double *F, int N,
		       int dim);
int fish_setfluid(fish_state *S, int fluid);
int fish_setriemannsolver(fish_state *S, int riemannsolver);
int fish_setreconstruction(fish_state *S, int reconstruction);
int fish_setplmtheta(fish_state *S, double plmtheta);

#ifdef FISH_PRIVATE_DEFS
struct fish_state {
  int fluid;
  int riemannsolver;
  int reconstruction;
  double plmtheta;
} ;

/* http://en.wikipedia.org/wiki/Bitwise_operation#NOT */
#define BITWISENOT(x) (-(x) - 1)

#endif // FISH_PRIVATE_DEFS
#endif // FISH_HEADER_INCLUDED
