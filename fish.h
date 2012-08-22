
#ifndef FISH_HEADER_INCLUDED
#define FISH_HEADER_INCLUDED

#define FISH_FLAGSALL          ((1<<30) - 1)

#include "cow.h"
#include "fluids.h"

struct fish_state;
typedef struct fish_state fish_state;

fish_state *fish_new(void);
int fish_del(fish_state *S);

int fish_setstate(fish_state *S, cow_domain *d);
int fish_evolve(fish_state *S, double dt);
int fish_intercellflux(fish_state *S, fluid_state **fluid, double *F, int N);


#ifdef FISH_PRIVATE_DEFS
struct fish_state {

} ;

/* http://en.wikipedia.org/wiki/Bitwise_operation#NOT */
#define BITWISENOT(x) (-(x) - 1)

#endif // FISH_PRIVATE_DEFS

#endif // FISH_HEADER_INCLUDED
