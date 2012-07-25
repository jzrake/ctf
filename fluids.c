
#define FLUIDS_PRIVATE_DEFS
#include "fluids.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static void _alloc_state(fluid_state *S, long modes);
static void _dealloc_state(fluid_state *S, long modes);



fluid_state *fluids_new()
{
  fluid_state *S = (fluid_state*) malloc(sizeof(fluid_state));
  fluid_state state = {
    .fluid = FLUIDS_SCALAR_ADVECTION,
    .eos = FLUIDS_EOS_GAMMALAW,
    .coordsystem = FLUIDS_COORD_CARTESIAN,
    .nwaves = 0,
    .npassive = 0,
    .location = NULL,
    .conserved = NULL,
    .primitive = NULL,
    .flux0 = NULL,
    .flux1 = NULL,
    .flux2 = NULL,
    .passive = NULL,
    .magnetic = NULL,
    .fourvelocity = NULL,
    .eigenvalues = NULL,
    .leigenvectors = NULL,
    .reigenvectors = NULL,
    .soundspeedsquared = 0.0,
    .temperature = 0.0,
    .specificenthalpy = 0.0,
    .specificinternal = 0.0,
  } ;
  *S = state;
  return S;
}

int fluids_del(fluid_state *S)
{
  _dealloc_state(S, FLUIDS_FLAGSALL);
  free(S);
  return 0;
}

int fluids_setfluid(fluid_state *S, int fluid)
{
  switch (fluid) {
  case FLUIDS_SCALAR_ADVECTION:
    S->fluid = fluid;
    S->nwaves = 1;
    break;
  case FLUIDS_SCALAR_BURGERS:
    S->fluid = fluid;
    S->nwaves = 1;
    break;
  case FLUIDS_SHALLOW_WATER:
    S->fluid = fluid;
    S->nwaves = 4;
    break;
  case FLUIDS_NRHYD:
    S->fluid = fluid;
    S->nwaves = 5;
    break;
  }
  _alloc_state(S, FLUIDS_CONSERVED | FLUIDS_PRIMITIVE);
  return 0;
}

int fluids_setnpassive(fluid_state *S, int n)
{
  if (n < 0) {
    return FLUIDS_ERROR_BADARG;
  }
  else {
    S->npassive = n;
    return 0;
  }
}



#define __GETSETV(a,s)							\
  int fluids_get##a(fluid_state *S, double *x)				\
  {                                                                     \
    if (S->a == NULL) {                                                 \
      return FLUIDS_ERROR_BADREQUEST;                                   \
    }                                                                   \
    else {                                                              \
      memcpy(x, S->a, (s) * sizeof(double));                            \
      return 0;                                                         \
    }                                                                   \
  }                                                                     \
  int fluids_set##a(fluid_state *S, double *x)				\
  {                                                                     \
    S->a = (double*) realloc(S->a, (s) * sizeof(double));               \
    memcpy(S->a, x, (s) * sizeof(double));                              \
    return 0;                                                           \
  }                                                                     \

#define __GETSETS(a)							\
  int fluids_get##a(fluid_state *S, double *x)				\
  {                                                                     \
    *x = S->a;								\
    return 0;								\
  }                                                                     \
  int fluids_set##a(fluid_state *S, double *x)				\
  {									\
    S->a = *x;								\
    return 0;                                                           \
  }                                                                     \


__GETSETV(location, 3)
__GETSETV(passive, S->npassive)
__GETSETV(conserved, S->nwaves)
__GETSETV(primitive, S->nwaves)
__GETSETV(flux0, S->nwaves)
__GETSETV(flux1, S->nwaves)
__GETSETV(flux2, S->nwaves)
__GETSETV(magnetic, 3)
__GETSETV(fourvelocity, 4)
__GETSETV(eigenvalues, S->nwaves)
__GETSETV(leigenvectors, S->nwaves * S->nwaves)
__GETSETV(reigenvectors, S->nwaves * S->nwaves)
__GETSETS(soundspeedsquared)
__GETSETS(temperature)
__GETSETS(specificenthalpy)
__GETSETS(specificinternal)





void _alloc_state(fluid_state *S, long modes)
{
#define A(a,s,m) if(modes&m)S->a=(double*)realloc(S->a,(s)*sizeof(double))
  A(location, 3, FLUIDS_LOCATION);
  A(passive, S->npassive, FLUIDS_PASSIVE);
  A(conserved, S->nwaves, FLUIDS_CONSERVED);
  A(primitive, S->nwaves, FLUIDS_PRIMITIVE);
  A(flux0, S->nwaves, FLUIDS_FLUX0);
  A(flux1, S->nwaves, FLUIDS_FLUX1);
  A(flux2, S->nwaves, FLUIDS_FLUX2);
  A(magnetic, 3, FLUIDS_MAGNETIC);
  A(fourvelocity, 4, FLUIDS_FOURVELOCITY);
  A(eigenvalues, 4, FLUIDS_EIGENVALUES);
  A(leigenvectors, 4, FLUIDS_EIGENVECTORS);
  A(reigenvectors, 4, FLUIDS_EIGENVECTORS);
#undef A
}

void _dealloc_state(fluid_state *S, long modes)
{
#define A(a,s,m) if(modes&m)S->a=(double*)realloc(S->a, 0)
  A(location, 3, FLUIDS_LOCATION);
  A(passive, S->npassive, FLUIDS_PASSIVE);
  A(conserved, S->nwaves, FLUIDS_CONSERVED);
  A(primitive, S->nwaves, FLUIDS_PRIMITIVE);
  A(flux0, S->nwaves, FLUIDS_FLUX0);
  A(flux1, S->nwaves, FLUIDS_FLUX1);
  A(flux2, S->nwaves, FLUIDS_FLUX2);
  A(magnetic, 3, FLUIDS_MAGNETIC);
  A(fourvelocity, 4, FLUIDS_FOURVELOCITY);
  A(eigenvalues, 4, FLUIDS_EIGENVALUES);
  A(leigenvectors, 4, FLUIDS_EIGENVECTORS);
  A(reigenvectors, 4, FLUIDS_EIGENVECTORS);
#undef A
}


