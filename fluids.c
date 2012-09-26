

#define FLUIDS_PRIVATE_DEFS
#define FLUIDS_INDEX_VARS
#include "fluids.h"
#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

int _getsetattrib(fluids_state *S, double *x, long flag, char op);


fluids_cache *fluids_cache_new(void)
{
  fluids_cache *C = (fluids_cache*) malloc(sizeof(fluids_cache));
  fluids_cache cache = {
    .flux = { NULL, NULL, NULL },
    .eigenvalues = { NULL, NULL, NULL },
    .leigenvectors = { NULL, NULL, NULL },
    .reigenvectors = { NULL, NULL, NULL },
    .jacobian = { NULL, NULL, NULL },
    .fourvelocity = NULL,
    .soundspeedsquared = 0.0,
    .temperature = 0.0,
    .specificenthalpy = 0.0,
    .specificinternal = 0.0,
    .state = NULL,
    .needsupdateflags = FLUIDS_FLAGSALL,
  } ;
  *C = cache;
  return C;
} ;

int fluids_cache_del(fluids_cache *C)
{
  if (C != NULL) {
    free(C);
  }
  return 0;
}

fluids_descr *fluids_descr_new(void)
{
  fluids_descr *D = (fluids_descr*) malloc(sizeof(fluids_descr));
  fluids_descr descr = {
    .fluid = FLUIDS_NRHYD,
    .eos = FLUIDS_EOS_GAMMALAW,
    .coordsystem = FLUIDS_COORD_CARTESIAN,
    .nprimitive = 0,
    .npassive = 0,
    .ngravity = 0,
    .nmagnetic = 0,
    .nlocation = 0,
    .gammalawindex = 1.4,
  } ;
  *D = descr;
  return D;
} ;

int fluids_descr_del(fluids_descr *D)
{
  if (D != NULL) {
    free(D);
  }
  return 0;
}

int fluids_descr_getfluid(fluids_descr *D, int *fluid)
{
  *fluid = D->fluid;
  return 0;
}
int fluids_descr_setfluid(fluids_descr *D, int fluid)
{
  D->fluid = fluid;
  switch (fluid) {
  case FLUIDS_SCADV:
    D->nprimitive = 1;
    D->npassive = 0;
    D->ngravity = 0;
    D->nmagnetic = 0;
    D->nlocation = 0;
    break;
  case FLUIDS_NRHYD:
    D->nprimitive = 5;
    D->npassive = 0;
    D->ngravity = 0;
    D->nmagnetic = 0;
    D->nlocation = 0;
    break;
  }
  return 0;
}
int fluids_descr_geteos(fluids_descr *D, int *eos)
{
  *eos = D->eos;
  return 0;
}
int fluids_descr_seteos(fluids_descr *D, int eos)
{
  D->eos = eos;
  return 0;
}
int fluids_descr_getcoordsystem(fluids_descr *D, int *coordsystem)
{
  *coordsystem = D->coordsystem;
  return 0;
}
int fluids_descr_setcoordsystem(fluids_descr *D, int coordsystem)
{
  D->coordsystem = coordsystem;
  return 0;
}
int fluids_descr_getgamma(fluids_descr *D, double *gam)
{
  *gam = D->gammalawindex;
  return 0;
}
int fluids_descr_setgamma(fluids_descr *D, double gam)
{
  D->gammalawindex = gam;
  return 0;
}

fluids_state *fluids_state_new(void)
{
  fluids_state *S = (fluids_state*) malloc(sizeof(fluids_state));
  fluids_state state = {
    .primitive = NULL,
    .gravity = NULL,
    .location = NULL,
    .passive = NULL,
    .cache = NULL,
    .descr = NULL,
  } ;
  *S = state;
  return S;
}

int fluids_state_del(fluids_state *S)
{
  if (S != NULL) {
    free(S->primitive);
    free(S->passive);
    free(S->gravity);
    free(S->magnetic);
    free(S->location);
    fluids_cache_del(S->cache);
    free(S);
  }
  return 0;
}

int fluids_state_setdescr(fluids_state *S, fluids_descr *D)
{
  int np = D->nprimitive;
  int ns = D->npassive;
  int ng = D->ngravity;
  int nm = D->nmagnetic;
  int nl = D->nlocation;
  S->primitive = (double*) realloc(S->primitive, np * sizeof(double));
  S->passive = (double*) realloc(S->passive, ns * sizeof(double));
  S->gravity = (double*) realloc(S->gravity, ng * sizeof(double));
  S->magnetic = (double*) realloc(S->magnetic, nm * sizeof(double));
  S->location = (double*) realloc(S->location, nl * sizeof(double));
  S->descr = D;
  return 0;
}

int fluids_state_getattr(fluids_state *S, double *x, long flag)
{
  _getsetattrib(S, x, flag, 'g');
  return 0;
}

int fluids_state_setattr(fluids_state *S, double *x, long flag)
{
  _getsetattrib(S, x, flag, 's');
  return 0;
}





int _getsetattrib(fluids_state *S, double *x, long flag, char op)
{
  double *a = NULL;
  int size = 0;
#define CASE(f,m,s)case FLUIDS_##f: a = m; size = s; break
  switch (flag) {
    CASE(PRIMITIVE, S->primitive, S->descr->nprimitive);
    CASE(PASSIVE, S->passive, S->descr->npassive);
    CASE(GRAVITY, S->gravity, S->descr->ngravity);
    CASE(MAGNETIC, S->magnetic, S->descr->nmagnetic);
    CASE(LOCATION, S->location, S->descr->nlocation);
  default:
    break;
  }
#undef CASE
  if (op == 'g') { // get
    if (a == NULL) {
      return FLUIDS_ERROR_BADREQUEST;
    }
    memcpy(x, a, size * sizeof(double));
  }
  else if (op == 's') { // set
    if (a == NULL) {
      return FLUIDS_ERROR_BADREQUEST;
    }
    memcpy(a, x, size * sizeof(double));
  }
  return 0;
}
