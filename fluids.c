

#define FLUIDS_PRIVATE_DEFS
#define FLUIDS_INDEX_VARS
#include "fluids.h"
#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


static int _getsetcacheattr(fluids_cache *C, double *x, long flag, char op);
static int _getsetstateattr(fluids_state *S, double *x, long flag, char op);
static int _alloc_cache(fluids_cache *C, int op, int np, long flags);

#define FLUID_PROTOTYPES(nm)                                            \
  static int nm##_c2p(fluids_state *S, double *U);                      \
  static int nm##_p2c(fluids_state *S);                                 \
  static int nm##_sources(fluids_state *S);                             \
  static int nm##_cs2(fluids_state *S, double *cs2);                    \
  static int nm##_flux(fluids_state *S, long modes);                    \
  static int nm##_eigenval(fluids_state *S, long modes);                \
  static int nm##_eigenvec(fluids_state *S, int dim, int L, int R);     \
  static int nm##_jacobian(fluids_state *S, int dim);                   \

FLUID_PROTOTYPES()
FLUID_PROTOTYPES(_nrhyd)
FLUID_PROTOTYPES(_gravs)
FLUID_PROTOTYPES(_gravp)
FLUID_PROTOTYPES(_srhyd)

fluids_cache *fluids_cache_new(void)
{
  fluids_cache *C = (fluids_cache*) malloc(sizeof(fluids_cache));
  fluids_cache cache = {
    .conserved = NULL,
    .sourceterms = NULL,
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
    .needsupdateflags = FLUIDS_FLAGSALL,
    .state = NULL,
  } ;
  *C = cache;
  return C;
} ;

int fluids_cache_del(fluids_cache *C)
{
  if (C != NULL) {
    _alloc_cache(C, DEALLOC, 0, FLUIDS_FLAGSALL);
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
    .cacheflags = 0,
    .gammalawindex = 1.4,
    .rhobar = 0.0,
    .cache = NULL,
  } ;
  *D = descr;
  return D;
} ;

int fluids_descr_del(fluids_descr *D)
{
  if (D != NULL) {
    fluids_cache_del(D->cache);
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
  D->nprimitive = 0;
  D->npassive = 0;
  D->ngravity = 0;
  D->nmagnetic = 0;
  D->nlocation = 0;
  D->cacheflags = FLUIDS_FLAGSALL;

  switch (fluid) {
  case FLUIDS_SCADV:
    D->nprimitive = 1;
    break;
  case FLUIDS_NRHYD:
    D->nprimitive = 5;
    break;
  case FLUIDS_GRAVS:
    D->nprimitive = 5;
    D->ngravity = 4;
    break;
  case FLUIDS_GRAVP:
    D->nprimitive = 5;
    D->ngravity = 4;
    break;
  case FLUIDS_SRHYD:
    D->nprimitive = 5;
    break;
  }
  D->cache = fluids_cache_new();
  _alloc_cache(D->cache, ALLOC, D->nprimitive, D->cacheflags);
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
int fluids_descr_getrhobar(fluids_descr *D, double *rhobar)
{
  *rhobar = D->rhobar;
  return 0;
}
int fluids_descr_setrhobar(fluids_descr *D, double rhobar)
{
  D->rhobar = rhobar;
  return 0;
}

int fluids_descr_getncomp(fluids_descr *D, long flag)
{
  switch (flag) {
  case FLUIDS_PRIMITIVE: return D->nprimitive;
  case FLUIDS_PASSIVE: return D->npassive;
  case FLUIDS_GRAVITY: return D->ngravity;
  case FLUIDS_MAGNETIC: return D->nmagnetic;
  case FLUIDS_LOCATION: return D->nlocation;
  default: return FLUIDS_ERROR_BADARG;
  }
}

fluids_state *fluids_state_new(void)
{
  fluids_state *S = (fluids_state*) malloc(sizeof(fluids_state));
  fluids_state state = {
    .primitive = NULL,
    .gravity = NULL,
    .location = NULL,
    .passive = NULL,
    .ownscache = 0,
    .ownsbufferflags = FLUIDS_FLAGSALL,
    .cache = NULL,
    .descr = NULL,
  } ;
  *S = state;
  return S;
}

int fluids_state_del(fluids_state *S)
{
  if (S != NULL) {
    if (S->ownscache) {
      fluids_cache_del(S->cache);
    }
    if (S->ownsbufferflags & FLUIDS_PRIMITIVE) {
      free(S->primitive);
    }
    if (S->ownsbufferflags & FLUIDS_PASSIVE) {
      free(S->passive);
    }
    if (S->ownsbufferflags & FLUIDS_GRAVITY) {
      free(S->gravity);
    }
    if (S->ownsbufferflags & FLUIDS_MAGNETIC) {
      free(S->magnetic);
    }
    if (S->ownsbufferflags & FLUIDS_LOCATION) {
      free(S->location);
    }
    free(S);
  }
  return 0;
}

int fluids_state_getdescr(fluids_state *S, fluids_descr **D)
{
  *D = S->descr;
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
  S->ownsbufferflags = FLUIDS_FLAGSALL;
  S->ownscache = 0;
  S->cache = D->cache;
  D->cache->state = S;
  return 0;
}

int fluids_state_cache(fluids_state *S, int operation)
{
  switch (operation) {
  case FLUIDS_CACHE_NOTOUCH:
    break;
  case FLUIDS_CACHE_CREATE: /* has no effect if the state already has its own
                               cache */
    if (!S->ownscache) {
      S->cache = fluids_cache_new();
      S->cache->state = S;
      S->ownscache = 1;
      _alloc_cache(S->cache, ALLOC, S->descr->nprimitive, S->descr->cacheflags);
    }
  case FLUIDS_CACHE_STEAL: /* has no effect if the state already has its own
                              cache or the descriptor's cache is already
                              pointing to it */
    if (!S->ownscache && S->descr->cache->state != S) {
      S->cache = S->descr->cache;
      S->cache->state = S;
      S->cache->needsupdateflags = FLUIDS_FLAGSALL;
    }
  case FLUIDS_CACHE_RESET: /* all fields go out-of-date */
    S->cache->needsupdateflags = FLUIDS_FLAGSALL;
    break;
  case FLUIDS_CACHE_ERASE:
    if (S->ownscache) {
      fluids_cache_del(S->cache);
      S->cache = S->descr->cache;
      S->ownscache = 0;
    }
  default:
    return FLUIDS_ERROR_BADARG;
  }
  return 0;
}

int fluids_state_mapbuffer(fluids_state *S, double *buffer, long flag)
{
  switch (flag) {
  case FLUIDS_PRIMITIVE:
    if (S->ownsbufferflags & FLUIDS_PRIMITIVE) {
      free(S->primitive);
    }
    S->primitive = buffer;
    break;
  case FLUIDS_PASSIVE:
    if (S->ownsbufferflags & FLUIDS_PASSIVE) {
      free(S->passive);
    }
    S->passive = buffer;
    break;
  case FLUIDS_GRAVITY:
    if (S->ownsbufferflags & FLUIDS_GRAVITY) {
      free(S->gravity);
    }
    S->gravity = buffer;
    break;
  case FLUIDS_MAGNETIC:
    if (S->ownsbufferflags & FLUIDS_MAGNETIC) {
      free(S->magnetic);
    }
    S->magnetic = buffer;
    break;
  case FLUIDS_LOCATION:
    if (S->ownsbufferflags & FLUIDS_LOCATION) {
      free(S->location);
    }
    S->location = buffer;
    break;
  }
  S->ownsbufferflags &= ~flag;
  return 0;
}

int fluids_state_getattr(fluids_state *S, double *x, long flag)
{
  return _getsetstateattr(S, x, flag, 'g');
}

int fluids_state_setattr(fluids_state *S, double *x, long flag)
{
  fluids_state_cache(S, FLUIDS_CACHE_RESET);
  return _getsetstateattr(S, x, flag, 's');
}

int fluids_state_fromcons(fluids_state *S, double *U, int cache)
{
  if (cache != FLUIDS_CACHE_NOTOUCH) {
    fluids_state_cache(S, FLUIDS_CACHE_RESET);
  }
  return _c2p(S, U);
}

int fluids_state_derive(fluids_state *S, double *x, long flags)
/*
 * If `x` is not NULL, then `flag` must represent only a single field, and that
 * field will be (deep) copied into the array `x`. If `x` is NULL then it will
 * not be used, `flag` may reference many fields which will all be updated,
 * provided the state has a cache. This function will NOT have the side-effect
 * of creating a cache if the state does not have one.
 *
 *
 * First `modes` is augmented with other modes which are dependencies. For
 * example, the sound speed is required for eigenvalues, so if `modes` includes
 * the latter and not the former, then the sound speed is tacked onto `modes`.
 *
 * Then `modes` is stripped of all modes which are already current, in other
 * words do not have their needsupdate bit enabled.
 *
 * modes:                    001101011
 * needsupdateflags:         000010001
 * modes becomes:            000000001 (modes &= C->needsupdateflags)
 *
 * After the update is finished, any bit in needsupdateflags which is enabled in
 * modes should be set to zero. In other words:
 *
 * needsupdateflags:         000010001
 * modes:                    000000001
 * needsupdateflags becomes: 000010000 (needsupdateflags &= ~modes)
 */
{
  if (!S->ownscache) {
    fluids_state_cache(S, FLUIDS_CACHE_STEAL);
  }
  if (!(S->ownsbufferflags & FLUIDS_FLAGSALL)) {
    fluids_state_cache(S, FLUIDS_CACHE_RESET);
  }

  fluids_cache *C = S->cache;
  long modes = flags;
  long err;

  if (modes & FLUIDS_FLAGSALL) {
    modes |= FLUIDS_CONSERVED; // conserved quantities are used for everything
  }
  if (modes & FLUIDS_EVALSALL) {
    modes |= FLUIDS_SOUNDSPEEDSQUARED; // cs is used for eigenvalues
  }
  modes &= C->needsupdateflags;

  if (modes & FLUIDS_CONSERVED) {
    err = _p2c(S); if (err) return err;
  }
  if (modes & FLUIDS_SOUNDSPEEDSQUARED) {
    err = _cs2(S, &C->soundspeedsquared); if (err) return err;
  }

  err = _flux(S, modes); if (err) return err;
  err = _eigenval(S, modes); if (err) return err;

  if (modes & (FLUIDS_LEVECS0 | FLUIDS_REVECS0)) {
    err = _eigenvec(S, 0, modes & FLUIDS_LEVECS0, modes & FLUIDS_REVECS0);
    if (err) return err;
  }
  if (modes & (FLUIDS_LEVECS1 | FLUIDS_REVECS1)) {
    if (err) return err;
    err = _eigenvec(S, 1, modes & FLUIDS_LEVECS1, modes & FLUIDS_REVECS1);
  }
  if (modes & (FLUIDS_LEVECS2 | FLUIDS_REVECS2)) {
    err = _eigenvec(S, 2, modes & FLUIDS_LEVECS2, modes & FLUIDS_REVECS2);
    if (err) return err;
  }
  if (modes & FLUIDS_JACOBIAN0) {
    err = _jacobian(S, 0);
    if (err) return err;
  }
  if (modes & FLUIDS_JACOBIAN1) {
    err = _jacobian(S, 1);
    if (err) return err;
  }
  if (modes & FLUIDS_JACOBIAN2) {
    err = _jacobian(S, 2);
    if (err) return err;
  }
  if (modes & FLUIDS_SOURCETERMS) {
    err = _sources(S);
    if (err) return err;
  }

  C->needsupdateflags &= ~modes;
  /* only makes sense when flags contains a single bit */
  if (x != NULL) {
    return _getsetcacheattr(S->cache, x, flags, 'g');
  }
  else {
    return 0;
  }
}

int fluids_state_getcached(fluids_state *S, double *x, long flag)
{
  return _getsetcacheattr(S->cache, x, flag, 'g');
}



int _getsetcacheattr(fluids_cache *C, double *x, long flag, char op)
{
  double *a = NULL;
  int size = 0;
  int np = C->state->descr->nprimitive;
#define CASE(f,m,s)case FLUIDS_##f: a = m; size = s; break
  switch (flag) {
    CASE(CONSERVED, C->conserved, np);
    CASE(SOURCETERMS, C->sourceterms, np);
    CASE(FOURVELOCITY, C->fourvelocity, 4);
    CASE(FLUX0, C->flux[0], np);
    CASE(FLUX1, C->flux[1], np);
    CASE(FLUX2, C->flux[2], np);
    CASE(EVAL0, C->eigenvalues[0], np);
    CASE(EVAL1, C->eigenvalues[1], np);
    CASE(EVAL2, C->eigenvalues[2], np);
    CASE(LEVECS0, C->leigenvectors[0], np*np);
    CASE(LEVECS1, C->leigenvectors[1], np*np);
    CASE(LEVECS2, C->leigenvectors[2], np*np);
    CASE(REVECS0, C->reigenvectors[0], np*np);
    CASE(REVECS1, C->reigenvectors[1], np*np);
    CASE(REVECS2, C->reigenvectors[2], np*np);
    CASE(JACOBIAN0, C->jacobian[0], np*np);
    CASE(JACOBIAN1, C->jacobian[1], np*np);
    CASE(JACOBIAN2, C->jacobian[2], np*np);
    CASE(SOUNDSPEEDSQUARED, &C->soundspeedsquared, 1);
    CASE(TEMPERATURE, &C->temperature, 1);
    CASE(SPECIFICENTHALPY, &C->specificenthalpy, 1);
    CASE(SPECIFICINTERNAL, &C->specificinternal, 1);
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

int _getsetstateattr(fluids_state *S, double *x, long flag, char op)
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

int _alloc_cache(fluids_cache *C, int op, int np, long flags)
{
#define A(a,s,m) do {                                           \
    if (flags & m) {                                            \
      if (op == ALLOC) {                                        \
        C->a = (double*) realloc(C->a, (s)*sizeof(double));     \
      }                                                         \
      else if (op == DEALLOC) {                                 \
        free(C->a);                                             \
        C->a = NULL;                                            \
      }                                                         \
    }                                                           \
  } while (0)
  A(conserved, np, FLUIDS_CONSERVED);
  A(sourceterms, np, FLUIDS_SOURCETERMS);
  A(fourvelocity, 4, FLUIDS_FOURVELOCITY);
  A(flux[0], np, FLUIDS_FLUX0);
  A(flux[1], np, FLUIDS_FLUX1);
  A(flux[2], np, FLUIDS_FLUX2);
  A(eigenvalues[0], np, FLUIDS_EVAL0);
  A(eigenvalues[1], np, FLUIDS_EVAL1);
  A(eigenvalues[2], np, FLUIDS_EVAL2);
  A(leigenvectors[0], np*np, FLUIDS_LEVECS0);
  A(leigenvectors[1], np*np, FLUIDS_LEVECS1);
  A(leigenvectors[2], np*np, FLUIDS_LEVECS2);
  A(reigenvectors[0], np*np, FLUIDS_REVECS0);
  A(reigenvectors[1], np*np, FLUIDS_REVECS1);
  A(reigenvectors[2], np*np, FLUIDS_REVECS2);
  A(jacobian[0], np*np, FLUIDS_JACOBIAN0);
  A(jacobian[1], np*np, FLUIDS_JACOBIAN1);
  A(jacobian[2], np*np, FLUIDS_JACOBIAN2);
#undef A
  return 0;
}





int _nrhyd_c2p(fluids_state *S, double *U)
{
  if (U[ddd] < 0.0) {
    return FLUIDS_ERROR_NEGATIVE_DENSITY_CONS;
  }
  if (U[tau] < 0.0) {
    return FLUIDS_ERROR_NEGATIVE_ENERGY;
  }

  double gm1 = S->descr->gammalawindex - 1.0;
  double *P = S->primitive;
  P[rho] =  U[ddd];
  P[pre] = (U[tau] - 0.5*(U[Sx]*U[Sx] + U[Sy]*U[Sy] + U[Sz]*U[Sz])/U[ddd])*gm1;
  P[vx]  =  U[Sx] / U[ddd];
  P[vy]  =  U[Sy] / U[ddd];
  P[vz]  =  U[Sz] / U[ddd];

  if (P[pre] < 0.0) {
    return FLUIDS_ERROR_NEGATIVE_PRESSURE;
  }
  if (P[rho] < 0.0) {
    return FLUIDS_ERROR_NEGATIVE_DENSITY_PRIM;
  }
  return 0;
}

int _nrhyd_p2c(fluids_state *S)
{
  double gm1 = S->descr->gammalawindex - 1.0;
  double *U = S->cache->conserved;
  double *P = S->primitive;

  if (P[pre] < 0.0) {
    return FLUIDS_ERROR_NEGATIVE_PRESSURE;
  }
  if (P[rho] < 0.0) {
    return FLUIDS_ERROR_NEGATIVE_DENSITY_PRIM;
  }

  U[ddd] = P[rho];
  U[tau] = P[rho] * 0.5*(P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz]) + P[pre]/gm1;
  U[Sx]  = P[rho] * P[vx];
  U[Sy]  = P[rho] * P[vy];
  U[Sz]  = P[rho] * P[vz];
  return 0;
}

int _nrhyd_sources(fluids_state *S)
{
  double *T = S->cache->sourceterms;
  T[ddd] = 0.0;
  T[tau] = 0.0;
  T[Sx] = 0.0;
  T[Sy] = 0.0;
  T[Sz] = 0.0;
  return 0;
}

int _nrhyd_cs2(fluids_state *S, double *cs2)
{
  double gm = S->descr->gammalawindex;
  *cs2 = gm * S->primitive[pre] / S->primitive[rho];
  return 0;
}

int _nrhyd_flux(fluids_state *S, long modes)
{
  fluids_cache *C = S->cache;
  double *P = S->primitive;
  double *U = C->conserved;
  if (modes & FLUIDS_FLUX0) {
    C->flux[0][ddd] = U[ddd] * P[vx];
    C->flux[0][tau] = (U[tau] + P[pre]) * P[vx];
    C->flux[0][Sx] = U[Sx] * P[vx] + P[pre];
    C->flux[0][Sy] = U[Sy] * P[vx];
    C->flux[0][Sz] = U[Sz] * P[vx];
  }
  if (modes & FLUIDS_FLUX1) {
    C->flux[1][ddd] = U[ddd] * P[vy];
    C->flux[1][tau] = (U[tau] + P[pre]) * P[vy];
    C->flux[1][Sx] = U[Sx] * P[vy];
    C->flux[1][Sy] = U[Sy] * P[vy] + P[pre];
    C->flux[1][Sz] = U[Sz] * P[vy];
  }
  if (modes & FLUIDS_FLUX2) {
    C->flux[2][ddd] = U[ddd] * P[vz];
    C->flux[2][tau] = (U[tau] + P[pre]) * P[vz];
    C->flux[2][Sx] = U[Sx] * P[vz];
    C->flux[2][Sy] = U[Sy] * P[vz];
    C->flux[2][Sz] = U[Sz] * P[vz] + P[pre];
  }
  return 0;
}

int _nrhyd_eigenval(fluids_state *S, long modes)
{
  fluids_cache *C = S->cache;
  double *P = S->primitive;
  double a = sqrt(C->soundspeedsquared);
  if (modes & FLUIDS_EVAL0) {
    C->eigenvalues[0][0] = P[vx] - a;
    C->eigenvalues[0][1] = P[vx];
    C->eigenvalues[0][2] = P[vx];
    C->eigenvalues[0][3] = P[vx];
    C->eigenvalues[0][4] = P[vx] + a;
  }
  if (modes & FLUIDS_EVAL1) {
    C->eigenvalues[1][0] = P[vy] - a;
    C->eigenvalues[1][1] = P[vy];
    C->eigenvalues[1][2] = P[vy];
    C->eigenvalues[1][3] = P[vy];
    C->eigenvalues[1][4] = P[vy] + a;
  }
  if (modes & FLUIDS_EVAL2) {
    C->eigenvalues[2][0] = P[vz] - a;
    C->eigenvalues[2][1] = P[vz];
    C->eigenvalues[2][2] = P[vz];
    C->eigenvalues[2][3] = P[vz];
    C->eigenvalues[2][4] = P[vz] + a;
  }
  return 0;
}

int _nrhyd_eigenvec(fluids_state *S, int dim, int Lf, int Rt)
{
  double gm = S->descr->gammalawindex;
  double *P_ = S->primitive;
  double *U = S->cache->conserved;
  double *L = S->cache->leigenvectors[dim];
  double *R = S->cache->reigenvectors[dim];

  const double T[3][5][5] =
  // Tx
    {{{1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 1, 0, 0},
      {0, 0, 0, 1, 0},
      {0, 0, 0, 0, 1}},
  // Ty
     {{1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 0, 1, 0},
      {0, 0,-1, 0, 0},
      {0, 0, 0, 0, 1}},
  // Tz
     {{1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 0, 0, 1},
      {0, 0, 0, 1, 0},
      {0, 0,-1, 0, 0}}};

  const double V[3][5][5] = // T^{-1}
  // Vx
    {{{1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 1, 0, 0},
      {0, 0, 0, 1, 0},
      {0, 0, 0, 0, 1}},
  // Vy
     {{1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 0,-1, 0},
      {0, 0, 1, 0, 0},
      {0, 0, 0, 0, 1}},
  // Vz
     {{1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 0, 0,-1},
      {0, 0, 0, 1, 0},
      {0, 0, 1, 0, 0}}};

  double P[5];
  matrix_vector_product(T[dim][0], P_, P, 5, 5);

  const double gm1 = gm - 1.0;
  const double u = P[vx];
  const double v = P[vy];
  const double w = P[vz];
  const double V2 = u*u + v*v + w*w;
  const double a = sqrt(gm * P[pre] / P[rho]);
  const double H = (U[tau] + P[pre]) / P[rho];

  // --------------------------------------------------------------------------
  // Toro Equation 3.82 (rows are permuted to deal with Mara's convention on
  // the conserved quantities)
  // --------------------------------------------------------------------------
  const double RR[5][5] =
    {{       1,       1,      1,      0,      0 },  // rho
     { H - u*a, H + u*a, 0.5*V2,      v,      w },  // nrg
     {     u-a,     u+a,      u,      0,      0 },  // px
     {       v,       v,      v,      1,      0 },  // py
     {       w,       w,      w,      0,      1 }}; // pz
  // --------------------------------------------------------------------------
  // Toro Equation 3.83 up to (gam - 1) / (2*a^2)
  // --------------------------------------------------------------------------
  const double LL[5][5] =
    {{    H + (a/gm1)*(u-a),   1, -(u+a/gm1),        -v,        -w },
     {    H - (a/gm1)*(u+a),   1, -(u-a/gm1),        -v,        -w },
     { -2*H + (4/gm1)*(a*a),  -2,        2*u,       2*v,       2*w },
     {         -2*v*a*a/gm1,   0,          0, 2*a*a/gm1,         0 },
     {         -2*w*a*a/gm1,   0,          0,         0, 2*a*a/gm1 }};
  // --------------------------------------------------------------------------
  //                    rho, nrg,         px,        py,        pz
  // --------------------------------------------------------------------------

  matrix_matrix_product(V[dim][0], RR[0], R, 5, 5, 5);
  matrix_matrix_product(LL[0], T[dim][0], L, 5, 5, 5);

  // Replace the term in eqn 3.83 : (gam - 1) / (2*a^2)
  // ---------------------------------------------------------------------------
  const double norm = gm1 / (2*a*a);
  for (int i=0; i<25; ++i) {
    L[i] *= norm;
  }
  return 0;
}

int _nrhyd_jacobian(fluids_state *S, int dim)
{
  double gm = S->descr->gammalawindex;
  double g1 = gm - 1.0;
  double g2 = gm - 2.0;

  // inputs are Mara convention: {D,E,px,py,pz}
  double *U = S->cache->conserved;
  double D = U[0];
  double E = U[1];
  double u = U[2]/D;
  double v = U[3]/D;
  double w = U[4]/D;

  double nx = (dim == 0);
  double ny = (dim == 1);
  double nz = (dim == 2);
  double vn = u*nx + v*ny + w*nz;
  double ek = 0.5*(u*u + v*v + w*w);
  double p0 = (E-D*ek)*g1; // pressure
  double a2 = gm*p0/D;     // sound speed
  double h0 = a2 / g1 + ek;

  // output is A := dF{D,px,py,pz,E}/d{D,px,py,pz,E}, Toro's convention
  /*
  double A[5][5] = { { 0, nx, ny, nz, 0 },
                     { g1*ek*nx - u*vn,
                       1*vn - g2*u*nx,
                       u*ny - g1*v*nx,
                       u*nz - g1*w*nx, g1*nx },
                     { g1*ek*ny - v*vn,
                       v*nx - g1*u*ny,
                       1*vn - g2*v*ny,
                       v*nz - g1*w*ny, g1*ny },
                     { g1*ek*nz - w*vn,
                       w*nx - g1*u*nz,
                       w*ny - g1*v*nz,
                       1*vn - g2*w*nz, g1*nz },
                     { (g1*ek-h0)*vn,
                       h0*nx - g1*u*vn,
                       h0*ny - g1*v*vn,
                       h0*nz - g1*w*vn, gm*vn } };
  */
  // output is A := dF{D,E,px,py,pz}/d{D,E,px,py,pz}, Mara's convention
  double A[5][5] = { { 0, 0, nx, ny, nz },
                     { (g1*ek-h0)*vn, gm*vn,
                       h0*nx - g1*u*vn,
                       h0*ny - g1*v*vn,
                       h0*nz - g1*w*vn },
                     { g1*ek*nx - u*vn, g1*nx,
                       1*vn - g2*u*nx,
                       u*ny - g1*v*nx,
                       u*nz - g1*w*nx },
                     { g1*ek*ny - v*vn, g1*ny,
                       v*nx - g1*u*ny,
                       1*vn - g2*v*ny,
                       v*nz - g1*w*ny },
                     { g1*ek*nz - w*vn, g1*nz,
                       w*nx - g1*u*nz,
                       w*ny - g1*v*nz,
                       1*vn - g2*w*nz } };

  memcpy(S->cache->jacobian[dim], A[0], 25*sizeof(double));
  return 0;
}

int _gravs_c2p(fluids_state *S, double *U)
{
  return _nrhyd_c2p(S, U);
}
int _gravs_p2c(fluids_state *S)
{
  return _nrhyd_p2c(S);
}
int _gravs_sources(fluids_state *S)
{
  double *P = S->primitive;
  double *G = S->gravity;
  double *T = S->cache->sourceterms;
  double fx = -G[gph+0]; // grad_phi
  double fy = -G[gph+1];
  double fz = -G[gph+2];
  T[rho] = 0.0;
  T[tau] = P[rho] * (fx*P[vx] + fy*P[vy] + fz*P[vz]);
  T[Sx]  = P[rho] * fx;
  T[Sy]  = P[rho] * fy;
  T[Sz]  = P[rho] * fz;
  return 0;
}
int _gravs_cs2(fluids_state *S, double *cs2)
{
  return _nrhyd_cs2(S, cs2);
}
int _gravs_flux(fluids_state *S, long modes)
{
  return _nrhyd_flux(S, modes);
}
int _gravs_eigenval(fluids_state *S, long modes)
{
  return _nrhyd_eigenval(S, modes);
}
int _gravs_eigenvec(fluids_state *S, int dim, int L, int R)
{
  return _nrhyd_eigenvec(S, dim, L, L);
}
int _gravs_jacobian(fluids_state *S, int dim)
{
  return _nrhyd_jacobian(S, dim);
}


int _gravp_c2p(fluids_state *S, double *U)
{
  return _nrhyd_c2p(S, U);
}
int _gravp_p2c(fluids_state *S)
{
  return _nrhyd_p2c(S);
}
int _gravp_sources(fluids_state *S)
{
  double *P = S->primitive;
  double *G = S->gravity;
  double *T = S->cache->sourceterms;
  double fx = -G[gph+0]; // grad_phi
  double fy = -G[gph+1];
  double fz = -G[gph+2];
  T[rho] = 0.0;
  T[tau] = P[rho] * (fx*P[vx] + fy*P[vy] + fz*P[vz]);
  T[Sx]  = 0.0;
  T[Sy]  = 0.0;
  T[Sz]  = 0.0;
  return 0;
}
int _gravp_cs2(fluids_state *S, double *cs2)
{
  return _nrhyd_cs2(S, cs2);
}
int _gravp_flux(fluids_state *S, long modes)
{
  fluids_cache *C = S->cache;
  double *G = S->gravity;
  double gph2 = G[gph+0] * G[gph+0] + G[gph+1] * G[gph+1] + G[gph+2] * G[gph+2];
  double rhobar = S->descr->rhobar;

  _nrhyd_flux(S, modes);

  if (modes & FLUIDS_FLUX0) {
    C->flux[0][Sx] += G[gph+0] * G[gph+0] + rhobar*G[phi] - 0.5 * gph2;
    C->flux[0][Sy] += G[gph+0] * G[gph+1] + rhobar*G[phi];
    C->flux[0][Sz] += G[gph+0] * G[gph+2] + rhobar*G[phi];
  }
  if (modes & FLUIDS_FLUX1) {
    C->flux[1][Sx] += G[gph+1] * G[gph+0] + rhobar*G[phi];
    C->flux[1][Sy] += G[gph+1] * G[gph+1] + rhobar*G[phi] - 0.5 * gph2;
    C->flux[1][Sz] += G[gph+1] * G[gph+2] + rhobar*G[phi];
  }
  if (modes & FLUIDS_FLUX2) {
    C->flux[2][Sx] += G[gph+2] * G[gph+0] + rhobar*G[phi];
    C->flux[2][Sy] += G[gph+2] * G[gph+1] + rhobar*G[phi];
    C->flux[2][Sz] += G[gph+2] * G[gph+2] + rhobar*G[phi] - 0.5 * gph2;
  }
  return 0;
}
int _gravp_eigenval(fluids_state *S, long modes)
{
  return _nrhyd_eigenval(S, modes);
}
int _gravp_eigenvec(fluids_state *S, int dim, int L, int R)
{
  return _nrhyd_eigenvec(S, dim, L, R);
}
int _gravp_jacobian(fluids_state *S, int dim)
{
  return _nrhyd_jacobian(S, dim);
}

int _srhyd_c2p(fluids_state *S, double *U)
{
  static const double ERROR_TOLR = 1e-8;
  static const double BIGW = 1e12;
  static const int NEWTON_MAX_ITER = 50;

  if (U[ddd] < 0.0) {
    return FLUIDS_ERROR_NEGATIVE_DENSITY_CONS;
  }
  if (U[tau] < 0.0) {
    return FLUIDS_ERROR_NEGATIVE_ENERGY;
  }

  double gm = S->descr->gammalawindex;
  double *P = S->primitive;
  double D = U[ddd];
  double Tau = U[tau];
  double S2 = U[Sx]*U[Sx] + U[Sy]*U[Sy] + U[Sz]*U[Sz];

  int soln_found = 0;
  int n_iter = 0;
  double f,g,W_soln=1.0,p=P[pre];

  while (!soln_found) {
    double v2 = S2 / pow(Tau + D + p, 2);
    double W2 = 1.0 / (1.0 - v2);
    double W = sqrt(W2);
    double e = (Tau + D*(1.0 - W) + p*(1.0 - W2)) / (D*W);
    double Rho = D / W;
    double h = 1.0 + e + p/Rho;
    double cs2 = gm * p / (Rho * h);

    f = Rho * e * (gm - 1.0) - p;
    g = v2*cs2 - 1.0;
    p -= f/g;

    if (fabs(f) < ERROR_TOLR) {
      W_soln = W;
      soln_found = 1;
    }
    if (n_iter++ == NEWTON_MAX_ITER) {
      return FLUIDS_ERROR_C2P_MAXITER;
    }
  }
  P[rho] = D/W_soln;
  P[pre] = p;
  P[vx] = U[Sx] / (Tau + D + p);
  P[vy] = U[Sy] / (Tau + D + p);
  P[vz] = U[Sz] / (Tau + D + p);

  if (P[pre] < 0.0) {
    return FLUIDS_ERROR_NEGATIVE_PRESSURE;
  }
  if (P[rho] < 0.0) {
    return FLUIDS_ERROR_NEGATIVE_DENSITY_PRIM;
  }
  if (W_soln != W_soln || W_soln > BIGW) {
    return FLUIDS_ERROR_SUPERLUMINAL;
  }
  return 0;
}

int _srhyd_p2c(fluids_state *S)
{
  double gm = S->descr->gammalawindex;
  double *P = S->primitive;
  double *U = S->cache->conserved;
  double V2 = P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz];
  double W2 = 1.0 / (1.0 - V2);
  double W = sqrt(W2);
  double e = P[pre] / (P[rho] * (gm - 1.0));
  double h = 1.0 + e + P[pre]/P[rho];

  if (P[pre] < 0.0) {
    return FLUIDS_ERROR_NEGATIVE_PRESSURE;
  }
  if (P[rho] < 0.0) {
    return FLUIDS_ERROR_NEGATIVE_DENSITY_PRIM;
  }
  if (V2 >= 1.0) {
    return FLUIDS_ERROR_SUPERLUMINAL;
  }

  U[ddd] = P[rho]*W;
  U[tau] = P[rho]*h*W2 - P[pre] - U[ddd];
  U[Sx] = P[rho]*h*W2*P[vx];
  U[Sy] = P[rho]*h*W2*P[vy];
  U[Sz] = P[rho]*h*W2*P[vz];

  if (U[tau] < 0.0) {
    return FLUIDS_ERROR_NEGATIVE_ENERGY;
  }
  return 0;
}

int _srhyd_sources(fluids_state *S)
{
  return _nrhyd_sources(S);
}

int _srhyd_cs2(fluids_state *S, double *cs2)
{
  double gm = S->descr->gammalawindex;
  double *P = S->primitive;
  double e = P[pre] / (P[rho] * (gm - 1.0));
  *cs2 = gm * P[pre] / (P[pre] + P[rho] + P[rho]*e);
  return 0;
}

int _srhyd_flux(fluids_state *S, long modes)
{
  return _nrhyd_flux(S, modes);
}

int _srhyd_eigenval(fluids_state *S, long modes)
{
  fluids_cache *C = S->cache;
  double *P = S->primitive;
  double cs2 = C->soundspeedsquared;
  double vx2 = P[vx]*P[vx];
  double vy2 = P[vy]*P[vy];
  double vz2 = P[vz]*P[vz];
  double v2 = vx2 + vy2 + vz2;
  if (modes & FLUIDS_EVAL0) {
    double x = sqrt(cs2*(1-v2)*(1-v2*cs2-vx2*(1-cs2)));
    C->eigenvalues[0][0] = (P[vx]*(1-cs2) - x)/(1-v2*cs2);
    C->eigenvalues[0][1] = P[vx];
    C->eigenvalues[0][2] = P[vx];
    C->eigenvalues[0][3] = P[vx];
    C->eigenvalues[0][4] = (P[vx]*(1-cs2) + x)/(1-v2*cs2);
  }
  if (modes & FLUIDS_EVAL1) {
    double x = sqrt(cs2*(1-v2)*(1-v2*cs2-vy2*(1-cs2)));
    C->eigenvalues[1][0] = (P[vy]*(1-cs2) - x)/(1-v2*cs2);
    C->eigenvalues[1][1] = P[vy];
    C->eigenvalues[1][2] = P[vy];
    C->eigenvalues[1][3] = P[vy];
    C->eigenvalues[1][4] = (P[vy]*(1-cs2) + x)/(1-v2*cs2);
  }
  if (modes & FLUIDS_EVAL2) {
    double x = sqrt(cs2*(1-v2)*(1-v2*cs2-vz2*(1-cs2)));
    C->eigenvalues[2][0] = (P[vz]*(1-cs2) - x)/(1-v2*cs2);
    C->eigenvalues[2][1] = P[vz];
    C->eigenvalues[2][2] = P[vz];
    C->eigenvalues[2][3] = P[vz];
    C->eigenvalues[2][4] = (P[vz]*(1-cs2) + x)/(1-v2*cs2);
  }
  return 0;
}

int _srhyd_eigenvec(fluids_state *S, int dim, int Lf, int Rt)
/* -----------------------------------------------------------------------------
 *
 * Authors: Jonathan Zrake, Bez Laderman: NYU CCPP
 *
 * Date: May 7th, 2012
 *
 * This piece of code implements the left and right eigenvectors of the ideal
 * special relativistic hydrodynamics equations. The formulas are an exact
 * translation of those given in the literature:
 *
 * R. Donat, J.A. Font, J.M. Ibanez, & A. Marquina
 * JCP, 1998, 146, 58
 *
 * http:*www.sciencedirect.com/science/article/pii/S0021999198959551
 *
 *
 * Having these eigenvectors in a hydrodynamics code is good. They can be used
 * for any scheme which requires flux splitting with characteristic
 * decomposition, such as high order ENO or WENO schemes.
 *
 * -----------------------------------------------------------------------------
 */
{
  double gm = S->descr->gammalawindex;
  double *P_ = S->primitive;
  double *L = S->cache->leigenvectors[dim];
  double *R = S->cache->reigenvectors[dim];

  const double T[3][5][5] =
  // Tx
    {{{1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 1, 0, 0},
      {0, 0, 0, 1, 0},
      {0, 0, 0, 0, 1}},
  // Ty
     {{1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 0, 1, 0},
      {0, 0,-1, 0, 0},
      {0, 0, 0, 0, 1}},
  // Tz
     {{1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 0, 0, 1},
      {0, 0, 0, 1, 0},
      {0, 0,-1, 0, 0}}};

  const double V[3][5][5] = // T^{-1}
  // Vx
    {{{1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 1, 0, 0},
      {0, 0, 0, 1, 0},
      {0, 0, 0, 0, 1}},
  // Vy
     {{1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 0,-1, 0},
      {0, 0, 1, 0, 0},
      {0, 0, 0, 0, 1}},
  // Vz
     {{1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 0, 0,-1},
      {0, 0, 0, 1, 0},
      {0, 0, 1, 0, 0}}};

  double P[5];
  matrix_vector_product(T[dim][0], P_, P, 5, 5);

  const double D = P[ddd]; // rest mass density
  const double p = P[pre]; // pressure
  const double u = P[vx];  // vx (three velocity)
  const double v = P[vy];  // vy
  const double w = P[vz];  // vz

  const double sie = (p/D) / (gm - 1); // specific internal energy
  const double h = 1 + sie + p/D;                  // specific enthalpy
  const double cs2 = gm * p / (D*h);   // sound speed squared
  const double V2 = u*u + v*v + w*w;
  const double W = 1.0 / sqrt(1 - V2);             // Lorentz factor
  const double W2 = W*W;
  const double K = h;                              // for gamma-law only, K = h
  const double hW = h*W;

  // equations (14) and (15)
  const double lp = (u*(1-cs2) + sqrt(cs2*(1-V2)*(1-V2*cs2-u*u*(1-cs2))))/(1-V2*cs2);
  const double lm = (u*(1-cs2) - sqrt(cs2*(1-V2)*(1-V2*cs2-u*u*(1-cs2))))/(1-V2*cs2);

  const double Ap = (1 - u*u) / (1 - u*lp);
  const double Am = (1 - u*u) / (1 - u*lm);

  // NOTES
  // ---------------------------------------------------------------------------
  // (1) Donat describes the columns if the right eigenvector matrix
  // horizontally, which is how they are written below. So we take the transpose
  // at the end of the day.
  //
  // (2) Donat's notation is provided next to the formulas to make clear the
  // permutation into to Mara's convention (vector components last).
  // ---------------------------------------------------------------------------


  // ---------------------------------------------------------------------------
  // Right eigenvectors (transpose of), equations (17) through (20)
  // ---------------------------------------------------------------------------
  const double RT[5][5] =
    {{1, hW*Am - 1, hW*Am*lm, hW*v, hW*w},                              // R_{-}
     {1, hW*Ap - 1, hW*Ap*lp, hW*v, hW*w},                              // R_{+}
     {K/hW, 1-K/hW, u, v, w},                                           // R_{1}
     {W*v, 2*h*W2*v - W*v, 2*h*W2*u*v, h*(1+2*W2*v*v), 2*h*W2*v*w},     // R_{2}
     {W*w, 2*h*W2*w - W*w, 2*h*W2*u*w, 2*h*W2*v*w, h*(1+2*W2*w*w)}};    // R_{3}

  double RR[5][5]; // un-transpose them
  for (int n=0; n<5; ++n) {
    for (int m=0; m<5; ++m) {
      RR[n][m] = RT[m][n];
    }
  }

  // ---------------------------------------------------------------------------
  // Left eigenvectors
  // ---------------------------------------------------------------------------
  const double Delta = h*h*h*W*(K-1)*(1-u*u)*(Ap*lp - Am*lm); // equation (21)
  const double a = W / (K-1);
  const double b = 1 / (h*(1 - u*u));
  const double c = 1 / (h*(1 - u*u));
  const double d = -h*h / Delta;
  const double e = +h*h / Delta;

  const double LL[5][5] =
    {{e*(hW*Ap*(u-lp) - u - W2*(V2 - u*u)*(2*K - 1)*(u - Ap*lp) + K*Ap*lp),
      e*(-u - W2*(V2 - u*u)*(2*K - 1)*(u - Ap*lp) + K*Ap*lp),
      e*(1 + W2*(V2 - u*u)*(2*K - 1)*(1 - Ap) - K*Ap),
      e*(W2*v*(2*K - 1)*Ap*(u - lp)),
      e*(W2*w*(2*K - 1)*Ap*(u - lp))},            // L_{-} (negative eigenvalue)
     {d*(hW*Am*(u-lm) - u - W2*(V2 - u*u)*(2*K - 1)*(u - Am*lm) + K*Am*lm),
      d*(-u - W2*(V2 - u*u)*(2*K - 1)*(u - Am*lm) + K*Am*lm),
      d*(1 + W2*(V2 - u*u)*(2*K - 1)*(1 - Am) - K*Am),
      d*(W2*v*(2*K - 1)*Am*(u - lm)),
      d*(W2*w*(2*K - 1)*Am*(u - lm))},            // L_{+} (positive eigenvalue)
     {a*(h-W), -a*W, a*W*u, a*W*v, a*W*w},        // L_{1}
     {-b*v, -b*v, b*u*v, b*(1-u*u), 0},           // L_{2}
     {-c*w, -c*w, c*u*w, 0, c*(1-u*u)}};          // L_{3}

  matrix_matrix_product(V[dim][0], RR[0], R, 5, 5, 5);
  matrix_matrix_product(LL[0], T[dim][0], L, 5, 5, 5);
  return 0;
}

int _srhyd_jacobian(fluids_state *S, int dim)
{
  return FLUIDS_ERROR_NOT_IMPLEMENTED;
}

int _c2p(fluids_state *S, double *U)
{
  switch (S->descr->fluid) {
  case FLUIDS_NRHYD: return _nrhyd_c2p(S, U);
  case FLUIDS_GRAVS: return _gravs_c2p(S, U);
  case FLUIDS_GRAVP: return _gravp_c2p(S, U);
  case FLUIDS_SRHYD: return _srhyd_c2p(S, U);
  default: return FLUIDS_ERROR_BADARG;
  }
}
int _p2c(fluids_state *S)
{
  switch (S->descr->fluid) {
  case FLUIDS_NRHYD: return _nrhyd_p2c(S);
  case FLUIDS_GRAVS: return _gravs_p2c(S);
  case FLUIDS_GRAVP: return _gravp_p2c(S);
  case FLUIDS_SRHYD: return _srhyd_p2c(S);
  default: return FLUIDS_ERROR_BADARG;
  }
}
int _sources(fluids_state *S)
{
  switch (S->descr->fluid) {
  case FLUIDS_NRHYD: return _nrhyd_sources(S);
  case FLUIDS_GRAVS: return _gravs_sources(S);
  case FLUIDS_GRAVP: return _gravp_sources(S);
  case FLUIDS_SRHYD: return _srhyd_sources(S);
  default: return FLUIDS_ERROR_BADARG;
  }
}
int _cs2(fluids_state *S, double *cs2)
{
  switch (S->descr->fluid) {
  case FLUIDS_NRHYD: return _nrhyd_cs2(S, cs2);
  case FLUIDS_GRAVS: return _gravs_cs2(S, cs2);
  case FLUIDS_GRAVP: return _gravp_cs2(S, cs2);
  case FLUIDS_SRHYD: return _srhyd_cs2(S, cs2);
  default: return FLUIDS_ERROR_BADARG;
  }
}
int _flux(fluids_state *S, long modes)
{
  switch (S->descr->fluid) {
  case FLUIDS_NRHYD: return _nrhyd_flux(S, modes);
  case FLUIDS_GRAVS: return _gravs_flux(S, modes);
  case FLUIDS_GRAVP: return _gravp_flux(S, modes);
  case FLUIDS_SRHYD: return _srhyd_flux(S, modes);
  default: return FLUIDS_ERROR_BADARG;
  }
}
int _eigenval(fluids_state *S, long modes)
{
  switch (S->descr->fluid) {
  case FLUIDS_NRHYD: return _nrhyd_eigenval(S, modes);
  case FLUIDS_GRAVS: return _gravs_eigenval(S, modes);
  case FLUIDS_GRAVP: return _gravp_eigenval(S, modes);
  case FLUIDS_SRHYD: return _srhyd_eigenval(S, modes);
  default: return FLUIDS_ERROR_BADARG;
  }
}
int _eigenvec(fluids_state *S, int dim, int L, int R)
{
  switch (S->descr->fluid) {
  case FLUIDS_NRHYD: return _nrhyd_eigenvec(S, dim, L, R);
  case FLUIDS_GRAVS: return _gravs_eigenvec(S, dim, L, R);
  case FLUIDS_GRAVP: return _gravp_eigenvec(S, dim, L, R);
  case FLUIDS_SRHYD: return _srhyd_eigenvec(S, dim, L, R);
  default: return FLUIDS_ERROR_BADARG;
  }
}
int _jacobian(fluids_state *S, int dim)
{
  switch (S->descr->fluid) {
  case FLUIDS_NRHYD: return _nrhyd_jacobian(S, dim);
  case FLUIDS_GRAVS: return _gravs_jacobian(S, dim);
  case FLUIDS_GRAVP: return _gravp_jacobian(S, dim);
  case FLUIDS_SRHYD: return _srhyd_jacobian(S, dim);
  default: return FLUIDS_ERROR_BADARG;
  }
}
