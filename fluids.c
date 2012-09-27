

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
static void _alloc_cache(fluids_cache *C, int op);

static int _nrhyd_c2p(fluids_state *S, double *U);
static int _nrhyd_p2c(fluids_state *S);
static int _nrhyd_update(fluids_state *S, long flags);
static void _nrhyd_cs2(fluids_state *S, double *cs2);
static void _nrhyd_eigenvec(fluids_state *S, int dim, int doleft, int dorght);
static void _nrhyd_jacobian(fluids_state *S, int dim);


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
    .needsupdateflags = FLUIDS_FLAGSALL,
    .state = NULL,
  } ;
  *C = cache;
  return C;
} ;

int fluids_cache_del(fluids_cache *C)
{
  if (C != NULL) {
    _alloc_cache(C, DEALLOC);
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
    D->cacheflags = FLUIDS_FLAGSALL;
    break;
  case FLUIDS_NRHYD:
    D->nprimitive = 5;
    D->npassive = 0;
    D->ngravity = 0;
    D->nmagnetic = 0;
    D->nlocation = 0;
    D->cacheflags = FLUIDS_FLAGSALL;
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
    .cache = NULL,
    .descr = NULL,
  } ;
  *S = state;
  return S;
}

int fluids_state_del(fluids_state *S)
{
  if (S != NULL) {
    fluids_state_erasecache(S);
    free(S->primitive);
    free(S->passive);
    free(S->gravity);
    free(S->magnetic);
    free(S->location);
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
  return 0;
}

int fluids_state_resetcache(fluids_state *S)
{
  if (S->cache == NULL) {
    S->cache = fluids_cache_new();
    S->cache->state = S;
    _alloc_cache(S->cache, ALLOC);
  }
  S->cache->needsupdateflags = FLUIDS_FLAGSALL;
  return 0;
}

int fluids_state_erasecache(fluids_state *S)
{
  fluids_cache_del(S->cache);
  S->cache = NULL;
  return 0;
}

int fluids_state_getattr(fluids_state *S, double *x, long flag)
{
  _getsetstateattr(S, x, flag, 'g');
  return 0;
}

int fluids_state_setattr(fluids_state *S, double *x, long flag)
{
  fluids_state_resetcache(S);
  _getsetstateattr(S, x, flag, 's');
  return 0;
}

int fluids_state_fromcons(fluids_state *S, double *U, int cachebehavior)
/*
 * The fluid state `S` must be complete except for its primitive field, which
 * will be over-written, derived from the conserved state `U`. If the inversion
 * from conserved to primitive requires a rootfinder, then the existing
 * primitive in `S` may be used a guess value.
 *
 * `cachebehavior` is one of:
 *
 *     NOTOUCH ... will be left alone
 *     RESET   ... (recommended) then the state's cache will be reset
 *     ERASE   ... will be erased 
 */
{
  _nrhyd_c2p(S, U);
  switch (cachebehavior) {
  case FLUIDS_CACHE_NOTOUCH:
    break;    
  case FLUIDS_CACHE_RESET:
    fluids_state_resetcache(S);
    break;
  case FLUIDS_CACHE_ERASE:
    fluids_state_erasecache(S);
    break;
  default:
    fluids_state_resetcache(S);
    break;
  }
  return 0;
}

int fluids_state_derive(fluids_state *S, double *x, int flag)
/*
 * If `x` is not NULL, then `flag` must represent only a single field, and that
 * field will be (deep) copied into the array `x`. If `x` is NULL then it will
 * not be used, `flag` may reference many fields which will all be updated.
 */
{
  if (S->cache == NULL) {
    fluids_state_resetcache(S);
  }
  _nrhyd_update(S, flag);
  if (x != NULL) {
    _getsetcacheattr(S->cache, x, flag, 'g');
  }
  return 0;
}





int _getsetcacheattr(fluids_cache *C, double *x, long flag, char op)
{
  double *a = NULL;
  int size = 0;
  int np = C->state->descr->nprimitive;
#define CASE(f,m,s)case FLUIDS_##f: a = m; size = s; break
  switch (flag) {
    CASE(CONSERVED, C->conserved, np);
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

void _alloc_cache(fluids_cache *C, int op)
{
#define A(a,s,m) do {                                           \
    if (C->state->descr->cacheflags & m) {                      \
      if (op == ALLOC) {                                        \
        C->a = (double*) realloc(C->a, (s)*sizeof(double));     \
      }                                                         \
      else if (op == DEALLOC) {                                 \
        free(C->a);                                             \
        C->a = NULL;                                            \
      }                                                         \
    }                                                           \
  } while (0)
  int np = C->state->descr->nprimitive;
  A(conserved, np, FLUIDS_CONSERVED);
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
}










int _nrhyd_c2p(fluids_state *S, double *U)
{
  double gm1 = S->descr->gammalawindex - 1.0;
  double *P = S->primitive;
  P[rho] =  U[ddd];
  P[pre] = (U[tau] - 0.5*(U[Sx]*U[Sx] + U[Sy]*U[Sy] + U[Sz]*U[Sz])/U[ddd])*gm1;
  P[vx]  =  U[Sx] / U[ddd];
  P[vy]  =  U[Sy] / U[ddd];
  P[vz]  =  U[Sz] / U[ddd];
  return 0;
}

int _nrhyd_p2c(fluids_state *S)
{
  double gm1 = S->descr->gammalawindex - 1.0;
  double *U = S->cache->conserved;
  double *P = S->primitive;
  U[ddd] = P[rho];
  U[Sx]  = P[rho] * P[vx];
  U[Sy]  = P[rho] * P[vy];
  U[Sz]  = P[rho] * P[vz];
  U[tau] = P[rho] * 0.5*(P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz]) + P[pre]/gm1;
  return 0;
}

void _nrhyd_cs2(fluids_state *S, double *cs2)
{
  double gm = S->descr->gammalawindex;
  *cs2 = gm * S->primitive[pre] / S->primitive[rho];
}

int _nrhyd_update(fluids_state *S, long modes)
/*
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
 * needsupdateflags becomes: 000010000 (needsupdateflags &= !modes)
 */
{
  fluids_cache *C = S->cache;
  double *U = C->conserved;
  double *P = S->primitive;
  double a=0.0, cs2=0.0;

  if (modes & FLUIDS_FLAGSALL) {
    modes |= FLUIDS_CONSERVED; // conserved quantities are used for everything
  }
  if (modes & FLUIDS_EVALSALL) {
    modes |= FLUIDS_SOUNDSPEEDSQUARED; // cs is used for eigenvalues
  }
  modes &= C->needsupdateflags;

  /*
    printf("update conserved? %s\n", (modes & FLUIDS_CONSERVED) ? "yes" : "no");
    printf("update flux0? %s\n", (modes & FLUIDS_FLUX0) ? "yes" : "no");
  */

  if (modes & FLUIDS_CONSERVED) {
    _nrhyd_p2c(S);
  }

  if (modes & FLUIDS_FLUX0) {
    C->flux[0][rho] = U[rho] * P[vx];
    C->flux[0][tau] = (U[tau] + P[pre]) * P[vx];
    C->flux[0][Sx] = U[Sx] * P[vx] + P[pre];
    C->flux[0][Sy] = U[Sy] * P[vx];
    C->flux[0][Sz] = U[Sz] * P[vx];
  }
  if (modes & FLUIDS_FLUX1) {
    C->flux[1][rho] = U[rho] * P[vy];
    C->flux[1][tau] = (U[tau] + P[pre]) * P[vy];
    C->flux[1][Sx] = U[Sx] * P[vy];
    C->flux[1][Sy] = U[Sy] * P[vy] + P[pre];
    C->flux[1][Sz] = U[Sz] * P[vy];
  }
  if (modes & FLUIDS_FLUX2) {
    C->flux[2][rho] = U[rho] * P[vz];
    C->flux[2][tau] = (U[tau] + P[pre]) * P[vz];
    C->flux[2][Sx] = U[Sx] * P[vz];
    C->flux[2][Sy] = U[Sy] * P[vz];
    C->flux[2][Sz] = U[Sz] * P[vz] + P[pre];
  }

  if (modes & FLUIDS_SOUNDSPEEDSQUARED) {
    _nrhyd_cs2(S, &cs2);
    a = sqrt(cs2);
    C->soundspeedsquared = cs2;
  }

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

  if (modes & (FLUIDS_LEVECS0 | FLUIDS_REVECS0)) {
    _nrhyd_eigenvec(S, 0,
                    modes & FLUIDS_LEVECS0,
                    modes & FLUIDS_REVECS0);
  }
  if (modes & (FLUIDS_LEVECS1 | FLUIDS_REVECS1)) {
    _nrhyd_eigenvec(S, 1,
                    modes & FLUIDS_LEVECS1,
                    modes & FLUIDS_REVECS1);
  }
  if (modes & (FLUIDS_LEVECS2 | FLUIDS_REVECS2)) {
    _nrhyd_eigenvec(S, 2,
                    modes & FLUIDS_LEVECS2,
                    modes & FLUIDS_REVECS2);
  }

  if (modes & FLUIDS_JACOBIAN0) {
    _nrhyd_jacobian(S, 0);
  }
  if (modes & FLUIDS_JACOBIAN1) {
    _nrhyd_jacobian(S, 1);
  }
  if (modes & FLUIDS_JACOBIAN2) {
    _nrhyd_jacobian(S, 2);
  }

  C->needsupdateflags &= BITWISENOT(modes);
  return 0;
}



void _nrhyd_eigenvec(fluids_state *S, int dim, int doleft, int dorght)
{
  int v1=0, v2=0, v3=0;
  switch (dim) {
  case 0:
    v1=vx; v2=vy; v3=vz;
    break;
  case 1:
    v1=vy; v2=vz; v3=vx;
    break;
  case 2:
    v1=vz; v2=vx; v3=vy;
    break;
  }
  double *P = S->primitive;
  double *U = S->cache->conserved;
  double *L = S->cache->leigenvectors[dim];
  double *R = S->cache->reigenvectors[dim];
  double gm = S->descr->gammalawindex;
  double gm1 = gm - 1.0;
  double u = P[v1];
  double v = P[v2];
  double w = P[v3];
  double V2 = u*u + v*v + w*w;
  double a = sqrt(gm * P[pre] / P[rho]);
  double H = (U[tau] + P[pre]) / P[rho];

  // Toro Equation 3.82
  // ---------------------------------------------------------------------------
  double R_[5][5] =
    { {       1,      1,      0,      0,     1   },
      {     u-a,      u,      0,      0,     u+a },
      {       v,      v,      1,      0,     v   },
      {       w,      w,      0,      1,     w   },
      { H - u*a, 0.5*V2,      v,      w, H + u*a } };

  // Toro Equation 3.83 up to (gam - 1) / (2*a^2)
  // ---------------------------------------------------------------------------
  double L_[5][5] =
    { {    H + (a/gm1)*(u-a),  -(u+a/gm1),        -v,        -w,  1 },
      { -2*H + (4/gm1)*(a*a),         2*u,       2*v,       2*w, -2 },
      {         -2*v*a*a/gm1,           0, 2*a*a/gm1,         0,  0 },
      {         -2*w*a*a/gm1,           0,         0, 2*a*a/gm1,  0 },
      {    H - (a/gm1)*(u+a),  -(u-a/gm1),        -v,        -w,  1 } };

  // Permute the eigenvectors according to the direction:
  // ---------------------------------------------------------------------------
  // L' = L P
  // R' = P^{-1} R
  // ---------------------------------------------------------------------------
  double P1[5][5] =
    { { 1, 0, 0, 0, 0 },
      { 0, 1, 0, 0, 0 },
      { 0, 0, 1, 0, 0 },
      { 0, 0, 0, 1, 0 },
      { 0, 0, 0, 0, 1 } };

  double P2[5][5] =
    { { 1, 0, 0, 0, 0 },
      { 0, 0, 1, 0, 0 },
      { 0, 0, 0, 1, 0 },
      { 0, 1, 0, 0, 0 },
      { 0, 0, 0, 0, 1 } };

  double P3[5][5] =
    { { 1, 0, 0, 0, 0 },
      { 0, 0, 0, 1, 0 },
      { 0, 1, 0, 0, 0 },
      { 0, 0, 1, 0, 0 },
      { 0, 0, 0, 0, 1 } };

  switch (dim) {
  case 0:
    if (doleft) matrix_matrix_product(L_[0], P1[0], L, 5, 5, 5);
    if (dorght) matrix_matrix_product(P1[0], R_[0], R, 5, 5, 5);
    break;
  case 1:
    if (doleft) matrix_matrix_product(L_[0], P2[0], L, 5, 5, 5);
    if (dorght) matrix_matrix_product(P3[0], R_[0], R, 5, 5, 5);
    break;
  case 2:
    if (doleft) matrix_matrix_product(L_[0], P3[0], L, 5, 5, 5);
    if (dorght) matrix_matrix_product(P2[0], R_[0], R, 5, 5, 5);
    break;
  }

  // Replace the term in eqn 3.83 : (gam - 1) / (2*a^2)
  // ---------------------------------------------------------------------------
  double norm = gm1 / (2*a*a);
  for (int i=0; i<25; ++i) {
    L[i] *= norm;
  }
}

void _nrhyd_jacobian(fluids_state *S, int dim)
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

  memcpy(S->cache->jacobian[dim], A[0], 25*sizeof(double));
}
