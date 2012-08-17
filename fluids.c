

#define FLUIDS_PRIVATE_DEFS
#define FLUIDS_INDEX_VARS
#include "fluids.h"
#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* http://en.wikipedia.org/wiki/Bitwise_operation#NOT */
#define BITWISENOT(x) -(x) - 1
#define ALLOC 1
#define DEALLOC 0

static void _alloc_state(fluid_state *S, long modes, int op);
static int _getsetattrib(fluid_state *S, double *x, long flag, char op);
static int _nrhyd_c2p(fluid_state *S);
static int _nrhyd_p2c(fluid_state *S);
static double _nrhyd_cs2(fluid_state *S);
static int _nrhyd_update(fluid_state *S, long flags);
static void _nrhyd_eigenvec(fluid_state *S, int dim, int doleft, int dorght);
static void _nrhyd_jacobian(fluid_state *S, int dim);

fluid_state *fluids_new()
{
  fluid_state *S = (fluid_state*) malloc(sizeof(fluid_state));
  fluid_state state = {
    .fluid = FLUIDS_SCALAR_ADVECTION,
    .eos = FLUIDS_EOS_GAMMALAW,
    .coordsystem = FLUIDS_COORD_CARTESIAN,
    .nwaves = 0,
    .npassive = 0,
    .needsupdateflags = FLUIDS_FLAGSALL,
    .location = NULL,
    .passive = NULL,
    .conserved = NULL,
    .primitive = NULL,
    .magnetic = NULL,
    .fourvelocity = NULL,
    .flux = {NULL, NULL, NULL},
    .eigenvalues = {NULL, NULL, NULL},
    .leigenvectors = {NULL, NULL, NULL},
    .reigenvectors = {NULL, NULL, NULL},
    .soundspeedsquared = 0.0,
    .temperature = 0.0,
    .specificenthalpy = 0.0,
    .specificinternal = 0.0,
    .gammalawindex = 1.4,
  } ;
  *S = state;
  return S;
}

int fluids_del(fluid_state *S)
{
  _alloc_state(S, FLUIDS_FLAGSALL, DEALLOC);
  free(S);
  return 0;
}

int fluids_setfluid(fluid_state *S, int fluid)
{
  S->needsupdateflags = FLUIDS_FLAGSALL;
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
  default:
    return FLUIDS_ERROR_BADREQUEST;
  }
  return 0;
}

int fluids_seteos(fluid_state *S, int eos)
{
  S->needsupdateflags = FLUIDS_FLAGSALL;
  S->eos = eos;
  return 0;
}

int fluids_setcoordsystem(fluid_state *S, int coordsystem)
{
  S->needsupdateflags = FLUIDS_FLAGSALL;
  S->coordsystem = coordsystem;
  return 0;
}

int fluids_setnpassive(fluid_state *S, int n)
{
  S->needsupdateflags = FLUIDS_FLAGSALL;
  S->npassive = n;
  return 0;
}

int fluids_getattrib(fluid_state *S, double *x, long flag)
{
  fluids_update(S, flag);
  return _getsetattrib(S, x, flag, 'g');
}

int fluids_setattrib(fluid_state *S, double *x, long flag)
{
  S->needsupdateflags = FLUIDS_FLAGSALL; // most conservative invalidation rule
  return _getsetattrib(S, x, flag, 's');
}

int _getsetattrib(fluid_state *S, double *x, long flag, char op)
{
  double *a = NULL;
  int size = 0;

#define CASE(f,m,s)case FLUIDS_##f: a = m; size = s; break
  switch (flag) {
    CASE(LOCATION, S->location, 3);
    CASE(PASSIVE, S->passive, S->npassive);
    CASE(CONSERVED, S->conserved, S->nwaves);
    CASE(PRIMITIVE, S->primitive, S->nwaves);
    CASE(MAGNETIC, S->magnetic, 3);
    CASE(FOURVELOCITY, S->fourvelocity, 4);
    CASE(FLUX0, S->flux[0], S->nwaves);
    CASE(FLUX1, S->flux[1], S->nwaves);
    CASE(FLUX2, S->flux[2], S->nwaves);
    CASE(EVAL0, S->eigenvalues[0], S->nwaves);
    CASE(EVAL1, S->eigenvalues[1], S->nwaves);
    CASE(EVAL2, S->eigenvalues[2], S->nwaves);
    CASE(LEVECS0, S->leigenvectors[0], S->nwaves*S->nwaves);
    CASE(LEVECS1, S->leigenvectors[1], S->nwaves*S->nwaves);
    CASE(LEVECS2, S->leigenvectors[2], S->nwaves*S->nwaves);
    CASE(REVECS0, S->reigenvectors[0], S->nwaves*S->nwaves);
    CASE(REVECS1, S->reigenvectors[1], S->nwaves*S->nwaves);
    CASE(REVECS2, S->reigenvectors[2], S->nwaves*S->nwaves);
    CASE(JACOBIAN0, S->jacobian[0], S->nwaves*S->nwaves);
    CASE(JACOBIAN1, S->jacobian[1], S->nwaves*S->nwaves);
    CASE(JACOBIAN2, S->jacobian[2], S->nwaves*S->nwaves);
    CASE(SOUNDSPEEDSQUARED, &S->soundspeedsquared, 1);
    CASE(TEMPERATURE, &S->temperature, 1);
    CASE(SPECIFICENTHALPY, &S->specificenthalpy, 1);
    CASE(SPECIFICINTERNAL, &S->specificinternal, 1);
    CASE(GAMMALAWINDEX, &S->gammalawindex, 1);
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

void _alloc_state(fluid_state *S, long modes, int op)
{
#define A(a,s,m) do {							\
    if (modes & m) {							\
      if (op == ALLOC) {						\
	S->a = (double*) realloc(S->a,(s)*sizeof(double));		\
      }									\
      else if (op == DEALLOC) {						\
	free(S->a);							\
      }									\
    }									\
  } while (0)								\

  A(location, 3, FLUIDS_LOCATION);
  A(passive, S->npassive, FLUIDS_PASSIVE);
  A(conserved, S->nwaves, FLUIDS_CONSERVED);
  A(primitive, S->nwaves, FLUIDS_PRIMITIVE);
  A(magnetic, 3, FLUIDS_MAGNETIC);
  A(fourvelocity, 4, FLUIDS_FOURVELOCITY);
  A(flux[0], S->nwaves, FLUIDS_FLUX0);
  A(flux[1], S->nwaves, FLUIDS_FLUX1);
  A(flux[2], S->nwaves, FLUIDS_FLUX2);
  A(eigenvalues[0], S->nwaves, FLUIDS_EVAL0);
  A(eigenvalues[1], S->nwaves, FLUIDS_EVAL1);
  A(eigenvalues[2], S->nwaves, FLUIDS_EVAL2);
  A(leigenvectors[0], S->nwaves*S->nwaves, FLUIDS_LEVECS0);
  A(leigenvectors[1], S->nwaves*S->nwaves, FLUIDS_LEVECS1);
  A(leigenvectors[2], S->nwaves*S->nwaves, FLUIDS_LEVECS2);
  A(reigenvectors[0], S->nwaves*S->nwaves, FLUIDS_REVECS0);
  A(reigenvectors[1], S->nwaves*S->nwaves, FLUIDS_REVECS1);
  A(reigenvectors[2], S->nwaves*S->nwaves, FLUIDS_REVECS2);
  A(jacobian[0], S->nwaves*S->nwaves, FLUIDS_JACOBIAN0);
  A(jacobian[1], S->nwaves*S->nwaves, FLUIDS_JACOBIAN1);
  A(jacobian[2], S->nwaves*S->nwaves, FLUIDS_JACOBIAN2);
#undef A
}

int fluids_update(fluid_state *S, long flags)
{
  switch (S->fluid) {
  case FLUIDS_NRHYD:
    return _nrhyd_update(S, flags);
  default:
    return FLUIDS_ERROR_BADREQUEST;
  }
}

int fluids_resetcache(fluid_state *S)
{
  S->needsupdateflags = FLUIDS_FLAGSALL;
  return 0;
}

int fluids_getlastupdate(fluid_state *S, long *flags)
{
  *flags = S->lastupdatedflags;
  return 0;
}

int fluids_alloc(fluid_state *S, long flags)
{
  _alloc_state(S, flags, ALLOC);
  return 0;
}

int fluids_c2p(fluid_state *S)
{
  switch (S->fluid) {
  case FLUIDS_NRHYD:
    return _nrhyd_c2p(S);
  default:
    return FLUIDS_ERROR_BADREQUEST;
  }
}

int fluids_p2c(fluid_state *S)
{
  switch (S->fluid) {
  case FLUIDS_NRHYD:
    return _nrhyd_p2c(S);
  default:
    return FLUIDS_ERROR_BADREQUEST;
  }
}

int _nrhyd_c2p(fluid_state *S)
{
  double gm1 = S->gammalawindex - 1.0;
  double *U = S->conserved;
  double *P = S->primitive;
  P[rho] =  U[ddd];
  P[pre] = (U[tau] - 0.5*(U[Sx]*U[Sx] + U[Sy]*U[Sy] + U[Sz]*U[Sz])/U[ddd])*gm1;
  P[vx]  =  U[Sx] / U[ddd];
  P[vy]  =  U[Sy] / U[ddd];
  P[vz]  =  U[Sz] / U[ddd];
  return 0;
}

int _nrhyd_p2c(fluid_state *S)
{
  double gm1 = S->gammalawindex - 1.0;
  double *U = S->conserved;
  double *P = S->primitive;
  U[ddd] = P[rho];
  U[Sx]  = P[rho] * P[vx];
  U[Sy]  = P[rho] * P[vy];
  U[Sz]  = P[rho] * P[vz];
  U[tau] = P[rho] * 0.5*(P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz]) + P[pre]/gm1;
  return 0;
}

double _nrhyd_cs2(fluid_state *S)
{
  double gm = S->gammalawindex;
  return gm * S->primitive[pre] / S->primitive[rho];
}

int _nrhyd_update(fluid_state *S, long modes)
{
  /* It's only necessary to update the fields which were requested in `modes`,
   * but also have their needsupdate bit enabled.
   *
   * modes:                    001101011
   * needsupdateflags:         000010001
   * modes becomes:            000000001 (modes &= S->needsupdateflags)
   *
   * After the update is finished, any bit in needsupdateflags which is enabled
   * in modes should be set to zero. In other words:
   *
   * needsupdateflags:         000010001
   * modes:                    000000001
   * needsupdateflags becomes: 000010000 (needsupdateflags &= !modes)
   */
  modes &= S->needsupdateflags;

  double *U = S->conserved;
  double *P = S->primitive;
  double a=0.0, cs2=0.0;

  if (modes & FLUIDS_FLUX0) {
    S->flux[0][rho] = U[rho] * P[vx];
    S->flux[0][tau] = (U[tau] + P[pre]) * P[vx];
    S->flux[0][Sx] = U[Sx] * P[vx] + P[pre];
    S->flux[0][Sy] = U[Sy] * P[vx];
    S->flux[0][Sz] = U[Sz] * P[vx];
  }
  if (modes & FLUIDS_FLUX1) {
    S->flux[1][rho] = U[rho] * P[vy];
    S->flux[1][tau] = (U[tau] + P[pre]) * P[vy];
    S->flux[1][Sx] = U[Sx] * P[vy];
    S->flux[1][Sy] = U[Sy] * P[vy] + P[pre];
    S->flux[1][Sz] = U[Sz] * P[vy];
  }
  if (modes & FLUIDS_FLUX2) {
    S->flux[2][rho] = U[rho] * P[vz];
    S->flux[2][tau] = (U[tau] + P[pre]) * P[vz];
    S->flux[2][Sx] = U[Sx] * P[vz];
    S->flux[2][Sy] = U[Sy] * P[vz];
    S->flux[2][Sz] = U[Sz] * P[vz] + P[pre];
  }

  if (modes & (FLUIDS_EVALSALL | FLUIDS_SOUNDSPEEDSQUARED)) {
    cs2 = _nrhyd_cs2(S);
    a = sqrt(cs2);
  }

  if (modes & FLUIDS_SOUNDSPEEDSQUARED) {
    S->soundspeedsquared = cs2;
  }

  if (modes & FLUIDS_EVAL0) {
    S->eigenvalues[0][0] = P[vx] - a;
    S->eigenvalues[0][1] = P[vx];
    S->eigenvalues[0][2] = P[vx];
    S->eigenvalues[0][3] = P[vx];
    S->eigenvalues[0][4] = P[vx] + a;
  }
  if (modes & FLUIDS_EVAL1) {
    S->eigenvalues[1][0] = P[vy] - a;
    S->eigenvalues[1][1] = P[vy];
    S->eigenvalues[1][2] = P[vy];
    S->eigenvalues[1][3] = P[vy];
    S->eigenvalues[1][4] = P[vy] + a;
  }
  if (modes & FLUIDS_EVAL2) {
    S->eigenvalues[2][0] = P[vz] - a;
    S->eigenvalues[2][1] = P[vz];
    S->eigenvalues[2][2] = P[vz];
    S->eigenvalues[2][3] = P[vz];
    S->eigenvalues[2][4] = P[vz] + a;
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

  S->lastupdatedflags = modes;
  S->needsupdateflags &= BITWISENOT(modes);
  return 0;
}



void _nrhyd_eigenvec(fluid_state *S, int dim, int doleft, int dorght)
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
  double *U = S->conserved;
  double *P = S->primitive;
  double *L = S->leigenvectors[dim];
  double *R = S->reigenvectors[dim];
  double gm = S->gammalawindex;
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

void _nrhyd_jacobian(fluid_state *S, int dim)
{
  double gm = S->gammalawindex;
  double g1 = gm - 1.0;
  double g2 = gm - 2.0;

 // inputs are Mara convention: {D,E,px,py,pz}
  double *U = S->conserved;
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

  memcpy(S->jacobian[dim], A[0], 25*sizeof(double));
}

