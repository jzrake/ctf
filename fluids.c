
#define FLUIDS_PRIVATE_DEFS
#include "fluids.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

static void _alloc_state(fluid_state *S, long modes);
static void _dealloc_state(fluid_state *S, long modes);
static int _getsetattrib(fluid_state *S, double *x, long flag, char op);

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
  default:
    return FLUIDS_ERROR_BADREQUEST;
  }
  _alloc_state(S, FLUIDS_CONSERVED | FLUIDS_PRIMITIVE);
  return 0;
}
int fluids_seteos(fluid_state *S, int eos)
{
  S->eos = eos;
  return 0;
}
int fluids_setcoordsystem(fluid_state *S, int coordsystem)
{
  S->coordsystem = coordsystem;
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

int fluids_getattrib(fluid_state *S, double *x, long flag)
{
  return _getsetattrib(S, x, flag, 'g');
}
int fluids_setattrib(fluid_state *S, double *x, long flag)
{
  return _getsetattrib(S, x, flag, 's');
}


int _getsetattrib(fluid_state *S, double *x, long flag, char op)
{
  double **a = NULL;
  int size = 0;

  switch (flag) {
  case FLUIDS_LOCATION:
    a = &S->location;
    size = 1;
    break;
  case FLUIDS_PASSIVE:
    a = &S->passive;
    size = S->npassive;
    break;
  case FLUIDS_CONSERVED:
    a = &S->conserved;
    size = S->nwaves;
    break;
  case FLUIDS_PRIMITIVE:
    a = &S->primitive;
    size = S->nwaves;
    break;
  case FLUIDS_MAGNETIC:
    a = &S->magnetic;
    size = 3;
    break;
  case FLUIDS_FOURVELOCITY:
    a = &S->fourvelocity;
    size = 3;
    break;
  default: // leave a == NULL
    break;
  }

  if (op == 'g') { // get
    if (a == NULL) {
      return FLUIDS_ERROR_BADREQUEST;
    }
    memcpy(x, *a, size * sizeof(double));
  }
  else if (op == 's') { // set
    *a = (double*) realloc(*a, size * sizeof(double));
    memcpy(*a, x, size * sizeof(double));
  }
  return 0;
}

void _alloc_state(fluid_state *S, long modes)
{
#define A(a,s,m) if(modes&m)S->a=(double*)realloc(S->a,(s)*sizeof(double))
  A(location, 3, FLUIDS_LOCATION);
  A(passive, S->npassive, FLUIDS_PASSIVE);
  A(conserved, S->nwaves, FLUIDS_CONSERVED);
  A(primitive, S->nwaves, FLUIDS_PRIMITIVE);
  A(magnetic, 3, FLUIDS_MAGNETIC);
  A(fourvelocity, 4, FLUIDS_FOURVELOCITY);
  A(flux[0], S->nwaves, FLUIDS_FLUX0);
  A(flux[1], S->nwaves, FLUIDS_FLUX1);
  A(flux[2], S->nwaves, FLUIDS_FLUX2);
  A(eigenvalues[0], S->nwaves, FLUIDS_EIGENVALUES0);
  A(eigenvalues[1], S->nwaves, FLUIDS_EIGENVALUES1);
  A(eigenvalues[2], S->nwaves, FLUIDS_EIGENVALUES2);
  A(leigenvectors[0], S->nwaves*S->nwaves, FLUIDS_LEIGENVECTORS0);
  A(leigenvectors[1], S->nwaves*S->nwaves, FLUIDS_LEIGENVECTORS1);
  A(leigenvectors[2], S->nwaves*S->nwaves, FLUIDS_LEIGENVECTORS2);
  A(leigenvectors[0], S->nwaves*S->nwaves, FLUIDS_REIGENVECTORS0);
  A(leigenvectors[1], S->nwaves*S->nwaves, FLUIDS_REIGENVECTORS1);
  A(leigenvectors[2], S->nwaves*S->nwaves, FLUIDS_REIGENVECTORS2);
#undef A
}

void _dealloc_state(fluid_state *S, long modes)
{
#define A(a,s,m) if(modes&m)S->a=(double*)realloc(S->a, 0)
  A(location, 3, FLUIDS_LOCATION);
  A(passive, S->npassive, FLUIDS_PASSIVE);
  A(conserved, S->nwaves, FLUIDS_CONSERVED);
  A(primitive, S->nwaves, FLUIDS_PRIMITIVE);
  A(magnetic, 3, FLUIDS_MAGNETIC);
  A(fourvelocity, 4, FLUIDS_FOURVELOCITY);
  A(flux[0], S->nwaves, FLUIDS_FLUX0);
  A(flux[1], S->nwaves, FLUIDS_FLUX1);
  A(flux[2], S->nwaves, FLUIDS_FLUX2);
  A(eigenvalues[0], S->nwaves, FLUIDS_EIGENVALUES0);
  A(eigenvalues[1], S->nwaves, FLUIDS_EIGENVALUES1);
  A(eigenvalues[2], S->nwaves, FLUIDS_EIGENVALUES2);
  A(leigenvectors[0], S->nwaves*S->nwaves, FLUIDS_LEIGENVECTORS0);
  A(leigenvectors[1], S->nwaves*S->nwaves, FLUIDS_LEIGENVECTORS1);
  A(leigenvectors[2], S->nwaves*S->nwaves, FLUIDS_LEIGENVECTORS2);
  A(leigenvectors[0], S->nwaves*S->nwaves, FLUIDS_REIGENVECTORS0);
  A(leigenvectors[1], S->nwaves*S->nwaves, FLUIDS_REIGENVECTORS1);
  A(leigenvectors[2], S->nwaves*S->nwaves, FLUIDS_REIGENVECTORS2);
#undef A
}



enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
enum { rho, pre, vx, vy, vz };             // Primitive
static int _nrhyd_c2p(fluid_state *S);
static int _nrhyd_p2c(fluid_state *S);
static int _nrhyd_cs2(fluid_state *S);
static int _nrhyd_update(fluid_state *S, long flags);

int fluids_update(fluid_state *S, long flags)
{
  switch (S->fluid) {
  case FLUIDS_NRHYD:
    return _nrhyd_update(S, flags);
  default:
    return FLUIDS_ERROR_BADREQUEST;
  }
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
  if (S->conserved == NULL) {
    return FLUIDS_ERROR_BADREQUEST;
  }
  _alloc_state(S, FLUIDS_PRIMITIVE);
  double gm1 = S->gammalawindex - 1.0;
  double *U = S->conserved;
  double *P = S->primitive;
  P[rho] =  U[ddd];
  P[pre] = (U[tau] - 0.5*(U[Sx]*U[Sx] + U[Sy]*U[Sy] + U[Sz]*U[Sz])/U[ddd])*gm1;
  P[vx ] =  U[Sx ] / U[ddd];
  P[vy ] =  U[Sy ] / U[ddd];
  P[vz ] =  U[Sz ] / U[ddd];
  return 0;
}

int _nrhyd_p2c(fluid_state *S)
{
  if (S->primitive == NULL) {
    return FLUIDS_ERROR_BADREQUEST;
  }
  _alloc_state(S, FLUIDS_CONSERVED);
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

int _nrhyd_cs2(fluid_state *S)
{
  double gm = S->gammalawindex;
  S->soundspeedsquared = gm * S->primitive[pre] / S->primitive[rho];
  return 0;
}

int _nrhyd_update(fluid_state *S, long modes)
{
  _alloc_state(S, modes);

  double *U = S->conserved;
  double *P = S->primitive;
  double *F = S->flux[0];
  double *G = S->flux[1];
  double *H = S->flux[2];

  if (modes & FLUIDS_FLUX0) {
    F[rho]  =  U[rho] * P[vx];
    F[tau]  = (U[tau] + P[pre])*P[vx];
    F[Sx]   =  U[Sx]  * P[vx] + P[pre];
    F[Sy]   =  U[Sy]  * P[vx];
    F[Sz]   =  U[Sz]  * P[vx];
  }
  if (modes & FLUIDS_FLUX1) {
    G[rho]  =  U[rho] * P[vy];
    G[tau]  = (U[tau] + P[pre])*P[vy];
    G[Sx]   =  U[Sx]  * P[vy];
    G[Sy]   =  U[Sy]  * P[vy] + P[pre];
    G[Sz]   =  U[Sz]  * P[vy];
  }
  if (modes & FLUIDS_FLUX2) {
    H[rho]  =  U[rho] * P[vz];
    H[tau]  = (U[tau] + P[pre])*P[vz];
    H[Sx]   =  U[Sx]  * P[vz];
    H[Sy]   =  U[Sy]  * P[vz];
    H[Sz]   =  U[Sz]  * P[vz] + P[pre];
  }

  if (modes & FLUIDS_EIGENVALUES0) {
    _nrhyd_cs2(S);
    double a = sqrt(S->soundspeedsquared);
    S->eigenvalues[0][0] = P[vx] - a;
    S->eigenvalues[0][1] = P[vx];
    S->eigenvalues[0][2] = P[vx];
    S->eigenvalues[0][3] = P[vx];
    S->eigenvalues[0][4] = P[vx] + a;
  }
  if (modes & FLUIDS_EIGENVALUES1) {
    _nrhyd_cs2(S);
    double a = sqrt(S->soundspeedsquared);
    S->eigenvalues[1][0] = P[vy] - a;
    S->eigenvalues[1][1] = P[vy];
    S->eigenvalues[1][2] = P[vy];
    S->eigenvalues[1][3] = P[vy];
    S->eigenvalues[1][4] = P[vy] + a;
  }
  if (modes & FLUIDS_EIGENVALUES2) {
    _nrhyd_cs2(S);
    double a = sqrt(S->soundspeedsquared);
    S->eigenvalues[2][0] = P[vz] - a;
    S->eigenvalues[2][1] = P[vz];
    S->eigenvalues[2][2] = P[vz];
    S->eigenvalues[2][3] = P[vz];
    S->eigenvalues[2][4] = P[vz] + a;
  }
  return 0;
}
