
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
  long modes = 0;

  modes |= FLUIDS_CONSERVED;
  modes |= FLUIDS_PRIMITIVE;
  modes |= FLUIDS_FLUXALL;
  modes |= FLUIDS_EIGENVALUESALL;

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
  _alloc_state(S, modes);
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
    CASE(EIGENVALUES0, S->eigenvalues[0], S->nwaves);
    CASE(EIGENVALUES1, S->eigenvalues[1], S->nwaves);
    CASE(EIGENVALUES2, S->eigenvalues[2], S->nwaves);
    CASE(LEIGENVECTORS0, S->leigenvectors[0], S->nwaves*S->nwaves);
    CASE(LEIGENVECTORS1, S->leigenvectors[1], S->nwaves*S->nwaves);
    CASE(LEIGENVECTORS2, S->leigenvectors[2], S->nwaves*S->nwaves);
    CASE(REIGENVECTORS0, S->reigenvectors[0], S->nwaves*S->nwaves);
    CASE(REIGENVECTORS1, S->reigenvectors[1], S->nwaves*S->nwaves);
    CASE(REIGENVECTORS2, S->reigenvectors[2], S->nwaves*S->nwaves);
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
#define A(a,s,m) if(modes&m)free(S->a)//=(double*)realloc(S->a, 0)
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
static double _nrhyd_cs2(fluid_state *S);
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
  double *U = S->conserved;
  double *P = S->primitive;
  double cs2 = _nrhyd_cs2(S);
  double a = sqrt(cs2);

  if (modes & FLUIDS_SOUNDSPEEDSQUARED) {
    S->soundspeedsquared = cs2;
  }

  if (modes & FLUIDS_FLUX0) {
    S->flux[0][rho]  =  U[rho] * P[vx];
    S->flux[0][tau]  = (U[tau] + P[pre])*P[vx];
    S->flux[0][Sx]   =  U[Sx]  * P[vx] + P[pre];
    S->flux[0][Sy]   =  U[Sy]  * P[vx];
    S->flux[0][Sz]   =  U[Sz]  * P[vx];
  }
  if (modes & FLUIDS_FLUX1) {
    S->flux[1][rho]  =  U[rho] * P[vy];
    S->flux[1][tau]  = (U[tau] + P[pre])*P[vy];
    S->flux[1][Sx]   =  U[Sx]  * P[vy];
    S->flux[1][Sy]   =  U[Sy]  * P[vy] + P[pre];
    S->flux[1][Sz]   =  U[Sz]  * P[vy];
  }
  if (modes & FLUIDS_FLUX2) {
    S->flux[2][rho]  =  U[rho] * P[vz];
    S->flux[2][tau]  = (U[tau] + P[pre])*P[vz];
    S->flux[2][Sx]   =  U[Sx]  * P[vz];
    S->flux[2][Sy]   =  U[Sy]  * P[vz];
    S->flux[2][Sz]   =  U[Sz]  * P[vz] + P[pre];
  }

  if (modes & FLUIDS_EIGENVALUES0) {
    S->eigenvalues[0][0] = P[vx] - a;
    S->eigenvalues[0][1] = P[vx];
    S->eigenvalues[0][2] = P[vx];
    S->eigenvalues[0][3] = P[vx];
    S->eigenvalues[0][4] = P[vx] + a;
  }
  if (modes & FLUIDS_EIGENVALUES1) {
    S->eigenvalues[1][0] = P[vy] - a;
    S->eigenvalues[1][1] = P[vy];
    S->eigenvalues[1][2] = P[vy];
    S->eigenvalues[1][3] = P[vy];
    S->eigenvalues[1][4] = P[vy] + a;
  }
  if (modes & FLUIDS_EIGENVALUES2) {
    S->eigenvalues[2][0] = P[vz] - a;
    S->eigenvalues[2][1] = P[vz];
    S->eigenvalues[2][2] = P[vz];
    S->eigenvalues[2][3] = P[vz];
    S->eigenvalues[2][4] = P[vz] + a;
  }
  return 0;
}
