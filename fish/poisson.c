

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define FISH_PRIVATE_DEFS
#define FLUIDS_INDEX_VARS
#include "fish.h"

#ifdef USE_FFTW
#include <fftw3.h>
#endif // USE_FFTW
#ifndef M_PI
#define M_PI (4*atan(1.0))
#endif // M_PI

#define CHECK(c,m) do{if(!(c)){B->error=m;return -1;}B->error=NULL;}while(0)

int fish_block_solvepoisson(fish_block *B)
/* -----------------------------------------------------------------------------
 *
 * Solve the Poisson equation in 1d assuming periodic boundary conditions.
 *
 * -----------------------------------------------------------------------------
 */
{
#ifndef USE_FFTW
  CHECK(0, "FFTW must be enabled to solve a poisson equation with DFT's");
#else
  CHECK(B->allocated, "block must already be allocated");
  CHECK(B->descr, "block needs a fluid descriptor");
  CHECK(B->rank == 1, "block must be 1d, higher dimensions not yet supported");

  int fluid_id;
  fluids_descr_getfluid(B->descr, &fluid_id);
  CHECK(fluid_id == FLUIDS_GRAVS, "fluid id must be GRAVS");

  int Ng = B->guard;
  int Nx = B->size[0];

  fftw_complex *Rhox = (fftw_complex*) fftw_malloc(Nx * sizeof(fftw_complex));
  fftw_complex *Phix = (fftw_complex*) fftw_malloc(Nx * sizeof(fftw_complex));
  fftw_complex *Gphx = (fftw_complex*) fftw_malloc(Nx * sizeof(fftw_complex));
  fftw_complex *Rhok = (fftw_complex*) fftw_malloc(Nx * sizeof(fftw_complex));
  fftw_complex *Phik = (fftw_complex*) fftw_malloc(Nx * sizeof(fftw_complex));
  fftw_complex *Gphk = (fftw_complex*) fftw_malloc(Nx * sizeof(fftw_complex));

  fftw_plan fwd = fftw_plan_dft_1d(Nx, Rhox, Rhok, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan revphi = fftw_plan_dft_1d(Nx, Phik, Phix, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan revgph = fftw_plan_dft_1d(Nx, Gphk, Gphx, FFTW_BACKWARD, FFTW_ESTIMATE);

  for (int i=0; i<Nx; ++i) {
    double P[5];
    fluids_state_getattr(B->fluid[i+Ng], P, FLUIDS_PRIMITIVE);
    Rhox[i][0] = P[0];
    Rhox[i][1] = 0.0;
  }
  fftw_execute(fwd);

  double rhobar = Rhok[0][0] / Nx;
  fluids_descr_setrhobar(B->descr, rhobar);

  Phik[0][0] = 0.0; // phibar = 0.0
  Phik[0][1] = 0.0;
  Gphk[0][0] = 0.0; // gphbar = 0.0
  Gphk[0][1] = 0.0;

  /*
   * gph[k] = I * k * (phi[k].re + I * phi[k].im)
   */
  for (int i=1; i<Nx; ++i) {
    double k = 2*M_PI * (i < Nx/2 ? i : i-Nx);
    Phik[i][0] = -Rhok[i][0] / (k*k);
    Phik[i][1] = -Rhok[i][1] / (k*k);
    Gphk[i][0] = -Phik[i][1] * k;
    Gphk[i][1] = +Phik[i][0] * k;
  }
  fftw_execute(revphi);
  fftw_execute(revgph);

  for (int i=0; i<Nx; ++i) {
    double G[4];
    G[phi] = Phix[i][0] / Nx;
    G[gph] = Gphx[i][0] / Nx;
    fluids_state_setattr(B->fluid[i+Ng], G, FLUIDS_GRAVITY);
  }

  fftw_destroy_plan(fwd);
  fftw_destroy_plan(revphi);
  fftw_destroy_plan(revgph);
  fftw_free(Rhox);
  fftw_free(Phix);
  fftw_free(Gphx);
  fftw_free(Rhok);
  fftw_free(Phik);
  fftw_free(Gphk);

  return 0;
#endif // USE_FFTW
}

