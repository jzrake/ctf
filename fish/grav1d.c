

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include "fish.h"

static struct fluids_state **fluid;
static fluids_descr *descr;
static int TotalZones;
static int NumGhostZones;
static double DomainLength = 1.0;
static double dx;


static void solve_poisson(double *Rho, double *Phi, double *Gph, double *rhobar);
static void timederiv(double *L);

void fish_grav1d_init(int N)
{
  NumGhostZones = 3;
  TotalZones = N + 2*NumGhostZones;

  descr = fluids_descr_new();
  fluids_descr_setfluid(descr, FLUIDS_GRAVS);
  fluids_descr_setgamma(descr, 1.4);
  fluids_descr_seteos(descr, FLUIDS_EOS_GAMMALAW);

  fluid = (fluids_state**) malloc(TotalZones * sizeof(fluids_state*));
  dx = DomainLength / (TotalZones - NumGhostZones);

  for (int n=0; n<TotalZones; ++n) {
    fluid[n] = fluids_state_new();
    fluids_state_setdescr(fluid[n], descr);
  }
}

void fish_grav1d_finalize()
{
  for (int n=0; n<TotalZones; ++n) {
    fluids_state_del(fluid[n]);
  }
  fluids_descr_del(descr);
  free(fluid);
}

void fish_grav1d_advance(double dt)
{
  double *L  = (double*) malloc(TotalZones * 5 * sizeof(double));
  double *U0 = (double*) malloc(TotalZones * 5 * sizeof(double));

  for (int m=0; m < 5 * TotalZones; ++m) { U0[m] = 0.0; }
  for (int n=0; n<TotalZones; ++n) {
    fluids_state_derive(fluid[n], &U0[5*n], FLUIDS_CONSERVED);
  }

  timederiv(L);
  for (int n=0; n<TotalZones; ++n) {
    double U[5];
    fluids_state_derive(fluid[n], U, FLUIDS_CONSERVED);
    for (int q=0; q<5; ++q) {
      U[q] += L[5*n + q] * dt * 0.5;
    }
    fluids_state_fromcons(fluid[n], U, FLUIDS_CACHE_DEFAULT);
  }

  timederiv(L);
  for (int n=0; n<TotalZones; ++n) {
    double U[5];
    fluids_state_derive(fluid[n], U, FLUIDS_CONSERVED);
    for (int q=0; q<5; ++q) {
      U[q] = U0[5*n + q] + L[5*n + q] * dt;
    }
    fluids_state_fromcons(fluid[n], U, FLUIDS_CACHE_DEFAULT);
  }
  free(L);
  free(U0);
}

double fish_grav1d_maxwavespeed()
{
  double a = 0.0;
  for (int n=0; n<TotalZones; ++n) {
    double A[5];
    fluids_state_derive(fluid[n], A, FLUIDS_EVAL0);
    for (int q=0; q<5; ++q) {
      a = fabs(A[q]) > a ? fabs(A[q]) : a;
    }
  }
  return a;
}

void fish_grav1d_getprim(double *prim, double *grav)
{
  for (int n=0; n<TotalZones; ++n) {
    double P[5];
    double G[4];
    fluids_state_getattr(fluid[n], P, FLUIDS_PRIMITIVE);
    fluids_state_getattr(fluid[n], G, FLUIDS_GRAVITY);
    memcpy(&prim[5*n], P, 5 * sizeof(double));
    memcpy(&grav[4*n], G, 4 * sizeof(double));
  }
}

void fish_grav1d_setprim(double *prim)
{
  for (int n=0; n<TotalZones; ++n) {
    fluids_state_setattr(fluid[n], &prim[5*n], FLUIDS_PRIMITIVE);
  }
}

void fish_grav1d_mapbuffer(double *x, long flag)
{
  int nq = fluids_descr_getncomp(descr, flag);
  for (int n=0; n<TotalZones; ++n) {
    fluids_state_mapbuffer(fluid[n], &x[nq*n], flag);
  }
}

void solve_poisson(double *Rho, double *Phi, double *Gph, double *rhobar)
{
  int Ng = NumGhostZones;
  int N = TotalZones - Ng;

  fftw_complex *Rhox = fftw_alloc_complex(N);
  fftw_complex *Phix = fftw_alloc_complex(N);
  fftw_complex *Gphx = fftw_alloc_complex(N);
  fftw_complex *Rhok = fftw_alloc_complex(N);
  fftw_complex *Phik = fftw_alloc_complex(N);
  fftw_complex *Gphk = fftw_alloc_complex(N);

  fftw_plan fwd = fftw_plan_dft_1d(N, Rhox, Rhok, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan revphi = fftw_plan_dft_1d(N, Phik, Phix, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan revgph = fftw_plan_dft_1d(N, Gphk, Gphx, FFTW_BACKWARD, FFTW_ESTIMATE);

  for (int i=0; i<N; ++i) {
    Rhox[i][0] = Rho[i+Ng];
    Rhox[i][1] = 0.0;
  }
  fftw_execute(fwd);

  *rhobar = Rhok[0][0] / N;
  Phik[0][0] = 0.0;
  Phik[0][1] = 0.0;
  Gphk[0][0] = 0.0;
  Gphk[0][1] = 0.0;

  /*
   * gph[k] = I * k * (phi[k].re + I * phi[k].im)
   */
  for (int i=1; i<N; ++i) {
    double k = 2*M_PI * (i < N/2 ? i : i-N);
    Phik[i][0] = -Rhok[i][0] / (k*k);
    Phik[i][1] = -Rhok[i][1] / (k*k);
    Gphk[i][0] = -Phik[i][1] * k;
    Gphk[i][1] = +Phik[i][0] * k;
  }
  fftw_execute(revphi);
  fftw_execute(revgph);

  for (int i=0; i<N; ++i) {
    Phi[i+Ng] = Phix[i][0] / N;
    Gph[i+Ng] = Gphx[i][0] / N;
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
}

void timederiv(double *L)
{
  int Ng = NumGhostZones;
  int N = TotalZones;

  fish_state *S = fish_new();
  fish_setparami(S, FISH_PLM, FISH_RECONSTRUCTION);
  fish_setparami(S, FLUIDS_RIEMANN_HLLC, FISH_RIEMANN_SOLVER);
  fish_setparami(S, FISH_GODUNOV, FISH_SOLVER_TYPE);
  fish_setparamd(S, 2.0, FISH_PLM_THETA);

  double *Rho = (double*) malloc(N * sizeof(double));
  double *Phi = (double*) malloc(N * sizeof(double));
  double *Gph = (double*) malloc(N * sizeof(double));
  double P[5];
  double rhobar;

  for (int i=0; i<N; ++i) {
    fluids_state_getattr(fluid[i], P, FLUIDS_PRIMITIVE);
    Rho[i] = P[0];
  }
  solve_poisson(Rho, Phi, Gph, &rhobar);

  for (int i=0; i<Ng; ++i) {
    Phi[    i    ] = Phi[ N - 2*Ng + i];
    Phi[N - i - 1] = Phi[-1 + 2*Ng - i];
    Gph[    i    ] = Gph[ N - 2*Ng + i];
    Gph[N - i - 1] = Gph[-1 + 2*Ng - i];
  }

  for (int i=0; i<N; ++i) {
    double G[4];
    G[0] = Phi[i];
    G[1] = Gph[i];
    G[2] = 0.0;
    G[3] = 0.0;
    fluids_state_setattr(fluid[i], G, FLUIDS_GRAVITY);
  }
  fluids_descr_setrhobar(descr, rhobar);
  free(Rho);
  free(Phi);
  free(Gph);

  for (int m=0; m < 5 * TotalZones; ++m) {
    L[m] = 0.0;
  }
  fish_timederivative(S, fluid, 1, &N, &dx, L);

  double source[5];
  for (int i=0; i<TotalZones; ++i) {
    fluids_state_derive(fluid[i], source, FLUIDS_SOURCETERMS);
    for (int q=0; q<5; ++q) {
      L[5*i + q] += source[q];
    }
  }

  for (int i=0; i<Ng; ++i) {
    for (int q=0; q<5; ++q) {
      L[(    i    )*5 + q] = L[( N - 2*Ng + i)*5 + q];
      L[(N - i - 1)*5 + q] = L[(-1 + 2*Ng - i)*5 + q];
    }
  }
  fish_del(S);
}
