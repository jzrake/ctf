
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define FISH_PRIVATE_DEFS
#include "fish.h"

#define MAXQ 8 // for small statically-declared arrays

static void _pcm(fish_state *S, fluids_state **src, fluids_state *L,
		 fluids_state *R, long flag);
static void _plm(fish_state *S, fluids_state **src, fluids_state *L,
		 fluids_state *R, long flag);
static void _weno5(fish_state *S, fluids_state **src, fluids_state *L,
		   fluids_state *R, long flag);
static int _matrix_product(double *A, double *B, double *C,
			   int ni, int nj, int nk);
static int _intercell_godunov(fish_state *S, fluids_state **fluid, double *Fiph,
			      int N, int dim);
static int _intercell_spectral(fish_state *S, fluids_state **fluid,
			       double *Fiph, int N, int dim);
static const long FLUIDS_FLUX[3] = {FLUIDS_FLUX0, FLUIDS_FLUX1, FLUIDS_FLUX2};
static const long FLUIDS_EVAL[3] = {FLUIDS_EVAL0, FLUIDS_EVAL1, FLUIDS_EVAL2};
static const long FLUIDS_LEVECS[3] = {FLUIDS_LEVECS0,
				      FLUIDS_LEVECS1,
				      FLUIDS_LEVECS2};
static const long FLUIDS_REVECS[3] = {FLUIDS_REVECS0,
				      FLUIDS_REVECS1,
				      FLUIDS_REVECS2};

fish_state *fish_new(void)
{
  fish_state *S = (fish_state*) malloc(sizeof(fish_state));
  fish_state state = {
    .solver_type = FISH_GODUNOV,
    .riemann_solver = FLUIDS_RIEMANN_HLL,
    .reconstruction = FISH_PLM,
    .smoothness_indicator = FISH_ISK_JIANGSHU96,
    .plm_theta = 2.0,
    .shenzha10_param = 0.0,
  } ;
  *S = state;
  return S;
}
int fish_del(fish_state *S)
{
  free(S);
  return 0;
}

int fish_getparami(fish_state *S, int *param, long flag)
{
  switch (flag) {
  case FISH_SOLVER_TYPE: *param = S->solver_type; return 0;
  case FISH_RIEMANN_SOLVER: *param = S->riemann_solver; return 0;
  case FISH_RECONSTRUCTION: *param = S->reconstruction; return 0;
  case FISH_SMOOTHNESS_INDICATOR: *param = S->smoothness_indicator; return 0;
  }
  return FISH_ERROR_BADARG;
}
int fish_setparami(fish_state *S, int param, long flag)
{
  switch (flag) {
  case FISH_SOLVER_TYPE: S->solver_type = param; return 0;
  case FISH_RIEMANN_SOLVER: S->riemann_solver = param; return 0;
  case FISH_RECONSTRUCTION: S->reconstruction = param; return 0;
  case FISH_SMOOTHNESS_INDICATOR: S->smoothness_indicator = param; return 0;
  }
  return FISH_ERROR_BADARG;
}
int fish_getparamd(fish_state *S, double *param, long flag)
{
  switch (flag) {
  case FISH_PLM_THETA: *param = S->plm_theta; return 0;
  case FISH_SHENZHA10_PARAM: *param = S->shenzha10_param; return 0;
  }
  return FISH_ERROR_BADARG;
}
int fish_setparamd(fish_state *S, double param, long flag)
{
  switch (flag) {
  case FISH_PLM_THETA: S->plm_theta = param; return 0;
  case FISH_SHENZHA10_PARAM: S->shenzha10_param = param; return 0;
  }
  return FISH_ERROR_BADARG;
}

int fish_intercellflux(fish_state *S, fluids_state **fluid, double *F, int N,
                       int dim)
{
  switch (S->solver_type) {
  case FISH_GODUNOV: return _intercell_godunov(S, fluid, F, N, dim);
  case FISH_SPECTRAL: return _intercell_spectral(S, fluid, F, N, dim);
  default: return FISH_ERROR_BADARG;
  }
}

int fish_timederivative(fish_state *S, fluids_state **fluid,
			int ndim, int *shape, double *dx,
			double *U, double *L)
/*
 * NOTE: L is assumed to already be initialized to zero's.
*/
{
  int si, sj, sk, Q, numerr = 0;
  double *Fiph;
  fluids_state **slice;
  fluids_descr *D;
  fluids_state_getdescr(fluid[0], &D);
  Q = fluids_descr_getncomp(D, FLUIDS_PRIMITIVE);

  switch (ndim) {
    /* ---------------------------------- (1d) -------------------------------*/
  case 1:
    for (int n=0; n<shape[0]; ++n) {
      int e = fluids_state_fromcons(fluid[n], &U[Q*n], FLUIDS_CACHE_DEFAULT);
      numerr += (e != 0);
    }
    if (numerr == 0) {
      return numerr;
    }
    Fiph = (double*) malloc(shape[0] * Q * sizeof(double));
    fish_intercellflux(S, fluid, Fiph, shape[0], 0);
    for (int n=Q; n<shape[0]*Q; ++n) {
      L[n] -= (Fiph[n] - Fiph[n-Q]) / dx[0];
    }
    free(Fiph);
    break;
    /* ---------------------------------- (2d) -------------------------------*/
  case 2:
    si = shape[1];
    sj = 1;
    for (int n=0; n<shape[0]*shape[1]; ++n) {
      int e = fluids_state_fromcons(fluid[n], &U[Q*n], FLUIDS_CACHE_DEFAULT);
      numerr += (e != 0);
    }

    // ----------------------------
    // sweeps along the x-direction
    // ----------------------------
    slice = (fluids_state **) malloc(shape[0] * sizeof(fluids_state*));
    Fiph = (double*) malloc(shape[0] * Q *sizeof(double));
    for (int j=0; j<shape[1]; ++j) {
      for (int i=0; i<shape[0]; ++i) {
	slice[i] = fluid[i*si + j*sj];
      }
      fish_intercellflux(S, slice, Fiph, shape[0], 0);
      for (int i=0; i<shape[0]; ++i) {
	for (int q=0; q<Q; ++q) {
	  L[(i*si + j*sj)*Q + q] -= (Fiph[i*Q+q] - Fiph[(i-1)*Q+q]) / dx[0];
	}
      }
    }
    free(Fiph);
    free(slice);

    // ----------------------------
    // sweeps along the y-direction
    // ----------------------------
    slice = (fluids_state **) malloc(shape[1] * sizeof(fluids_state*));
    Fiph = (double*) malloc(shape[1] * Q *sizeof(double));
    for (int i=0; i<shape[0]; ++i) {
      for (int j=0; j<shape[1]; ++j) {
	slice[j] = fluid[i*si + j*sj];
      }
      fish_intercellflux(S, slice, Fiph, shape[1], 1);
      for (int j=0; j<shape[1]; ++j) {
	for (int q=0; q<Q; ++q) {
	  L[(i*si + j*sj)*Q + q] -= (Fiph[j*Q+q] - Fiph[(j-1)*Q+q]) / dx[1];
	}
      }
    }
    free(Fiph);
    free(slice);
    break;

  case 3:
    /* ---------------------------------- (3d) -------------------------------*/
    si = shape[2] * shape[1];
    sj = shape[2];
    sk = 1;
    break;
  }

  return 0;
}


int _intercell_godunov(fish_state *S, fluids_state **fluid, double *F, int N,
                       int dim)
{
  fluids_descr *D;
  fluids_state_getdescr(fluid[0], &D);

  int Q = fluids_descr_getncomp(D, FLUIDS_PRIMITIVE);
  fluids_state *S_ = fluids_state_new();
  fluids_state *SL = fluids_state_new();
  fluids_state *SR = fluids_state_new();

  /* prevent the use of uninitialized bytes */
  for (int n=0; n<N*Q; ++n) {
    F[n] = 0.0;
  }

  /* Assumes all states have the same descriptor, after all what sense does this
     make otherwise? */
  fluids_state_setdescr(S_, D);
  fluids_state_setdescr(SL, D);
  fluids_state_setdescr(SR, D);
  fluids_state_cache(S_, FLUIDS_CACHE_CREATE);
  fluids_state_cache(SL, FLUIDS_CACHE_CREATE);
  fluids_state_cache(SR, FLUIDS_CACHE_CREATE);

  fluids_riemn *R = fluids_riemn_new();
  fluids_riemn_setsolver(R, S->riemann_solver);
  fluids_riemn_setdim(R, dim);
  fluids_riemn_setstateL(R, SL);
  fluids_riemn_setstateR(R, SR);

  switch (S->reconstruction) {
  case FISH_PCM:
    for (int n=0; n<N-1; ++n) {
      _pcm(S, &fluid[n], SL, SR, FLUIDS_PRIMITIVE);
      _pcm(S, &fluid[n], SL, SR, FLUIDS_GRAVITY);
      fluids_riemn_execute(R);
      fluids_riemn_sample(R, S_, 0.0);
      fluids_state_derive(S_, &F[Q*n], FLUIDS_FLUX[dim]);
    }
    break;
  case FISH_PLM:
    for (int n=1; n<N-2; ++n) {
      _plm(S, &fluid[n], SL, SR, FLUIDS_PRIMITIVE);
      _plm(S, &fluid[n], SL, SR, FLUIDS_GRAVITY);
      fluids_riemn_execute(R);
      fluids_riemn_sample(R, S_, 0.0);
      fluids_state_derive(S_, &F[Q*n], FLUIDS_FLUX[dim]);
    }
    break;
  case FISH_WENO5:
    for (int n=2; n<N-3; ++n) {
      _weno5(S, &fluid[n], SL, SR, FLUIDS_PRIMITIVE);
      _weno5(S, &fluid[n], SL, SR, FLUIDS_GRAVITY);
      fluids_riemn_execute(R);
      fluids_riemn_sample(R, S_, 0.0);
      fluids_state_derive(S_, &F[Q*n], FLUIDS_FLUX[dim]);
    }
  default:
    break;
  }

  fluids_state_del(SL);
  fluids_state_del(SR);
  fluids_state_del(S_);
  fluids_riemn_del(R);
  return 0;
}


int _intercell_spectral(fish_state *S, fluids_state **fluid, double *Fiph,
			int N, int dim)
/* -----------------------------------------------------------------------------
 *
 * This function uses characteristic decomposition to find the intercell
 * fluxes. The fluid system given must provide the full set of left and right
 * eigenvectors. PCM, PLM, and WENO5 reconstruction may be used. There's
 * presently no support for fluid systems woth gravity, although adding it here
 * is only a matter of choosing a reconstruction type for the gravitational
 * field. Right now the function hard-codes Q=5 primitive variables, but that
 * can easily be changed.
 *
 *
 * (1) Compute max eigenvalue A, flux F, and conserved U in each zone
 *
 * (2) Arithmetically average L/R primitive states to form L/R face-centered
 *     eigenvectors
 *
 * (3) Get max wave-speed over local 6-zone stencil
 *
 * (4) Create F+ and F- fluxes over each of those 6 zones
 *
 * (5) Decompose F+ and F- into f+ and f-, characteristic right and left-going
 *     fluxes in each zone
 *
 * (6) Use reconstruction on the stencil to get a left and right going
 *     characteristic f
 *
 * (7) Rotate f into back into F
 *
 * -----------------------------------------------------------------------------
 */
{
  double Pl[5], Pr[5], Pface[5];
  double Liph[5][5], Riph[5][5];
  double fp[6][5];
  double fm[6][5];
  double fpT[5][6];
  double fmT[5][6];
  double Fp[5];
  double Fm[5];
  double f[5];

  fluids_state *face = fluids_state_new();
  fluids_descr *D;
  fluids_state_getdescr(fluid[0], &D);
  fluids_state_setdescr(face, D);
  fluids_state_cache(face, FLUIDS_CACHE_CREATE);

  int Q = fluids_descr_getncomp(D, FLUIDS_PRIMITIVE);
  double *A = (double*) malloc(N*1*sizeof(double)); // array of max eigenvalues
  double *F = (double*) malloc(N*Q*sizeof(double)); // array of fluxes
  double *U = (double*) malloc(N*Q*sizeof(double)); // array of conserved

  /* prevent the use of uninitialized bytes */
  for (int n=0; n<N*Q; ++n) {
    F[n] = 0.0;
  }

  for (int n=0; n<N; ++n) {
    /*--------------------------------- (1) --------------------------------- */
    double lam[5];
    long flags = FLUIDS_EVAL[dim] | FLUIDS_FLUX[dim] | FLUIDS_CONSERVED;

    fluids_state_derive(fluid[n], NULL, flags);
    fluids_state_getcached(fluid[n], lam, FLUIDS_EVAL[dim]);
    fluids_state_getcached(fluid[n], &F[Q*n], FLUIDS_FLUX[dim]);
    fluids_state_getcached(fluid[n], &U[Q*n], FLUIDS_CONSERVED);

    A[n] = 0.0;
    for (int q=0; q<Q; ++q) {
      if (fabs(lam[q]) > A[n]) {
       	A[n] = fabs(lam[q]);
      }
    }
  }

  for (int n=2; n<N-3; ++n) {

    /*--------------------------------- (2) --------------------------------- */
    fluids_state_getattr(fluid[n+0], Pl, FLUIDS_PRIMITIVE);
    fluids_state_getattr(fluid[n+1], Pr, FLUIDS_PRIMITIVE);
    for (int q=0; q<Q; ++q) {
      Pface[q] = 0.5*(Pl[q] + Pr[q]);
    }
    fluids_state_setattr(face, Pface, FLUIDS_PRIMITIVE);
    fluids_state_derive(face, NULL, FLUIDS_LEVECS[dim] | FLUIDS_REVECS[dim]);
    fluids_state_getcached(face, Liph[0], FLUIDS_LEVECS[dim]);
    fluids_state_getcached(face, Riph[0], FLUIDS_REVECS[dim]);

    /*--------------------------------- (3) --------------------------------- */
    double ml = 0.0;
    for (int j=0; j<6; ++j) {
      if (fabs(A[n+j-2]) > ml) {
	ml = fabs(A[n+j-2]);
      }
    }

    /*--------------------------------- (4,5) ------------------------------- */
    for (int j=0; j<6; ++j) {
      for (int q=0; q<Q; ++q) {
	/*
	 * local Lax-Friedrichs flux splitting
	 */
	int m = (n+j-2)*Q + q;
	Fp[q] = 0.5*(F[m] + ml*U[m]);
	Fm[q] = 0.5*(F[m] - ml*U[m]);
      }
      _matrix_product(Liph[0], Fp, fp[j], Q, 1, Q);
      _matrix_product(Liph[0], Fm, fm[j], Q, 1, Q);
    }
    for (int q=0; q<Q; ++q) {
      for (int j=0; j<6; ++j) {
	fpT[q][j] = fp[j][q];
	fmT[q][j] = fm[j][q];
      }
    }

    /*--------------------------------- (6,7) ------------------------------- */
    switch (S->reconstruction) {
    case FISH_PCM:
      for (int q=0; q<Q; ++q) {
	f[q] = (_reconstruct(S, &fpT[q][2], PCM_C2R) +
		_reconstruct(S, &fmT[q][3], PCM_C2L));
      }
      break;
    case FISH_PLM:
      for (int q=0; q<Q; ++q) {
	f[q] = (_reconstruct(S, &fpT[q][2], PLM_C2R) +
		_reconstruct(S, &fmT[q][3], PLM_C2L));
      }
      break;
    case FISH_WENO5:
      for (int q=0; q<Q; ++q) {
	f[q] = (_reconstruct(S, &fpT[q][2], WENO5_FD_C2R) +
		_reconstruct(S, &fmT[q][3], WENO5_FD_C2L));
      }
      break;
    }
    _matrix_product(Riph[0], f, &Fiph[n*Q], Q, 1, Q);
    double P_[5];
    fluids_state_getattr(face, P_, FLUIDS_PRIMITIVE);
  }
  fluids_state_del(face);
  free(A);
  free(F);
  free(U);
  return 0;
}

void _pcm(fish_state *S, fluids_state **src, fluids_state *L,
	  fluids_state *R, long flag)
{
  double Pl[MAXQ], Pr[MAXQ];
  fluids_descr *D;
  int fluid;
  fluids_state_getdescr(src[0], &D);
  fluids_descr_getfluid(D, &fluid);
  fluids_state_getattr(src[0], Pl, flag);
  fluids_state_getattr(src[1], Pr, flag);
  fluids_state_setattr(L, Pl, flag);
  fluids_state_setattr(R, Pr, flag);
}

void _plm(fish_state *S, fluids_state **src, fluids_state *L,
	  fluids_state *R, long flag)
{
  double P[4][MAXQ];
  double Pl[MAXQ], Pr[MAXQ];
  fluids_descr *D;
  int fluid;
  fluids_state_getdescr(src[0], &D);
  fluids_descr_getfluid(D, &fluid);
  int Q = fluids_descr_getncomp(D, flag);
  for (int j=-1; j<3; ++j) {
    fluids_state_getattr(src[j], P[j+1], flag);
  }
  for (int q=0; q<Q; ++q) {
    double v[4];
    for (int j=0; j<4; ++j) {
      v[j] = P[j][q];
    }
    Pl[q] = _reconstruct(S, &v[1], PLM_C2R);
    Pr[q] = _reconstruct(S, &v[2], PLM_C2L);
  }
  fluids_state_setattr(L, Pl, flag);
  fluids_state_setattr(R, Pr, flag);
}

void _weno5(fish_state *S, fluids_state **src, fluids_state *L,
	    fluids_state *R, long flag)
{
  double P[6][MAXQ];
  double Pl[MAXQ], Pr[MAXQ];
  fluids_descr *D;
  int fluid;
  fluids_state_getdescr(src[0], &D);
  fluids_descr_getfluid(D, &fluid);
  int Q = fluids_descr_getncomp(D, flag);
  for (int j=-2; j<4; ++j) {
    fluids_state_getattr(src[j], P[j+2], flag);
  }
  for (int q=0; q<Q; ++q) {
    double v[6];
    for (int j=0; j<6; ++j) {
      v[j] = P[j][q];
    }
    Pl[q] = _reconstruct(S, &v[2], WENO5_FD_C2R);
    Pr[q] = _reconstruct(S, &v[3], WENO5_FD_C2L);
  }
  fluids_state_setattr(L, Pl, flag);
  fluids_state_setattr(R, Pr, flag);
}


int _matrix_product(double *A, double *B, double *C, int ni, int nj, int nk)
{
  for (int i=0; i<ni; ++i) {
    for (int j=0; j<nj; ++j) {
      C[i*nj+j] = 0.0;
      for (int k=0; k<nk; ++k) {
	C[i*nj+j] += A[i*nk+k] * B[k*nj+j];
      }
    }
  }
  return 0;
}
