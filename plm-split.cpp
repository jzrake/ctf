
#include "plm-split.hpp"
#include "weno.h"
#include "logging.hpp"

#define MAXNQ 8 // Used for static array initialization
typedef MethodOfLinesSplit Deriv;

std::valarray<double> Deriv::dUdt(const std::valarray<double> &Uin)
{
  this->prepare_integration();

  std::valarray<double> U = Uin;
  std::valarray<double> L(U.size());
  std::valarray<double> &P = Mara->PrimitiveArray;

  ConsToPrim(U, P);
  DriveSweeps(P, L);

  return L;
}
void Deriv::DriveSweeps(const std::valarray<double> &P,
                        std::valarray<double> &L)
{
   switch (ND) {
   case 1: drive_sweeps_1d(&P[0], &L[0]); break;
   case 2: drive_sweeps_2d(&P[0], &L[0]); break;
   case 3: drive_sweeps_3d(&P[0], &L[0]); break;
   }
}
void Deriv::intercell_flux_sweep(const double *P, double *F, int dim)
{
  int i,S=stride[dim];
  double Pl[MAXNQ], Pr[MAXNQ];

  for (i=2*S; i<stride[0]-3*S; i+=NQ) {
    for (int q=0; q<NQ; ++q) {
      const int m = i + q;
      double v[6] = { P[m-2*S], P[m-S], P[m+0], P[m+1*S], P[m+2*S], P[m+3*S] };

      switch (GodunovOperator::reconstruct_method) {
      case RECONSTRUCT_PCM:
	Pl[q] = v[2];
	Pr[q] = v[3];
	break;
      case RECONSTRUCT_PLM:
	Pl[q] = reconstruct(&v[2], PLM_C2R);
	Pr[q] = reconstruct(&v[3], PLM_C2L);
	break;
      case RECONSTRUCT_WENO5:
	Pl[q] = reconstruct(&v[2], WENO5_FV_C2R);
	Pr[q] = reconstruct(&v[3], WENO5_FV_C2L);
	break;
      }
    }
    int error = Mara->riemann->IntercellFlux(Pl, Pr, 0, &F[i], 0.0, dim);
    if (error) {
      error = Mara->riemann->IntercellFlux(&P[i], &P[i+S], 0, &F[i], 0.0, dim);

      DebugLog.Warning(__FUNCTION__) << "Reverting to first order at zone interface... ";
      if (!error) DebugLog.Warning() << "Success!" << std::endl;
      else        DebugLog.Warning() << "Still failed!" << std::endl;
    }
  }
}

void Deriv::drive_sweeps_1d(const double *P, double *L)
{
  double *F = (double*) malloc(stride[0]*sizeof(double));

  int i,sx=stride[1];
  intercell_flux_sweep(P,F,1);

  for (i=sx; i<stride[0]; ++i) {
    L[i] = -(F[i]-F[i-sx])/dx;
  }
  free(F);
}
void Deriv::drive_sweeps_2d(const double *P, double *L)
{
  double *F = (double*) malloc(stride[0]*sizeof(double));
  double *G = (double*) malloc(stride[0]*sizeof(double));

  int i,sx=stride[1],sy=stride[2];
  intercell_flux_sweep(P,F,1);
  intercell_flux_sweep(P,G,2);

  Mara->fluid->ConstrainedTransport2d(F,G,stride);

  for (i=sx; i<stride[0]; ++i) {
    L[i] = -(F[i]-F[i-sx])/dx - (G[i]-G[i-sy])/dy;
  }
  free(F); free(G);
}
void Deriv::drive_sweeps_3d(const double *P, double *L)
{
  double *F = (double*) malloc(stride[0]*sizeof(double));
  double *G = (double*) malloc(stride[0]*sizeof(double));
  double *H = (double*) malloc(stride[0]*sizeof(double));

  int i,sx=stride[1],sy=stride[2],sz=stride[3];
  intercell_flux_sweep(P,F,1);
  intercell_flux_sweep(P,G,2);
  intercell_flux_sweep(P,H,3);

  Mara->fluid->ConstrainedTransport3d(F,G,H,stride);

  for (i=sx; i<stride[0]; ++i) {
    L[i] = -(F[i]-F[i-sx])/dx - (G[i]-G[i-sy])/dy - (H[i]-H[i-sz])/dz;
  }
  free(F); free(G); free(H);
}
