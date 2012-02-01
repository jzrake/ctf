
#include "plm-split.hpp"
#include "logging.hpp"

#define MAXNQ 8 // Used for static array initialization

typedef PlmMethodOfLinesSplit Deriv;
static double plm_theta = 2.0;

/*------------------------------------------------------------------------------
 *
 * Private inline functions
 *
 */
static inline double sign(double x)
{
  return (x>0)-(x<0);
}
static inline double min3(double a, double b, double c)
{
  const double ab=(a<b)?a:b;
  return (ab<c)?ab:c;
}
static inline double plm_minmod(double ul, double u0, double ur)
{
  const double a = plm_theta * (u0 - ul);
  const double b =     0.5   * (ur - ul);
  const double c = plm_theta * (ur - u0);
  return 0.25*fabs(sign(a)+sign(b))*(sign(a)+sign(c))*min3(fabs(a),fabs(b),fabs(c));
}


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
void Deriv::reconstruct_plm(const double *P0, double *Pl, double *Pr, int S)
{
  int i,T=2*S;
  for (i=0; i<NQ; ++i) {
    Pr[i] = P0[S+i] - 0.5 * plm_minmod(P0[ 0+i], P0[S+i], P0[T+i]);
    Pl[i] = P0[0+i] + 0.5 * plm_minmod(P0[-S+i], P0[0+i], P0[S+i]);
  }
}
void Deriv::intercell_flux_sweep(const double *P, double *F, int dim)
{
  int i,S=stride[dim];
  double Pl[MAXNQ], Pr[MAXNQ];

  for (i=S; i<stride[0]-2*S; i+=NQ) {
    reconstruct_plm(&P[i], Pl, Pr, S);
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
