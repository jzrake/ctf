
/*------------------------------------------------------------------------------
 * FILE: rmhd-c2p-eos.c
 *
 * AUTHOR: Paul Duffell and Jonathan Zrake, NYU CCPP
 *
 * REFERENCES:
 *
 * Anton, L., Zanotti, O., Miralles, J. A., Martı, J. M., Ibanez, J. M., Font,
 * J. A., & Pons, J. A. 2006, The Astrophysical Journal, 637, 296
 *
 * Noble, S. C., Gammie, C. F., McKinney, J. C., & Zanna, L. D.  2006, The
 * Astrophysical Journal, 641, 626
 *
 *
 * DESCRIPTION:
 *
 * CONVENTIONS:
 *
 * -----------------------------------------------------------------------------
 * W     : rho*h*LorentzFactor^2
 * K     : v^2 LorentzFactor^2
 * T     : temperature
 * qdotn : -(Tau + D) (negative of total energy)
 * qdotb : momentum dotted with magnetic field
 * Q2    : momentum squared
 * B2    : magnetic field squared
 * D     : rho*LorentzFactor
 * ------------------------------------------------------------------------------
 */


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include "hydro.hpp"
#include "rmhd-c2p.h"

enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
enum { rho, pre, vx, vy, vz };             // Primitive

static int WsolveWKT(double *W , double *gamma, double *T,
                     double qdotn, double qdotb, double Q2, double B2, double D);
static void fghWKT(double W, double K, double T,
                   double *deltaW, double *deltaK, double *deltaT,
                   double qdotn, double qdotb,
                   double Q2, double B2, double D,
                   double *ferr, double *gerr, double *herr);


static const int MaxIteration = 250;
static const double Tolerance = 1e-12;
static const double bigZ = +1.0e20;
static const double smlZ =  0.0;
static const double bigT = +10.0;
static const double smlT = -10.0;

static int Iterations = 0;


#define NR_ERROR (fabs(deltax/x) + deltay*deltay/y/y + deltaz*deltaz/z/z)
#define YMAX 1e14

using std::isnan;
using std::isinf;

static const int verbose = 0;
static const double NR_TOL = 1e-14;
static const EquationOfState *eos = NULL;

static double D,Tau;
static double S2,B2,BS,BS2;
static double Cons[8];
static double Z_start, W_start, T_start;


void rmhd_c2p_eos_set_eos(const void *eos_)
{
  eos = static_cast<const EquationOfState*>(eos_);
}
void rmhd_c2p_eos_new_state(const double *U)
{
  D    = U[ddd];
  Tau  = U[tau];
  S2   = U[Sx]*U[Sx] + U[Sy]*U[Sy] + U[Sz]*U[Sz];
  B2   = U[Bx]*U[Bx] + U[By]*U[By] + U[Bz]*U[Bz];
  BS   = U[Bx]*U[Sx] + U[By]*U[Sy] + U[Bz]*U[Sz];
  BS2  = BS*BS;

  memcpy(Cons, U, 8*sizeof(double));
}
int rmhd_c2p_eos_get_iterations()
{
  return Iterations;
}

// This estimate becomes exact for no magnetic field in the NR limit.
// -----------------------------------------------------------------------------
void rmhd_c2p_eos_estimate_from_cons()
{
  W_start = 1.0;
  Z_start = sqrt(S2 + D*D);
  T_start = eos->Temperature_u(D, Tau - 0.5*B2);
}
void rmhd_c2p_eos_set_starting_prim(const double *P)
{
  const double V2 = P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz];
  const double W2 = 1.0 / (1.0 - V2);

  W_start = sqrt(W2);
  T_start = eos->Temperature_p(P[rho], P[pre]);
  Z_start = (P[rho] + eos->Internal(P[rho], T_start) + P[pre])*W2;
}




// Solution based on Noble et. al. (2006), using Z = rho h W^2 and T as the
// unkowns.
// -----------------------------------------------------------------------------
int rmhd_c2p_eos_solve_noble2dzt(double *Pout)
{
  int bad_input = rmhd_c2p_check_cons(Cons);
  if (bad_input) {
    return bad_input;
  }

  // Starting values
  // ---------------------------------------------------------------------------
  Iterations = 0;

  double error = 1.0;
  double Z = Z_start;
  double T = T_start;

  while (error > Tolerance) {
    double Jp[2], Ju[2];

    const double Z2  = Z*Z;
    const double Z3  = Z*Z2;
    const double a   = S2*Z2 + BS2*(B2 + 2*Z);
    const double b   = (B2 + Z)*(B2 + Z)*Z2;
    const double ap  = 2*(S2*Z + BS2);          // da/dZ
    const double bp  = 2*Z*(B2 + Z)*(B2 + 2*Z); // db/dZ
    const double V2  = a / b;
    const double W2  = 1.0 / (1.0 - V2);
    const double W   = sqrt(W2);
    const double W3  = W*W2;
    const double W4  = W*W3;

    if (V2 > 1.0) {
      //      printf("internal superluminal: %e %e %e %e %e %e\n", Z, T, W2, V2, a, b);
    }

    const double p = eos->Derivatives_p(D/W, T, Jp);
    const double u = eos->Derivatives_u(D/W, T, Ju);

    const double dpdRho = Jp[0];
    const double dudRho = Ju[0];

    const double dpdT = Jp[1];
    const double dudT = Ju[1];

    const double dv2dZ   = (ap*b - bp*a) / (b*b); // (a'b - b'a) / b^2
    const double dW2dZ   =   W4 * dv2dZ;
    const double dRhodZ  = -0.5 * (D/W3) * dW2dZ;
    const double dpdZ    = dpdRho * dRhodZ;
    const double dudZ    = dudRho * dRhodZ;
    const double dRhohdZ = dRhodZ + dudZ + dpdZ;
    const double dRhohdT = (dudT + dpdT)*W2;

    const double df0dZ = -0.5*B2*dv2dZ - BS2/Z3 - 1.0 + dpdZ;
    const double df0dT = dpdT;

    const double df1dZ = 1.0 - (dRhohdZ*W2 + (D/W + u + p)*dW2dZ);
    const double df1dT = -dRhohdT;

    double f[2];
    double J[4], G[4];

    f[0] = Tau + D - 0.5*B2*(1+V2) + 0.5*BS2/Z2 - Z + p; // equation (29)
    f[1] = Z - (D/W + u + p)*W2;

    J[0] = df0dZ; J[1] = df0dT;
    J[2] = df1dZ; J[3] = df1dT;

    // G in the inverse Jacobian
    // -------------------------------------------------------------------------
    const double det = J[0]*J[3] - J[1]*J[2];
    G[0] =  J[3]/det; G[1] = -J[1]/det;
    G[2] = -J[2]/det; G[3] =  J[0]/det;      // G = J^{-1}

    const double dZ = -(G[0]*f[0] + G[1]*f[1]); // Matrix multiply, dx = -G . f
    const double dT = -(G[2]*f[0] + G[3]*f[1]);

    // Bracketing the root
    // -------------------------------------------------------------------------
    double Z_new = Z + dZ;
    double T_new = T + dT;

    Z_new = (Z_new > smlZ) ? Z_new : -Z_new;
    Z_new = (Z_new < bigZ) ? Z_new :  Z;

    T_new = (T_new > smlT) ? T_new : smlT;
    T_new = (T_new < bigT) ? T_new : bigT;

    Z = Z_new;
    T = T_new;

    // -------------------------------------------------------------------------
    error = fabs(dZ/Z) + fabs(dT/T);
    ++Iterations;

    if (Iterations == MaxIteration) {
      return RMHD_C2P_MAXITER;
    }
  }

  // Recover the W value from the converged Z value.
  // -------------------------------------------------------------------------
  const double Z2  = Z*Z;
  const double a   = S2*Z2 + BS2*(B2 + 2*Z);
  const double b   = (B2 + Z)*(B2 + Z)*Z2;
  const double V2  = a / b;
  const double W2  = 1.0 / (1.0 - V2);
  const double W   = sqrt(W2);
  const double Rho = D/W;
  const double b0  = BS * W / Z;

  double P[8];
  const double *U = Cons;

  P[rho] =  Rho;
  P[pre] =  eos->Pressure(Rho, T);
  P[vx ] = (U[Sx] + b0*U[Bx]/W) / (Z+B2);
  P[vy ] = (U[Sy] + b0*U[By]/W) / (Z+B2);
  P[vz ] = (U[Sz] + b0*U[Bz]/W) / (Z+B2);
  P[Bx ] =  U[Bx];
  P[By ] =  U[By];
  P[Bz ] =  U[Bz];

  int prim_error = rmhd_c2p_check_prim(P);
  if (prim_error) return prim_error;

  // Don't actually modify the user's prim state until we're sure the state is
  // good.
  // ---------------------------------------------------------------------------
  memcpy(Pout, P, 8*sizeof(double));
  return RMHD_C2P_SUCCESS;
}




int rmhd_c2p_eos_solve_duffell3d(double *Pout)
{
  double Z = Z_start;
  double T = T_start;
  double W = W_start;

  if (WsolveWKT(&Z, &W, &T, -(D+Tau), BS, S2, B2, D)) {
    return RMHD_C2P_MAXITER;
  }

  const double Rho = D/W;
  const double b0 = BS * W / Z;

  double P[8];
  const double *U = Cons;

  P[rho] =  Rho;
  P[pre] =  eos->Pressure(Rho, T);
  P[vx ] = (U[Sx] + b0*U[Bx]/W) / (Z+B2);
  P[vy ] = (U[Sy] + b0*U[By]/W) / (Z+B2);
  P[vz ] = (U[Sz] + b0*U[Bz]/W) / (Z+B2);
  P[Bx ] =  U[Bx];
  P[By ] =  U[By];
  P[Bz ] =  U[Bz];

  int prim_error = rmhd_c2p_check_prim(P);
  if (prim_error) return prim_error;

  // Don't actually modify the user's prim state until we're sure the state is
  // good.
  // ---------------------------------------------------------------------------
  memcpy(Pout, P, 8*sizeof(double));
  return RMHD_C2P_SUCCESS;
}




void fghWKT(double W, double K, double T,
            double *deltaW, double *deltaK, double *deltaT,
            double qdotn, double qdotb,
            double Q2, double B2, double D,
            double *ferr, double *gerr, double *herr)
{
  double gamm = sqrt(fabs(K+1.0));
  double v2 = K/(K+1.0);
  double Rho = D/gamm;
  double Ju[2], Jp[2];

  double u = eos->Derivatives_u(Rho, T, Ju);
  double P = eos->Derivatives_p(Rho, T, Jp);

  double dudrho = Ju[0];
  double dudT   = Ju[1];
  double dPdrho = Jp[0];
  double dPdT   = Jp[1];

  double W2 = W*W;
  double X  = B2 + W;

  double f,g,h;
  double dfdW,dgdW,dhdW;
  double dfdK,dgdK,dhdK;
  double dfdT,dgdT,dhdT;

  f = W - P + qdotn - .5*qdotb*qdotb/W2 + .5*B2*(2.0*K+1.0)/(K+1.0);
  g = (X*X*W2*v2 - Q2*W2 - qdotb*qdotb*(B2+2.0*W));
  h = W - gamm*D - (K+1.0)*u - (K+1.0)*P;

  dfdW  =  1.0 + qdotb*qdotb/W/W2;
  dgdW  =  2.0 * (W*X*(B2+2.0*W)*v2 - Q2*W - qdotb*qdotb);
  dhdW  =  1.0;

  dfdK  =  0.5*B2/pow(K+1.0,2.0) + .5*D*dPdrho/pow(gamm,3.0);
  dgdK  =  W2*X*X/pow(K+1.0,2.0);
  dhdK  = -u - P + .5*D*(dudrho + dPdrho - 1.0)/gamm;

  dfdT  = -dPdT;
  dgdT  =  0.0;
  dhdT  = -(K+1.0)*(dudT + dPdT);

  double detJ = dfdT*(dgdW*dhdK - dgdK) + dhdT*(dfdW*dgdK - dfdK*dgdW);

  *deltaW = - ( f*(dgdK*dhdT       ) - g*(dfdK*dhdT - dfdT*dhdK) + h*(          - dfdT*dgdK) )/detJ;
  *deltaK = - (-f*(dgdW*dhdT       ) + g*(dfdW*dhdT - dfdT     ) - h*(          - dfdT*dgdW) )/detJ;
  *deltaT = - ( f*(dgdW*dhdK - dgdK) - g*(dfdW*dhdK - dfdK     ) + h*(dfdW*dgdK - dfdK*dgdW) )/detJ;

  if (verbose) {
    printf("Rho=%e T=%e: got derivatives dudD=%e dudT=%e dpdD=%e dpdT=%e ",
           Rho, T, Ju[0], Ju[1], Jp[0], Jp[1]);
    printf("dfdT=%e dhdT=%e detJ=%e\n", dfdT, dhdT, detJ);
  }

  *ferr = f/(W + .5*B2);
  *gerr = g/(X*X*W2*v2);
  *herr = h/W;
}

int WsolveWKT(double *W, double *gamma, double *T,
              double qdotn, double qdotb, double Q2, double B2, double D)
{
  Iterations = 0;

  int fail   = 0;
  int count  = 0;
  int Nmax   = 30;
  int Nmax2  = 30;
  int Nmaxsafe = 100;
  int Nextra = 5;

  double x = *W;
  double y = pow(*gamma,2.0)-1.0;
  double z = *T;
  double Told = z;
  double err1,err2,err3;

  double error = 1.0;
  double deltax,deltay,deltaz;

  while (error>NR_TOL && fail==0) {

    fghWKT(x,y,z,&deltax,&deltay,&deltaz,qdotn,qdotb,Q2,B2,D,&err1,&err2,&err3);
    count++;
    error = NR_ERROR;

    if (verbose) {
      printf("W=%e K=%e T=%e dx=%e dy=%e dz=%e error=%e\n",
             x,y,z,deltax,deltay,deltaz,error);
    }

    if (isnan(deltax) || isinf(deltax)) { deltax=0.0; error=1.0; }
    if (isnan(deltay) || isinf(deltay)) { deltay=0.0; error=1.0; }
    if (isnan(deltaz) || isinf(deltaz)) { deltaz=0.0; error=1.0; }

    x += deltax;
    y += deltay;
    z += deltaz;

    if (y < 0.0) y = 0.0;

    if (y>YMAX || count>=Nmax) {
      fail = 1;
    }
    ++Iterations;
  }

  if (fail==1) {
    if (verbose) {
      printf("Resetting...\n");
    }

    count = 0;
    fail  = 0;
    error = 1.0;
    y = (qdotb*qdotb + 2.0*Q2*D)/(B2*D*D + 2.0*D*D*D);
    z = Told;
    double a =  qdotn + .5*B2*(y+2.0)/(y+1.0);
    double b = -qdotb;
    x = - (b + a*a*a)/(pow(b*b,1./3.)+a*a);

    while (error>NR_TOL && fail==0) {

      fghWKT(x,y,z,&deltax,&deltay,&deltaz,qdotn,qdotb,Q2,B2,D,&err1,&err2,&err3);
      count++;
      error = NR_ERROR;

      if (verbose) {
        printf("W=%e K=%e T=%e dx=%e dy=%e dz = %e error = %e\n",
               x,y,z,deltax,deltay,deltaz,error);
      }

      if (isnan(deltax) || isinf(deltax)) { deltax=0.0; error=1.0; }
      if (isnan(deltay) || isinf(deltay)) { deltay=0.0; error=1.0; }
      if (isnan(deltaz) || isinf(deltaz)) { deltaz=0.0; error=1.0; }

      x += deltax;
      y += deltay;
      z += deltaz;

      if (y < 0.0) y = 0.0;

      if (y>YMAX || count>=Nmax2) {
        fail = 1;
      }
      ++Iterations;
    }
  }

  if (fail==1) {

    if (verbose) {
      printf("Resetting...\n");
    }

    count = 0;
    fail  = 0;
    z = eos->Temperature_u(D / *gamma, -qdotn-D-.5*B2);
    double P = eos->Pressure(D, z);

    x = -qdotn + P - .5*B2;
    y = 1.0e8;

    while (error>NR_TOL && fail==0) {

      fghWKT(x,y,z,&deltax,&deltay,&deltaz,qdotn,qdotb,Q2,B2,D,&err1,&err2,&err3);
      count++;
      error = NR_ERROR;

      if (verbose) {
        printf("W=%e K=%e T=%e dx=%e dy=%e dz = %e error = %e\n",
               x,y,z,deltax,deltay,deltaz,error);
      }

      if (isnan(deltax) || isinf(deltax)) { deltax=0.0; error=1.0; }
      if (isnan(deltay) || isinf(deltay)) { deltay=0.0; error=1.0; }
      if (isnan(deltaz) || isinf(deltaz)) { deltaz=0.0; error=1.0; }

      x += deltax;
      y += deltay;
      z += deltaz;

      if (y < 0.0) y = 0.0;

      if (y>YMAX || count>=Nmaxsafe) {
        fail = 1;
      }
    }
    ++Iterations;
  }

  for (count=0; count<Nextra; count++) {

    fghWKT(x,y,z,&deltax,&deltay,&deltaz,qdotn,qdotb,Q2,B2,D,&err1,&err2,&err3);
    double test = deltax + deltay + deltaz;

    if (isnan(test) || isinf(test)) {
      deltax = 0.0; deltay = 0.0; deltaz = 0.0;
    }

    x += deltax;
    y += deltay;
    z += deltaz;

    x = fabs(x);
    y = fabs(y);
    z = fabs(z);
  }

  if ((x < 0.0) || (y < 0.0)) { // || (z < 0.0)) {
    fail = 1;
  }
  if (fail==1) {
    if (verbose) {
      printf("c2p failed: W=%e K=%e T=%e qdotn=%e qdotb=%e Q2=%e B2=%e, D=%e, count=%d\n",
             x,y,z,qdotn,qdotb,Q2,B2,D,count);
    }
    return RMHD_C2P_MAXITER;
  }

  *W = x;
  *T = z;
  *gamma = sqrt(y+1.0);

  if (verbose) {
    printf("success! Final K=%e\n", y);
  }

  return 0;
}
