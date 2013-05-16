
/*------------------------------------------------------------------------------
 * FILE: rmhd-c2p.c
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
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
 * This module solves for the primitive state (rho, pre, v, B) given the
 * conserved state (ddd, tau, S, B). In other words it solves for the 5 unknowns
 * rho, pre, vx, vy, vz. An adiabatic equation of state is assumed throughout,
 * with the index being provided by rmhd_c2p_set_gamma(). The solve functions
 * promise not to modify the pointer to result primitives unless the execution
 * is successful.
 *
 * ------------------------------------------------------------------------------
 */


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "rmhd-c2p.h"


enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
enum { rho, pre, vx, vy, vz };             // Primitive


static const int MaxIteration = 250;
static const double Tolerance = 1e-12;
static const double bigZ = 1e20;
static const double bigW = 1e12;
static const double smlZ = 0.0;
static const double smlW = 1.0;

static int Iterations = 0;
static double AdiabaticGamma = 1.4;
static double gamf;
static double D,Tau;
static double S2,B2,BS,BS2;
static double Cons[8];
static double Z_start, W_start;
static double PressureFloor = -1.0;



// Set method for the adiabatic index (defaults to 1.4)
// -----------------------------------------------------------------------------
void rmhd_c2p_set_gamma(double adiabatic_gamma)
{
  AdiabaticGamma = adiabatic_gamma;
  gamf = (AdiabaticGamma - 1.0) / AdiabaticGamma;
}

// Get method for the iterations on the last execution
// -----------------------------------------------------------------------------
int rmhd_c2p_get_iterations()
{
  return Iterations;
}

void rmhd_c2p_set_pressure_floor(double pf)
{
  PressureFloor = pf;
}

// Provide a new conserved state in memory. Before the solver is executed, the
// user must provide a guess for the initial primitive state, by calling either
// estimate_from_cons() or set_starting_prim(P).
// -----------------------------------------------------------------------------
void rmhd_c2p_new_state(const double *U)
{
  D    = U[ddd];
  Tau  = U[tau];
  S2   = U[Sx]*U[Sx] + U[Sy]*U[Sy] + U[Sz]*U[Sz];
  B2   = U[Bx]*U[Bx] + U[By]*U[By] + U[Bz]*U[Bz];
  BS   = U[Bx]*U[Sx] + U[By]*U[Sy] + U[Bz]*U[Sz];
  BS2  = BS*BS;

  memcpy(Cons, U, 8*sizeof(double));
}

// This estimate becomes exact for no magnetic field in the NR limit.
// -----------------------------------------------------------------------------
void rmhd_c2p_estimate_from_cons()
{
  Z_start = sqrt(S2 + D*D);
  W_start = Z_start / D;
}
void rmhd_c2p_get_starting_prim(double *P)
{
  rmhd_c2p_reconstruct_prim(Z_start, W_start, P);
}

// Explicitly provide a guess to be used at the initial iteration.
// -----------------------------------------------------------------------------
void rmhd_c2p_set_starting_prim(const double *P)
{
  const double V2 = P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz];
  const double W2 = 1.0 / (1.0 - V2);
  const double e = P[pre] / (P[rho] * (AdiabaticGamma - 1.0));
  const double h = 1.0 + e + P[pre] / P[rho];

  Z_start = P[rho] * h * W2;
  W_start = sqrt(W2);
}

// Using Z=rho*h*W^2, and W, get the primitive variables.
// -----------------------------------------------------------------------------
int rmhd_c2p_reconstruct_prim(double Z, double W, double *Pout)
{
  double P[8]; // Place result into temporary prim state for now.
  const double b0 = BS * W / Z;

  P[rho] =  D/W;
  P[pre] = (D/W) * (Z/(D*W) - 1.0) * gamf;
  P[vx ] = (Cons[Sx] + b0*Cons[Bx]/W) / (Z+B2);
  P[vy ] = (Cons[Sy] + b0*Cons[By]/W) / (Z+B2);
  P[vz ] = (Cons[Sz] + b0*Cons[Bz]/W) / (Z+B2);
  P[Bx ] =  Cons[Bx];
  P[By ] =  Cons[By];
  P[Bz ] =  Cons[Bz];

  if (P[pre] < PressureFloor && PressureFloor > 0.0) {
    //    printf("setting pressure floor, p=%4.3e -> %2.1e\n",
    //	   P[pre], PressureFloor);
    P[pre] = PressureFloor;    
  }

  int prim_error = rmhd_c2p_check_prim(P);
  if (prim_error) return prim_error;

  // Don't actually modify the user's prim state until we're sure the state is
  // good.
  // ---------------------------------------------------------------------------
  memcpy(Pout, P, 8*sizeof(double));
  return RMHD_C2P_SUCCESS;
}

// Routines to verify the health of conserved or primitive states
// -----------------------------------------------------------------------------
int rmhd_c2p_check_cons(const double *U)
{
  int i;
  if (U[ddd] < 0.0) return RMHD_C2P_CONS_NEGATIVE_DENSITY;
  if (U[tau] < 0.0) return RMHD_C2P_CONS_NEGATIVE_ENERGY;
  for (i=0; i<8; ++i) {
    if (isnan(U[i])) {
      return RMHD_C2P_CONS_CONTAINS_NAN;
    }
  }
  return RMHD_C2P_SUCCESS;
}
int rmhd_c2p_check_prim(const double *P)
{
  int i;
  const double v2 = P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz];
  if (v2   >=  1.0) return RMHD_C2P_PRIM_SUPERLUMINAL;
  if (P[pre] < 0.0) return RMHD_C2P_PRIM_NEGATIVE_PRESSURE;
  if (P[rho] < 0.0) return RMHD_C2P_PRIM_NEGATIVE_RESTMASS;
  for (i=0; i<8; ++i) {
    if (isnan(P[i])) {
      return RMHD_C2P_PRIM_CONTAINS_NAN;
    }
  }
  return RMHD_C2P_SUCCESS;
}

// Solution based on Anton & Zanotti (2006), equations 84 and 85.
// -----------------------------------------------------------------------------
int rmhd_c2p_solve_anton2dzw(double *P)
{
  int bad_input = rmhd_c2p_check_cons(Cons);
  if (bad_input) {
    return bad_input;
  }
  // Starting values
  // ---------------------------------------------------------------------------
  Iterations = 0;
  double error = 1.0;
  double W = W_start;
  double Z = Z_start;

  while (error > Tolerance) {

    const double Z2 = Z*Z;
    const double Z3 = Z*Z2;
    const double W2 = W*W;
    const double W3 = W*W2;
    const double Pre = (D/W) * (Z/(D*W) - 1.0) * gamf;

    const double df0dZ = 2*(B2+Z)*(BS2*W2 + (W2-1)*Z3) / (W2*Z3);
    const double df0dW = 2*(B2+Z)*(B2+Z) / W3;
    const double df1dZ = 1.0 + BS2/Z3 - gamf/W2;
    const double df1dW = B2/W3 + (2*Z - D*W)/W3 * gamf;

    double f[2];
    double J[4], G[4];


    // Evaluation of the function, and its Jacobian
    // -------------------------------------------------------------------------
    f[0] = -S2  + (Z+B2)*(Z+B2)*(W2-1)/W2 - (2*Z+B2)*BS2/Z2;      // eqn (84)
    f[1] = -Tau +  Z+B2 - Pre - 0.5*B2/W2 -      0.5*BS2/Z2 - D;  // eqn (85)

    J[0] = df0dZ; J[1] = df0dW;
    J[2] = df1dZ; J[3] = df1dW;

    // G in the inverse Jacobian
    // -------------------------------------------------------------------------
    const double det = J[0]*J[3] - J[1]*J[2];
    G[0] =  J[3]/det; G[1] = -J[1]/det;
    G[2] = -J[2]/det; G[3] =  J[0]/det;      // G = J^{-1}

    const double dZ = -(G[0]*f[0] + G[1]*f[1]); // Matrix multiply, dx = -G . f
    const double dW = -(G[2]*f[0] + G[3]*f[1]);

    // Bracketing the root
    // -------------------------------------------------------------------------
    double Z_new = Z + dZ;
    double W_new = W + dW;

    Z_new = (Z_new > smlZ) ? Z_new : -Z_new;
    Z_new = (Z_new < bigZ) ? Z_new :  Z;

    W_new = (W_new > smlW) ? W_new : smlW;
    W_new = (W_new < bigW) ? W_new : bigW;

    Z = Z_new;
    W = W_new;

    // -------------------------------------------------------------------------
    error = fabs(dZ/Z) + fabs(dW/W);
    ++Iterations;

    if (Iterations == MaxIteration) {
      return RMHD_C2P_MAXITER;
    }
  }

  return rmhd_c2p_reconstruct_prim(Z, W, P);
}


// Solution based on Noble et. al. (2006), using Z = rho h W^2 as the single
// unkown. Unfortunately, Noble uses 'W' for what I call Z. This function should
// really be called '1dz', but I use this name to reflect the name of the
// section in which it appears.
// -----------------------------------------------------------------------------
int rmhd_c2p_solve_noble1dw(double *P)
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

  double f, g;

  while (error > Tolerance) {

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
    const double Pre = (D/W) * (Z/(D*W) - 1.0) * gamf;

    const double dv2dZ    = (ap*b - bp*a) / (b*b); // (a'b - b'a) / b^2
    const double delPdelZ = gamf/W2;
    const double delPdelW = gamf*(D/W2 - 2*Z/W3);
    const double dWdv2    = 0.5*W3;
    const double dPdZ     = delPdelW * dWdv2 * dv2dZ + delPdelZ;

    f = Tau + D - 0.5*B2*(1+V2) + 0.5*BS2/Z2 - Z + Pre; // equation (29)
    g = -0.5*B2*dv2dZ - BS2/Z3 - 1.0 + dPdZ;

    const double dZ = -f/g;

    // -------------------------------------------------------------------------
    double Z_new = Z + dZ;

    Z_new = (Z_new > smlZ) ? Z_new : -Z_new;
    Z_new = (Z_new < bigZ) ? Z_new :  Z;

    Z = Z_new;

    error = fabs(dZ/Z);
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

  return rmhd_c2p_reconstruct_prim(Z, W, P);
}

char *rmhd_c2p_get_error(int error)
{
  switch (error) {
  case RMHD_C2P_SUCCESS:
    return "Successful inversion from conserved to primitive.";
  case RMHD_C2P_CONS_CONTAINS_NAN:
    return "Input conserved state contained nan's.";
  case RMHD_C2P_CONS_NEGATIVE_DENSITY:
    return "Input conserved state has negative density.";
  case RMHD_C2P_CONS_NEGATIVE_ENERGY:
    return "Input conserved state has negative energy.";
  case RMHD_C2P_PRIM_CONTAINS_NAN:
    return "Derived primitive state contains nan's.";
  case RMHD_C2P_PRIM_NEGATIVE_PRESSURE:
    return "Derived primitive state has negative pressure.";
  case RMHD_C2P_PRIM_NEGATIVE_RESTMASS:
    return "Derived primitive state contains negative density.";
  case RMHD_C2P_PRIM_SUPERLUMINAL:
    return "Derived primitive state contains superluminal velocity.";
  case RMHD_C2P_MAXITER:
    return "Rootfinder hit maximum iterations.";
  default:
    return "Unkown error.";
  }
}
