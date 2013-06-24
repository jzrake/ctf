
#include "magnetar.hpp"

enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
enum { rho, pre, vx, vy, vz };             // Primitive

FluxSourceTermsMagnetar::FluxSourceTermsMagnetar() :
  magnetar_radius(0.1),
  field_strength(24.0),
  light_cylinder(1.0) { }

void FluxSourceTermsMagnetar::
AddIntercellFlux(double x[3], int dim, double *F)
{
  double P[8];
  double z = x[2];
  double r = sqrt(x[0]*x[0] + x[1]*x[1]); // cylindrical r
  double R = sqrt(r*r + z*z);             // spherical r
  double L = magnetar_radius;
  double C = light_cylinder;

  double Omega = 1.0 / (1.0 + pow(R/L, 4));
  double vphi = Omega * (r / C);
  double phihat[3] = { -x[1]/r, x[0]/r, 0.0 };
  double m[3] = { 0.0, 0.0, 0.0005 };
  double rhat[3] = { x[0]/R, x[1]/R, x[2]/R };
  double rhatdotm = rhat[0]*m[0] + rhat[1]*m[1] + rhat[2]*m[2];

  P[rho] = 1.0;
  P[pre] = 1.0;
  P[vx] = vphi * phihat[0];
  P[vy] = vphi * phihat[1];
  P[vz] = 0.0;

  if (1) {
    P[Bx] = 0.0;
    P[By] = 0.0;
    P[Bz] = field_strength;
  }
  else if (0) {
    P[Bx] = rhat[0] / (R*R);
    P[By] = rhat[1] / (R*R);
    P[Bz] = rhat[2] / (R*R) * (z < 0 ? -1 : 1);
  }
  else if (0) {
    P[Bx] = (3*rhatdotm*rhat[0] - m[0]) / (R*R*R);
    P[By] = (3*rhatdotm*rhat[1] - m[1]) / (R*R*R);
    P[Bz] = (3*rhatdotm*rhat[2] - m[2]) / (R*R*R);
  }

  /*
  double U[8];
  Mara->fluid->PrimToCons(P, U);
  Mara->fluid->FluxAndEigenvalues(U, P, F, NULL, NULL, dim);
  */

  switch (dim) {
  case 1:
    F[Bx] += 0.0;
    F[By] += P[vx]*P[By] - P[vy]*P[Bx];
    F[Bz] += P[vx]*P[Bz] - P[vz]*P[Bx];
    break;
  case 2:
    F[Bx] += P[vy]*P[Bx] - P[vx]*P[By];
    F[By] += 0.0;
    F[Bz] += P[vy]*P[Bz] - P[vz]*P[By];
    break;
  case 3:
    F[Bx] += P[vz]*P[Bx] - P[vx]*P[Bz];
    F[By] += P[vz]*P[By] - P[vy]*P[Bz];
    F[Bz] += 0.0;
    break;
  }
}
