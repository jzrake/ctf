
#include <cstdio>
#include <cstring>
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

  double Omega = 1.0 - R/L > 0.0 ? 1.0 - R/L : 0.0;
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


void SourceTermsMagnetar::intercell_flux(double x[3], int dim, double *F)
{
  double P[8];
  double z = x[2];
  double r = sqrt(x[0]*x[0] + x[1]*x[1]); // cylindrical r
  double R = sqrt(r*r + z*z);             // spherical r
  double L = magnetar_radius;
  double C = light_cylinder;

  double Omega = 1.0 - R/L > 0.0 ? 1.0 - R/L : 0.0;
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

SourceTermsMagnetar::SourceTermsMagnetar() :
  magnetar_radius(0.1),
  field_strength(24.0),
  light_cylinder(1.0) { }

std::valarray<double> SourceTermsMagnetar::dUdt(const std::valarray<double> &Uin)
/*
 * Assumes that ConsToPrim has already been called by the real Godunov operator,
 * so that Mara->PrimitiveArray has the correct data.
 */
{
  this->prepare_integration();

  std::valarray<double> L(Uin.size());
  std::valarray<double> &P = Mara->PrimitiveArray;

  double *F = (double*) malloc(stride[0]*sizeof(double));
  double *G = (double*) malloc(stride[0]*sizeof(double));
  double *H = (double*) malloc(stride[0]*sizeof(double));
  int sx=stride[1],sy=stride[2],sz=stride[3];

  for (int i=0; i<stride[0]; ++i) {
    F[i] = 0.0;
    G[i] = 0.0;
    H[i] = 0.0;
  }

  for (int i=0; i<stride[0]; i+=NQ) {
    int N[3];
    absolute_index_to_3d(i/NQ, N);

    double x[3] = {
      Mara->domain->x_at(N[0]),
      Mara->domain->y_at(N[1]),
      Mara->domain->z_at(N[2]) };

    double xx[3] = {x[0], x[1], x[2]};
    double xy[3] = {x[0], x[1], x[2]};
    double xz[3] = {x[0], x[1], x[2]};

    xx[0] += 0.5*dx;
    xy[1] += 0.5*dy;
    xz[2] += 0.5*dz;

    this->intercell_flux(xx, 1, &F[i]);
    this->intercell_flux(xy, 2, &G[i]);
    this->intercell_flux(xz, 3, &H[i]);
  }
  Mara->fluid->ConstrainedTransport3d(F, G, H, stride);

  for (int i=sx; i<stride[0]; i+=NQ) {

    double dt = 1e-6; // fiducial time-step for approximating time-derivative
    double U0[8]; // conserved before source terms
    double U1[8]; // conserved after source terms
    double P0[8]; // primitive before and after source terms

    memcpy(P0, &P[i], 8 * sizeof(double));
    Mara->fluid->PrimToCons(P0, U0);

    P0[Bx] -= dt*((F[i+Bx]-F[i+Bx-sx])/dx + (G[i+Bx]-G[i+Bx-sy])/dy + (H[i+Bx]-H[i+Bx-sz])/dz);
    P0[By] -= dt*((F[i+By]-F[i+By-sx])/dx + (G[i+By]-G[i+By-sy])/dy + (H[i+By]-H[i+By-sz])/dz);
    P0[Bz] -= dt*((F[i+Bz]-F[i+Bz-sx])/dx + (G[i+Bz]-G[i+Bz-sy])/dy + (H[i+Bz]-H[i+Bz-sz])/dz);

    Mara->fluid->PrimToCons(P0, U1);

    for (int q=0; q<8; ++q) {
      L[i+q] += (U1[q] - U0[q]) / dt;
    }
  }

  free(F);
  free(G);
  free(H);

  return L;
}


SourceTermsWind::SourceTermsWind() { }

std::valarray<double> SourceTermsWind::dUdt(const std::valarray<double> &Uin)
{
  this->prepare_integration();

  std::valarray<double> Lsrc(Uin.size());
  //  std::valarray<double> &P = Mara->PrimitiveArray;

  for (int i=0; i<stride[0]; i+=NQ) {

    int N[3];
    absolute_index_to_3d(i/NQ, N);

    double x[3] = {
      Mara->domain->x_at(N[0]),
      Mara->domain->y_at(N[1]),
      Mara->domain->z_at(N[2]) };

    double L0[8];
    double L = 0.05;
    double R = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

    L0[ddd] = 1e1 * (1.0 - R/L > 0.0 ? 1.0 - R/L : 0.0);
    L0[tau] = 1e3 * (1.0 - R/L > 0.0 ? 1.0 - R/L : 0.0);

    for (int q=0; q<NQ; ++q) {
      Lsrc[i+q] += L0[q];
    }
  }

  return Lsrc;
}



#include <iostream>
#include <map>
#include "valman.hpp"
#include "rmhd.hpp"
#include "rmhd-c2p.h"

MagneticBubbleBoundary::MagneticBubbleBoundary()
  : rotation_profile(RIGID_ROTATION),
    rotation_radius(0.2)
{

}
int MagneticBubbleBoundary::set_rotation_radius(double r0)
{
  if (r0 < 0.0) {
    return 1;
  }
  else {
    rotation_radius = r0;
    return 0;
  }
}
int MagneticBubbleBoundary::set_rotation_profile(const char *profile_string)
{
  std::map<std::string, RotationProfile> profs;
  profs["rigid_rotation"        ] = RIGID_ROTATION;
  profs["rigid_in_keplarian_out"] = RIGID_IN_KEPLARIAN_OUT;
  std::map<std::string, RotationProfile>::iterator prof = profs.find(profile_string);
  if (prof == profs.end()) {
    fprintf(stderr, "[MagneticBubbleBoundary] error: unknown rotation profile"
            ": %s\n", profile_string);
    return 1;
  }
  else {
    return 0;
  }
}
void MagneticBubbleBoundary::set_bc_z0_wall(std::valarray<double> &U) const
{
  const int Nx = Mara->domain->get_N(1);
  const int Ny = Mara->domain->get_N(2);
  const int Ng = Mara->domain->get_Ng();

  const double r0 = rotation_radius;

  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  for (int i=0; i<Nx+2*Ng; ++i) {
    for (int j=0; j<Ny+2*Ng; ++j) {
      for (int k=0; k<Ng; ++k) {

        double x = Mara->domain->x_at(i);
        double y = Mara->domain->y_at(j);
        double r = sqrt(x*x + y*y);
        double Omega;

        switch (rotation_profile) {
        case RIGID_ROTATION:
          if (r < r0) {
            Omega = 1.0;
          }
          else {
            Omega = 0.0;
          }
          break;
        case RIGID_IN_KEPLARIAN_OUT:
          if (r < 0.5*r0) {
            Omega = 1.0;
          }
          else if (r > 2*r0) {
            Omega = 0.0;
          }
          else {
            Omega = pow(2*r/r0, -1.5);
          }
          break;
        }

        std::valarray<double> U1 = U[ M(i,j,2*Ng-k-1) ]; // reflected zone
        std::valarray<double> P1(8); // reflected zone primitive variables
        std::valarray<double> U0(8);
        std::valarray<double> P0(8);

        if (receive_primitive) {
          P1 = U1; // U was actually P (confusing, I know)
        }
        else {
          int err = Mara->fluid->ConsToPrim(&U1[0], &P1[0]);
          if (err) {
            /*
              fprintf(stderr,
              "[MagneticBubbleBoundary] unphysical boundary value\n");
              std::cerr << Mara->fluid->PrintCons(&U1[0]) << std::endl;
              std::cerr << rmhd_c2p_get_error(err) << std::endl;
            */
            //      exit(1);
            continue;
          }
        }

        /* not reflecting vz ensures that v < 1 in the guard zones */
        P0[0] =  P1[0];
        P0[1] =  P1[1];
        P0[2] = -y/r0 * Omega;
        P0[3] =  x/r0 * Omega;
        P0[4] =  0.0;//-P1[4];
        P0[5] = -P1[5];
        P0[6] = -P1[6];
        P0[7] =  P1[7];

        if (receive_primitive) {
          U0 = P0; // U really is P
          int err = rmhd_c2p_check_prim(&P0[0]);
          if (err) {
            std::cerr << rmhd_c2p_get_error(err) << std::endl;
            exit(1);
          }
        }
        else {
          int err = Mara->fluid->PrimToCons(&P0[0], &U0[0]);
          if (err) {
            std::cerr << rmhd_c2p_get_error(err) << std::endl;
            exit(1);
          }
        }
        U[ M(i,j,k) ] = U0;
      }
    }
  }
}


void WindInflowBoundary::set_bc_x0_wall(std::valarray<double> &U) const
{
  const int Ny = Mara->domain->get_N(2);
  const int Ng = Mara->domain->get_Ng();

  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  for (int i=0; i<Ng; ++i)
    for (int j=0; j<Ny+2*Ng; ++j)
      {
        U[ M(i,j) ] = U[ M(Ng,j) ];
      }
}
void WindInflowBoundary::set_bc_x1_wall(std::valarray<double> &U) const
{
  const int Nx = Mara->domain->get_N(1);
  const int Ny = Mara->domain->get_N(2);
  const int Ng = Mara->domain->get_Ng();

  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  for (int i=Nx+Ng; i<Nx+2*Ng; ++i)
    for (int j=0; j<Ny+2*Ng; ++j)
      {
        U[ M(i,j) ] = U[ M(Nx+Ng-1,j) ];
      }
}
void WindInflowBoundary::set_bc_y0_wall(std::valarray<double> &U) const
{
  const int Nx = Mara->domain->get_N(1);
  const int Ng = Mara->domain->get_Ng();

  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  for (int i=0; i<Nx+2*Ng; ++i)
    for (int j=0; j<Ng; ++j)
      {
        U[ M(i,j) ] = U[ M(i,Ng) ];
      }
}
void WindInflowBoundary::set_bc_y1_wall(std::valarray<double> &U) const
{
  const int Nx = Mara->domain->get_N(1);
  const int Ny = Mara->domain->get_N(2);
  const int Ng = Mara->domain->get_Ng();

  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  for (int i=0; i<Nx+2*Ng; ++i)
    for (int j=Ny+Ng; j<Ny+2*Ng; ++j)
      {
        U[ M(i,j) ] = U[ M(i,Ny+Ng-1) ];
      }
}
