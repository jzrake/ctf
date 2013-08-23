
/*------------------------------------------------------------------------------
 * FILE: boundary.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * REFERENCES:
 *
 *------------------------------------------------------------------------------
 */

#include "boundary.hpp"
#include "valman.hpp"


void OutflowBoundary::ApplyBoundaries(std::valarray<double> &U) const
{
  BoundaryConditions *boundary=NULL;
  switch (Mara->domain->get_Nd()) {
  case 1: boundary = new OutflowBoundary1d; break;
  case 2: boundary = new OutflowBoundary2d; break;
  case 3: boundary = new OutflowBoundary3d; break;
  }
  boundary->ApplyBoundaries(U);
  delete boundary;
}


// OutflowBoundary1d
// -----------------------------------------------------------------------------
void OutflowBoundary1d::ApplyBoundaries(std::valarray<double> &U) const
{
  Mara->domain->Synchronize(U);

  if (Mara->domain->GetSubgridIndex(0) == 0)
    set_bc_x0_wall(U);
  if (Mara->domain->GetSubgridIndex(0) == Mara->domain->GetSubgridSizes(0)-1)
    set_bc_x1_wall(U);
}
void OutflowBoundary1d::set_bc_x0_wall(std::valarray<double> &U) const
{
  const int Ng = Mara->domain->get_Ng();

  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  for (int i=0; i<Ng; ++i)
    {
      U[ M(i) ] = U[ M(Ng) ];
    }
}
void OutflowBoundary1d::set_bc_x1_wall(std::valarray<double> &U) const
{
  const int Nx = Mara->domain->get_N(1);
  const int Ng = Mara->domain->get_Ng();

  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  for (int i=Nx+Ng; i<Nx+2*Ng; ++i)
    {
      U[ M(i) ] = U[ M(Nx+Ng-1) ];
    }
}


// OutflowBoundary2d
// -----------------------------------------------------------------------------
void OutflowBoundary2d::ApplyBoundaries(std::valarray<double> &U) const
{
  Mara->domain->Synchronize(U);

  if (Mara->domain->GetSubgridIndex(0) == 0)
    set_bc_x0_wall(U);
  if (Mara->domain->GetSubgridIndex(0) == Mara->domain->GetSubgridSizes(0)-1)
    set_bc_x1_wall(U);

  if (Mara->domain->GetSubgridIndex(1) == 0)
    set_bc_y0_wall(U);
  if (Mara->domain->GetSubgridIndex(1) == Mara->domain->GetSubgridSizes(1)-1)
    set_bc_y1_wall(U);
}

void OutflowBoundary2d::set_bc_x0_wall(std::valarray<double> &U) const
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
void OutflowBoundary2d::set_bc_x1_wall(std::valarray<double> &U) const
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
void OutflowBoundary2d::set_bc_y0_wall(std::valarray<double> &U) const
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
void OutflowBoundary2d::set_bc_y1_wall(std::valarray<double> &U) const
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


// PeriodicXOutflowY2d
// -----------------------------------------------------------------------------
void PeriodicXOutflowY2d::ApplyBoundaries(std::valarray<double> &U) const
{
  Mara->domain->Synchronize(U);

  if (Mara->domain->GetSubgridIndex(0) == 0)
    set_bc_x0_wall(U);
  if (Mara->domain->GetSubgridIndex(0) == Mara->domain->GetSubgridSizes(0)-1)
    set_bc_x1_wall(U);

  if (Mara->domain->GetSubgridIndex(1) == 0)
    set_bc_y0_wall(U);
  if (Mara->domain->GetSubgridIndex(1) == Mara->domain->GetSubgridSizes(1)-1)
    set_bc_y1_wall(U);
}

void PeriodicXOutflowY2d::set_bc_x0_wall(std::valarray<double> &U) const
{
  // periodic in x
}
void PeriodicXOutflowY2d::set_bc_x1_wall(std::valarray<double> &U) const
{
  // periodic in x
}
void PeriodicXOutflowY2d::set_bc_y0_wall(std::valarray<double> &U) const
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
void PeriodicXOutflowY2d::set_bc_y1_wall(std::valarray<double> &U) const
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


// OutflowBoundary3d
// -----------------------------------------------------------------------------
void OutflowBoundary3d::ApplyBoundaries(std::valarray<double> &U) const
{
  Mara->domain->Synchronize(U);

  if (Mara->domain->GetSubgridIndex(0) == 0)
    set_bc_x0_wall(U);
  if (Mara->domain->GetSubgridIndex(0) == Mara->domain->GetSubgridSizes(0)-1)
    set_bc_x1_wall(U);

  if (Mara->domain->GetSubgridIndex(1) == 0)
    set_bc_y0_wall(U);
  if (Mara->domain->GetSubgridIndex(1) == Mara->domain->GetSubgridSizes(1)-1)
    set_bc_y1_wall(U);

  if (Mara->domain->GetSubgridIndex(2) == 0)
    set_bc_z0_wall(U);
  if (Mara->domain->GetSubgridIndex(2) == Mara->domain->GetSubgridSizes(2)-1)
    set_bc_z1_wall(U);
}

void OutflowBoundary3d::set_bc_x0_wall(std::valarray<double> &U) const
{
  const int Ny = Mara->domain->get_N(2);
  const int Nz = Mara->domain->get_N(3);
  const int Ng = Mara->domain->get_Ng();

  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  for (int i=0; i<Ng; ++i)
    for (int j=0; j<Ny+2*Ng; ++j)
      for (int k=0; k<Nz+2*Ng; ++k)
        {
          U[ M(i,j,k) ] = U[ M(Ng,j,k) ];
        }
}
void OutflowBoundary3d::set_bc_x1_wall(std::valarray<double> &U) const
{
  const int Nx = Mara->domain->get_N(1);
  const int Ny = Mara->domain->get_N(2);
  const int Nz = Mara->domain->get_N(3);
  const int Ng = Mara->domain->get_Ng();

  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  for (int i=Nx+Ng; i<Nx+2*Ng; ++i)
    for (int j=0; j<Ny+2*Ng; ++j)
      for (int k=0; k<Nz+2*Ng; ++k)
        {
          U[ M(i,j,k) ] = U[ M(Nx+Ng-1,j,k) ];
        }
}
void OutflowBoundary3d::set_bc_y0_wall(std::valarray<double> &U) const
{
  const int Nx = Mara->domain->get_N(1);
  const int Nz = Mara->domain->get_N(3);
  const int Ng = Mara->domain->get_Ng();

  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  for (int i=0; i<Nx+2*Ng; ++i)
    for (int j=0; j<Ng; ++j)
      for (int k=0; k<Nz+2*Ng; ++k)
        {
          U[ M(i,j,k) ] = U[ M(i,Ng,k) ];
        }
}
void OutflowBoundary3d::set_bc_y1_wall(std::valarray<double> &U) const
{
  const int Nx = Mara->domain->get_N(1);
  const int Ny = Mara->domain->get_N(2);
  const int Nz = Mara->domain->get_N(3);
  const int Ng = Mara->domain->get_Ng();

  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  for (int i=0; i<Nx+2*Ng; ++i)
    for (int j=Ny+Ng; j<Ny+2*Ng; ++j)
      for (int k=0; k<Nz+2*Ng; ++k)
        {
          U[ M(i,j,k) ] = U[ M(i,Ny+Ng-1,k) ];
        }
}
void OutflowBoundary3d::set_bc_z0_wall(std::valarray<double> &U) const
{
  const int Nx = Mara->domain->get_N(1);
  const int Ny = Mara->domain->get_N(2);
  const int Ng = Mara->domain->get_Ng();

  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  for (int i=0; i<Nx+2*Ng; ++i)
    for (int j=0; j<Ny+2*Ng; ++j)
      for (int k=0; k<Ng; ++k)
        {
          U[ M(i,j,k) ] = U[ M(i,j,Ng) ];
        }
}
void OutflowBoundary3d::set_bc_z1_wall(std::valarray<double> &U) const
{
  const int Nx = Mara->domain->get_N(1);
  const int Ny = Mara->domain->get_N(2);
  const int Nz = Mara->domain->get_N(3);
  const int Ng = Mara->domain->get_Ng();

  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  for (int i=0; i<Nx+2*Ng; ++i)
    for (int j=0; j<Ny+2*Ng; ++j)
      for (int k=Nz+Ng; k<Nz+2*Ng; ++k)
        {
          U[ M(i,j,k) ] = U[ M(i,j,Nz+Ng-1) ];
        }
}


// ReflectingBoundary2d
// -----------------------------------------------------------------------------
ReflectingBoundary2d::ReflectingBoundary2d(int IndexReverseX, int IndexReverseY)
  : IndexReverseX(IndexReverseX), IndexReverseY(IndexReverseY)
{

}
void ReflectingBoundary2d::ApplyBoundaries(std::valarray<double> &U) const
{
  Mara->domain->Synchronize(U);

  if (Mara->domain->GetSubgridIndex(0) == 0)
    set_bc_x0_wall(U);
  if (Mara->domain->GetSubgridIndex(0) == Mara->domain->GetSubgridSizes(0)-1)
    set_bc_x1_wall(U);

  if (Mara->domain->GetSubgridIndex(1) == 0)
    set_bc_y0_wall(U);
  if (Mara->domain->GetSubgridIndex(1) == Mara->domain->GetSubgridSizes(1)-1)
    set_bc_y1_wall(U);
}
void ReflectingBoundary2d::set_bc_x0_wall(std::valarray<double> &U) const
{
  const int Ny = Mara->domain->get_N(2);
  const int Ng = Mara->domain->get_Ng();

  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  for (int i=0; i<Ng; ++i)
    for (int j=0; j<Ny+2*Ng; ++j)
      {
        std::valarray<double> U0 = U[ M(Ng,j) ];
        U0[IndexReverseX] *= -1.0;
        U[ M(i,j) ] = U0;
      }
}
void ReflectingBoundary2d::set_bc_x1_wall(std::valarray<double> &U) const
{
  const int Nx = Mara->domain->get_N(1);
  const int Ny = Mara->domain->get_N(2);
  const int Ng = Mara->domain->get_Ng();

  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  for (int i=Nx+Ng; i<Nx+2*Ng; ++i)
    for (int j=0; j<Ny+2*Ng; ++j)
      {
        std::valarray<double> U0 = U[ M(Nx+Ng-1,j) ];
        U0[IndexReverseX] *= -1.0;
        U[ M(i,j) ] = U0;
      }
}
void ReflectingBoundary2d::set_bc_y0_wall(std::valarray<double> &U) const
{
  const int Nx = Mara->domain->get_N(1);
  const int Ng = Mara->domain->get_Ng();

  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  for (int i=0; i<Nx+2*Ng; ++i)
    for (int j=0; j<Ng; ++j)
      {
        std::valarray<double> U0 = U[ M(i,Ng) ];
        U0[IndexReverseY] *= -1.0;
        U[ M(i,j) ] = U0;
      }
}
void ReflectingBoundary2d::set_bc_y1_wall(std::valarray<double> &U) const
{
  const int Nx = Mara->domain->get_N(1);
  const int Ny = Mara->domain->get_N(2);
  const int Ng = Mara->domain->get_Ng();

  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  for (int i=0; i<Nx+2*Ng; ++i)
    for (int j=Ny+Ng; j<Ny+2*Ng; ++j)
      {
        std::valarray<double> U0 = U[ M(i,Ny+Ng-1) ];
        U0[IndexReverseY] *= -1.0;
        U[ M(i,j) ] = U0;
      }
}

#include <cstdio>
#include <iostream>
#include <map>
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
