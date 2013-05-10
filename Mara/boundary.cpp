
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

#include "rmhd.hpp"
void MagneticBubbleBoundary::set_bc_z0_wall(std::valarray<double> &U) const
{
  const int Nx = Mara->domain->get_N(1);
  const int Ny = Mara->domain->get_N(2);
  const int Ng = Mara->domain->get_Ng();

  double r0 = 0.05;

  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  for (int i=0; i<Nx+2*Ng; ++i) {
    for (int j=0; j<Ny+2*Ng; ++j) {
      for (int k=0; k<Ng; ++k) {

	double x = Mara->domain->x_at(i);
	double y = Mara->domain->y_at(j);
	double r = sqrt(x*x + y*y) / r0;

	double Omega = r*r / (1 + r*r*r*r);

        std::valarray<double> U1 = U[ M(i,j,2*Ng-k-1) ]; // reflected zone
	std::valarray<double> P1(8); // reflected zone primitive variables
	std::valarray<double> U0(8);
	std::valarray<double> P0(8);

	Mara->fluid->ConsToPrim(&U1[0], &P1[0]);

	P0[0] =  P1[0];
	P0[1] =  P1[1];
	P0[2] = -y/r0 * Omega;
	P0[3] =  x/r0 * Omega;
	P0[4] = -P1[4];
	P0[5] = -P1[5];
	P0[6] = -P1[6];
	P0[7] =  P1[7];

	//	double v = sqrt(P0[2]*P0[2] + P0[3]*P0[3]);
	//	if (v > 0.5) printf("%f\n", v);

	Mara->fluid->PrimToCons(&P0[0], &U0[0]);
	U[ M(i,j,k) ] = U0;
      }
    }
  }
}
