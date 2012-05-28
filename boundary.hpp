
/*------------------------------------------------------------------------------
 * FILE: boundary.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * REFERENCES:
 *
 *------------------------------------------------------------------------------
 */


#ifndef __Boundary_HEADER__
#define __Boundary_HEADER__

#include <valarray>
#include "hydro.hpp"


// OutflowBoundary1d
// -----------------------------------------------------------------------------
class OutflowBoundary1d : public BoundaryConditions
{
public:
  void ApplyBoundaries(std::valarray<double> &U) const;

protected:
  void set_bc_x0_wall(std::valarray<double> &U) const;
  void set_bc_x1_wall(std::valarray<double> &U) const;
} ;


// OutflowBoundary2d
// -----------------------------------------------------------------------------
class OutflowBoundary2d : public BoundaryConditions
{
public:
  void ApplyBoundaries(std::valarray<double> &U) const;

private:
  void set_bc_x0_wall(std::valarray<double> &U) const;
  void set_bc_x1_wall(std::valarray<double> &U) const;
  void set_bc_y0_wall(std::valarray<double> &U) const;
  void set_bc_y1_wall(std::valarray<double> &U) const;
} ;


// PeriodicXOutflowY2d
// -----------------------------------------------------------------------------
class PeriodicXOutflowY2d : public BoundaryConditions
{
public:
  void ApplyBoundaries(std::valarray<double> &U) const;

private:
  void set_bc_x0_wall(std::valarray<double> &U) const;
  void set_bc_x1_wall(std::valarray<double> &U) const;
  void set_bc_y0_wall(std::valarray<double> &U) const;
  void set_bc_y1_wall(std::valarray<double> &U) const;
} ;


// OutflowBoundary3d
// -----------------------------------------------------------------------------
class OutflowBoundary3d : public BoundaryConditions
{
public:
  void ApplyBoundaries(std::valarray<double> &U) const;

private:
  void set_bc_x0_wall(std::valarray<double> &U) const;
  void set_bc_x1_wall(std::valarray<double> &U) const;
  void set_bc_y0_wall(std::valarray<double> &U) const;
  void set_bc_y1_wall(std::valarray<double> &U) const;
  void set_bc_z0_wall(std::valarray<double> &U) const;
  void set_bc_z1_wall(std::valarray<double> &U) const;
} ;


class OutflowBoundary : public BoundaryConditions
{
public:
  void ApplyBoundaries(std::valarray<double> &U) const;
} ;


// ReflectingBoundary2d
// -----------------------------------------------------------------------------
class ReflectingBoundary2d : public BoundaryConditions
{
public:
  ReflectingBoundary2d(int IndexReverseX, int IndexReverseY);
  void ApplyBoundaries(std::valarray<double> &U) const;

protected:
  int IndexReverseX, IndexReverseY;

  void set_bc_x0_wall(std::valarray<double> &U) const;
  void set_bc_x1_wall(std::valarray<double> &U) const;
  void set_bc_y0_wall(std::valarray<double> &U) const;
  void set_bc_y1_wall(std::valarray<double> &U) const;
} ;

class PeriodicBoundary : public BoundaryConditions
{
public:
  void ApplyBoundaries(std::valarray<double> &U) const
  {
    Mara->domain->Synchronize(U);
  }
} ;

#endif // __Boundary_HEADER__
