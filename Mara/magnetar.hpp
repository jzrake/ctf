
/*------------------------------------------------------------------------------
 * FILE: magnetar.hpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * REFERENCES:
 *
 *------------------------------------------------------------------------------
 */

#ifndef __Magnetar_HEADER__
#define __Magnetar_HEADER__

#include "hydro.hpp"
#include "boundary.hpp"

// MagneticBubbleBoundary
// -----------------------------------------------------------------------------
class MagneticBubbleBoundary : public OutflowBoundary3d
{
public:
  MagneticBubbleBoundary();
  int set_rotation_radius(double r0);
  int set_rotation_profile(const char *profile_string);
private:
  enum RotationProfile {RIGID_ROTATION, RIGID_IN_KEPLARIAN_OUT};
  RotationProfile rotation_profile;
  double rotation_radius;
  virtual void set_bc_z0_wall(std::valarray<double> &U) const;
} ;

class WindInflowBoundary : public OutflowBoundary2d
{
private:
  virtual void set_bc_x0_wall(std::valarray<double> &U) const;
  virtual void set_bc_x1_wall(std::valarray<double> &U) const;
  virtual void set_bc_y0_wall(std::valarray<double> &U) const;
  virtual void set_bc_y1_wall(std::valarray<double> &U) const;
} ;


class FluxSourceTermsMagnetar : public FluxSourceTermsModule
// -----------------------------------------------------------------------------
{
private:
  double magnetar_radius;
  double field_strength;
  double light_cylinder;
public:
  FluxSourceTermsMagnetar();
  virtual ~FluxSourceTermsMagnetar() { }
  virtual void AddIntercellFlux(double x[3], int dim, double *F);
  void set_magnetar_radius(double L) { magnetar_radius = L; }
  void set_field_strength(double B) { field_strength = B; }
  void set_light_cylinder(double C) { light_cylinder = C; }
} ;

class SourceTermsMagnetar : public SourceTermsModule
// -----------------------------------------------------------------------------
{
private:
  double magnetar_radius;
  double field_strength;
  double light_cylinder;
  void intercell_flux(double x[3], int dim, double *F);
public:
  SourceTermsMagnetar();
  std::valarray<double> dUdt(const std::valarray<double> &Uin);
  void set_magnetar_radius(double L) { magnetar_radius = L; }
  void set_field_strength(double B) { field_strength = B; }
  void set_light_cylinder(double C) { light_cylinder = C; }
} ;

class SourceTermsWind : public SourceTermsModule
// -----------------------------------------------------------------------------
{
public:
  SourceTermsWind();
  std::valarray<double> dUdt(const std::valarray<double> &Uin);
} ;

#endif // __Magnetar_HEADER__
