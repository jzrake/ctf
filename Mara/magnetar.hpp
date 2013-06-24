
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

#endif // __Magnetar_HEADER__
