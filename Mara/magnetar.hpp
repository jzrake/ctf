
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
public:
  virtual ~FluxSourceTermsMagnetar() { }
  virtual void AddIntercellFlux(double x[3], int dim, double *F);
} ;

#endif // __Magnetar_HEADER__
