
#ifndef __RungeKutta_HEADER__
#define __RungeKutta_HEADER__

#include "hydro.hpp"

class RungeKuttaSingleStep : public RungeKuttaIntegration
{
public:
  void AdvanceState(std::valarray<double> &U, double dt) const
  {
    Mara->godunov->SetTimeStepDt(dt);
    U += dt * Mara->godunov->dUdt(U);
  }
} ;
class RungeKuttaRk2Tvd : public RungeKuttaIntegration
{
public:
  void AdvanceState(std::valarray<double> &U, double dt) const
  {
    std::valarray<double> U1(U.size());
    Mara->godunov->SetTimeStepDt(dt);

    U1 =      U +      dt*Mara->godunov->dUdt(U);
    U  = 0.5*(U + U1 + dt*Mara->godunov->dUdt(U1));
  }
} ;
class RungeKuttaShuOsherRk3 : public RungeKuttaIntegration
{
public:
  void AdvanceState(std::valarray<double> &U, double dt) const
  {
    std::valarray<double> U1(U.size());
    Mara->godunov->SetTimeStepDt(dt);

    U1 =      U +                  dt * Mara->godunov->dUdt(U );
    U1 = 3./4*U + 1./4*U1 + 1./4 * dt * Mara->godunov->dUdt(U1);
    U  = 1./3*U + 2./3*U1 + 2./3 * dt * Mara->godunov->dUdt(U1);
  }
} ;

#endif // __RungeKutta_HEADER__
