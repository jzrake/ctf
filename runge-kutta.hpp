
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
class RungeKuttaClassicRk4 : public RungeKuttaIntegration
{
public:
  void AdvanceState(std::valarray<double> &U, double dt) const
  {
    Mara->godunov->SetTimeStepDt(dt);

    std::valarray<double> L1 = dt * Mara->godunov->dUdt(U);
    std::valarray<double> L2 = dt * Mara->godunov->dUdt(U + 0.5*L1);
    std::valarray<double> L3 = dt * Mara->godunov->dUdt(U + 0.5*L2);
    std::valarray<double> L4 = dt * Mara->godunov->dUdt(U + 1.0*L3);

    U += (1.0/6.0) * (L1 + 2.0*L2 + 2.0*L3 + L4);
  }
} ;

#endif // __RungeKutta_HEADER__
