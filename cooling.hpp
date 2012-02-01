
#ifndef __CoolingModule_HEADER__
#define __CoolingModule_HEADER__

#include "hydro.hpp"

class CoolingModuleT4 : public CoolingModule
{
private:
  const double Tref, t0;
  double energy_removed;

public:
  CoolingModuleT4(double Tref, double t0);
  void Cool(std::valarray<double> &P, double dt);
  double EnergyRemoved();
} ;

class CoolingModuleE4 : public CoolingModule
{
private:
  const double eref, t0;
  double energy_removed;

public:
  CoolingModuleE4(double eref, double t0);
  void Cool(std::valarray<double> &P, double dt);
  double EnergyRemoved();
} ;


#endif // __CoolingModule_HEADER__
