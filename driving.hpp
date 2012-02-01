


#ifndef __DrivingProcedure_HEADER__
#define __DrivingProcedure_HEADER__

#include "hydro.hpp"
#include "ou-field.hpp"

class DrivingProcedure : public DrivingModule
{
protected:
  const int num_dims;
  StochasticVectorField *field;
  std::valarray<double> Fx, Fy, Fz;
  double InjectionRate;
  double CoolingRatio;

public:
  DrivingProcedure(StochasticVectorField *field);
  ~DrivingProcedure();
  StochasticVectorField *GetField();
  double AveragePowerInFields() const;
  double EnergyInjectionRate() const;
  void ResampleField();
  void Drive(std::valarray<double> &P, double dt);

private:
  void Drive_eulers(std::valarray<double> &P, double dt);
  void Drive_srhd  (std::valarray<double> &P, double dt);
  void Drive_rmhd  (std::valarray<double> &P, double dt);
} ;

#endif // __DrivingProcedure_HEADER__
