

#ifndef __PlmCtuHancockOperator_HEADER__
#define __PlmCtuHancockOperator_HEADER__

#include "hydro.hpp"

class PlmCtuHancockOperator : public GodunovOperator
{
private:
  void reconstruct_plm(const double *P0, double *Pl, double *Pr, int S);
  void ctu_hancock(const double *P, const double *U, double *L, double dt);
  double TimeStepDt;

public:
  std::valarray<double> dUdt(const std::valarray<double> &Uin);
  void SetTimeStepDt(double dt) { TimeStepDt = dt; }
  void SetPlmTheta(double plm);
  void SetSafetyLevel(int level);
} ;

#endif // __PlmCtuHancockOperator_HEADER__
