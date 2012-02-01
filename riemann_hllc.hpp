
#ifndef __HllcRiemannSolver_HEADER__
#define __HllcRiemannSolver_HEADER__

#include <typeinfo>
#include "eulers.hpp"
#include "srhd.hpp"
#include "rmhd.hpp"


class HllcRiemannSolver : public RiemannSolver
{
public:
  int IntercellFlux(const double *pl, const double *pr, double *U,
		    double *F, double s, int dim);
} ;

class HllcEulersRiemannSolver : public RiemannSolver
{
public:
  int IntercellFlux(const double *pl, const double *pr, double *U,
		    double *F, double s, int dim);
} ;

class HllcSrhdRiemannSolver : public RiemannSolver
{
public:
  int IntercellFlux(const double *pl, const double *pr, double *U,
                    double *F, double s, int dim);
} ;

class HllcRmhdRiemannSolver : public RiemannSolver
{
public:
  int IntercellFlux(const double *pl, const double *pr, double *U,
		    double *F, double s, int dim);
} ;


#endif // __HllcRmhdRiemannSolver_HEADER__
