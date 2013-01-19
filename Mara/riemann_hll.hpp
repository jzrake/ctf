
#ifndef __HllRiemanSolver_HEADER__
#define __HllRiemanSolver_HEADER__

#include "hydro.hpp"

class HllRiemannSolver : public RiemannSolver
{
public:
  int IntercellFlux(const double *pl, const double *pr, double *U,
		    double *F, double s, int dim);
} ;


#endif // __HllRiemanSolver_HEADER__
