/*------------------------------------------------------------------------------
 * FILE: quartic.c
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP (adapted from Mathematica code)
 *
 * DESCRIPTION: Implementation of exact solution to the quartic polynomial
 *
 *------------------------------------------------------------------------------
 */


#include <math.h>
#include <complex.h>

typedef double complex Complex;

int m2_solve_quartic_equation(double d4, double d3,
			      double d2, double d1, double d0,
			      double roots[4])
{
  Complex a3 = d3 / d4;
  Complex a2 = d2 / d4;
  Complex a1 = d1 / d4;
  Complex a0 = d0 / d4;
  Complex W = 1./3;
  Complex X = 12*a0 - 3*a1*a3 + a2*a2;
  Complex P = -72*a0*a2 - 9*a1*a2*a3 + 27*a1*a1 + 27*a0*a3*a3 + 2*a2*a2*a2;
  Complex Q = cpow(X,3);
  Complex S = csqrt(-4*Q + cpow(P,2));
  Complex T = -8*a1 + 4*a2*a3 - a3*a3*a3;
  Complex B = cabs(P + S) < 1e-15 ? 0.0 : (cpow(2,W)*X)/(3.*cpow(P + S,W));
  Complex U = (-2*a2)/3. + (a3*a3)/4. + B;
  Complex C = csqrt(U + cpow(P + S,W)/(3.*cpow(2,W)))/2.;
  Complex D = cpow(P + S,W)/(3.*cpow(2,W));
  Complex E = T/(4.*csqrt(U + D));
  Complex F = csqrt((-4*a2)/3. + (a3*a3)/2. - B - D - E)/2.;
  Complex G = csqrt((-4*a2)/3. + (a3*a3)/2. - B - D + E)/2.;
  Complex r0 = -a3/4. - C - F;
  Complex r1 = -a3/4. - C + F;
  Complex r2 = -a3/4. + C - G;
  Complex r3 = -a3/4. + C + G;

  roots[0] = creal(r0);
  roots[1] = creal(r1);
  roots[2] = creal(r2);
  roots[3] = creal(r3);

  if (roots[0] != roots[0] ||
      roots[1] != roots[1] ||
      roots[2] != roots[2] ||
      roots[3] != roots[3]) {
    /* set breakpoint */
    return 0;
  }

  /* int nr = 0; */
  /* if (fabs(cimag(r0)) < 1e-10) ++nr; */
  /* if (fabs(cimag(r1)) < 1e-10) ++nr; */
  /* if (fabs(cimag(r2)) < 1e-10) ++nr; */
  /* if (fabs(cimag(r3)) < 1e-10) ++nr; */

  /* the check for realness of roots is hard to make robust */
  return 4;
}
