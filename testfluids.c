
#include <stdio.h>
#include "fluids.h"

int main()
{
  fluid_state *S = fluids_new();
  double x[5] = {1, 1, 1, 1, 1};
  double y[5];

  fluids_setfluid(S, FLUIDS_NRHYD);
  fluids_setprimitive(S, x);

  fluids_getprimitive(S, y);
  printf("%f %f %f %f %f\n", y[0], y[1], y[2], y[3], y[4]);

  fluids_del(S);

  return 0;
}

