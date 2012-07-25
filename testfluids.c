
#include <stdio.h>
#include "fluids.h"

int main()
{
  fluid_state *S = fluids_new();
  double x[5] = {1, 1, 1, 1, 1};
  double y[5];

  fluids_setfluid(S, FLUIDS_NRHYD);
  fluids_setattrib(S, x, FLUIDS_PRIMITIVE);
  fluids_getattrib(S, y, FLUIDS_PRIMITIVE);

  printf("%f %f %f %f %f\n", y[0], y[1], y[2], y[3], y[4]);

  fluids_del(S);

  return 0;
}

