
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#define FLUIDS_PRIVATE_DEFS
#include "fluids.h"
#include "matrix.h"

int test1()
{
  int counter = 0;
  int N = 64*64*64;
  fluid_state **states = (fluid_state**) malloc(N * sizeof(fluid_state*));
  printf("sizeof(fluid_state) = %ld\n", sizeof(fluid_state));
  for (int n=0; n<N; ++n) {
    states[n] = fluids_new();
    fluids_setfluid(states[n], FLUIDS_NRHYD);
    fluids_alloc(states[n], FLUIDS_PRIMITIVE|FLUIDS_CONSERVED);
  }
  while (counter++ < 1000000000) {
    double x[5] = {1, 1, 1, 1, 1};
    fluids_setattrib(states[0], x, FLUIDS_PRIMITIVE);
    fluids_setattrib(states[0], x, FLUIDS_CONSERVED);
  }
  for (int n=0; n<N; ++n) {
    fluids_del(states[n]);
  }
  return 0;
}

int main()
{
  test1();
  return 0;
}
