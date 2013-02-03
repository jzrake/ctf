

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define FISH_PRIVATE_DEFS
#include "fish.h"


fish_block *fish_block_new()
{
  fish_block *B = (fish_block*) malloc(sizeof(fish_block));
  fish_block block = {
    .rank = 1,
    .size = { 0, 0, 0 },
    .fluid = NULL,
    .error = NULL,
  } ;
  *B = block;
  return B;
}

int fish_block_del(fish_block *B)
{
  free(B);
  return 0;
}

char *fish_block_geterror(fish_block *B)
{
  return B->error;
}

int fish_block_getsize(fish_block *B, int dim)
{
  if (dim < B->rank) {
    return B->size[dim];
  }
  else {
    B->error = "argument 'dim' must be smaller than the rank of the block";
    return FISH_ERROR_BADARG;
  }
}

int fish_block_setsize(fish_block *B, int dim, int size)
{
  if (dim < B->rank) {
    B->size[dim] = size;
    return 0;
  }
  else {
    B->error = "argument 'dim' must be smaller than the rank of the block";
    return FISH_ERROR_BADARG;
  }
}

int fish_block_getrank(fish_block *B)
{
  return B->rank;
}

int fish_block_setrank(fish_block *B, int rank)
{
  if (rank >= 1 && rank <= 3) {
    B->rank = rank;
    return 0;
  }
  else {
    B->error = "rank must be 1, 2, or 3";
    return FISH_ERROR_BADARG;
  }
}
