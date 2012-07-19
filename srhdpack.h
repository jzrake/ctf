

#ifndef SRHDPACK_HEADER_INCLUDED
#define SRHDPACK_HEADER_INCLUDED
#include "cow.h"

void srhdpack_relativelorentzpairs(cow_dfield *vel,
				   cow_histogram *histpro,
				   cow_histogram *histlab,
				   int nbatch,
				   int nperbatch,
				   int seed);

#endif // SRHDPACK_HEADER_INCLUDED
