#include "cow.h"
#include "srmhdpack.h"

void srmhdpack_sample2dslice(cow_dfield *prim, int axis, int index, double *P)
{
  double *samp_result;
  cow_domain *domain = cow_dfield_getdomain(prim);
  int i1,i2,q,n=0;
  int a0 = (axis + 0) % 3;
  int a1 = (axis + 1) % 3;
  int a2 = (axis + 2) % 3;
  int N[3];
  int I[3];
  N[0] = cow_domain_getnumglobalzones(domain, a0);
  N[1] = cow_domain_getnumglobalzones(domain, a1);
  N[2] = cow_domain_getnumglobalzones(domain, a2);
  for (i1=0; i1<N[1]; ++i1) {
    for (i2=0; i2<N[2]; ++i2) {
      I[a0] = index;
      I[a1] = i1;
      I[a2] = i2;
      cow_dfield_sampleglobalind(prim, I[0], I[1], I[2], &samp_result, NULL);
      for (q=0; q<8; ++q) {
	P[n++] = samp_result[q];
      }
    }
  }
}
