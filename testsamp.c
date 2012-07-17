
#include <stdio.h>
#include "cow.h"
#if (COW_MPI)
#include <mpi.h>
#endif

int main(int argc, char **argv)
{
#if (COW_MPI)
  {
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) freopen("/dev/null", "w", stdout);
    printf("was compiled with MPI support\n");
  }
#endif

  cow_domain *domain = cow_domain_new();
  cow_dfield *data = cow_dfield_new(domain, "data");

  cow_domain_setndim(domain, 3);
  cow_domain_setguard(domain, 2);
  cow_domain_setsize(domain, 0, 42);
  cow_domain_setsize(domain, 1, 42);
  cow_domain_setsize(domain, 2, 42);
  cow_domain_commit(domain);

  cow_dfield_addmember(data, "d1");
  cow_dfield_addmember(data, "d2");
  cow_dfield_addmember(data, "d3");
  cow_dfield_commit(data);

  double *A = (double*) cow_dfield_getbuffer(data);
  for (int i=0; i<cow_domain_getnumlocalzonesincguard(domain, COW_ALL_DIMS);
       ++i) {
    A[3*i + 0] = 0.1;
    A[3*i + 1] = 0.2;
    A[3*i + 2] = 0.3;
  }
  cow_dfield_syncguard(data);

  int N = 10 * cow_domain_getcartrank(domain);
  double *r0 = (double*) malloc(N * 3 * sizeof(double));
  double *r1 = (double*) malloc(N * 3 * sizeof(double));
  double *sample = (double*) malloc(N * 3 * sizeof(double));

  for (int n=0; n<3*N; ++n) {
    r0[n] = (double) rand() / RAND_MAX;
  }
  char fname[256];
  cow_dfield_sampleglobalpos(data, r0, N, r1, sample, COW_SAMPLE_LINEAR);
  snprintf(fname, 256, "samp-%02d.dat", cow_domain_getcartrank(domain));
  FILE *fout = fopen(fname, "w");
  for (int n=0; n<N; ++n) {
    fprintf(fout, "%f %f %f\n", sample[3*n+0], sample[3*n+1], sample[3*n+2]);
  }
  fclose(fout);
  free(r0);
  free(r1);
  free(sample);

  double P[3];
  cow_dfield_sampleglobalind(data, 12, 12, 12, P);
  printf("%f %f %f\n", P[0], P[1], P[2]);

  cow_dfield_del(data);
  cow_domain_del(domain);

#if (COW_MPI)
  MPI_Finalize();
#endif
  return 0;
}
