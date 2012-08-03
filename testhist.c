
#include <stdio.h>
#include "cow.h"
#if (COW_MPI)
#include <mpi.h>
#endif

static void histcb(double *result, double **args, int **s, void *u)
{
  printf("%f %f %f\n", args[0][0], args[0][1], args[0][2]);
  result[0] = args[0][0];
}

cow_dfield *cow_dfield_new2(cow_domain *domain, char *name)
{
  cow_dfield *f = cow_dfield_new();
  cow_dfield_setdomain(f, domain);
  cow_dfield_setname(f, name);
  return f;
}

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
  cow_dfield *data = cow_dfield_new2(domain, "data");

  cow_domain_setndim(domain, 3);
  cow_domain_setguard(domain, 2);
  cow_domain_setsize(domain, 0, 4);
  cow_domain_setsize(domain, 1, 4);
  cow_domain_setsize(domain, 2, 4);
  cow_domain_commit(domain);

  cow_dfield_addmember(data, "d1");
  cow_dfield_addmember(data, "d2");
  cow_dfield_addmember(data, "d3");
  cow_dfield_commit(data);

  double *A = (double*) cow_dfield_getdatabuffer(data);
  for (int i=0; i<cow_domain_getnumlocalzonesincguard(domain, COW_ALL_DIMS);
       ++i) {
    A[3*i + 0] = 0.1;
    A[3*i + 1] = 0.2;
    A[3*i + 2] = 0.3;
  }
  cow_dfield_syncguard(data);

  cow_histogram *hist = cow_histogram_new();
  cow_histogram_setlower(hist, 0, -1.0);
  cow_histogram_setupper(hist, 0, +1.0);
  cow_histogram_setnbins(hist, 0, 200);
  cow_histogram_commit(hist);
  cow_histogram_setnickname(hist, "myhist");

  cow_histogram_populate(hist, data, histcb);

  for (int n=0; n<10000; ++n) {
    double samp = 2.0 * ((double) rand() / RAND_MAX - 0.5);
    double weight = (double) rand() / RAND_MAX;
    cow_histogram_addsample1(hist, samp, weight);
  }
  // test writing it to an ASCII file
  cow_histogram_dumpascii(hist, "thehist.dat");
  // test writing to HDF5 file
  cow_histogram_dumphdf5(hist, "thehist.h5", "");
  // test writing to arbitrary group location
  cow_histogram_dumphdf5(hist, "thehist.h5", "/G1/G2/G3");
  cow_histogram_del(hist);

  cow_dfield_del(data);
  cow_domain_del(domain);

#if (COW_MPI)
  MPI_Finalize();
#endif
  return 0;
}
