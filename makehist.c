
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cow.h"
#if (COW_MPI)
#include <mpi.h>
#endif
#define KILOBYTES (1<<10)
#define MEGABYTES (1<<20)
#define GETENVINT(a,dflt) (getenv(a) ? atoi(getenv(a)) : dflt)
#define GETENVDBL(a,dflt) (getenv(a) ? atof(getenv(a)) : dflt)

int ncalls = 0;
void take_elem0(double *result, double **args, int **s, void *u)
{
  ++ncalls;
  *result = args[0][0];
}
void take_mag3(double *result, double **args, int **s, void *u)
{
  ++ncalls;
  double *m = args[0];
  *result = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
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
  char *finp = argc > 1 ? argv[1] : NULL;
  char *grou = argc > 2 ? argv[2] : NULL;
  char *memb = argc > 3 ? argv[3] : NULL;
  char *fout = argc > 4 ? argv[4] : NULL;

  if (finp && grou && memb && fout) {
    printf("opening data set '%s/%s' from input file '%s'\n", grou, memb, finp);
  }
  else {
    printf("usage: $> makehist infile.h5 group member outfile.h5\n");
    goto done;
  }
  char wholedset1[1024];
  char wholedset2[1024];
  snprintf(wholedset1, 1024, "%s/%s", grou, memb);
  snprintf(wholedset2, 1024, "%s-hist", grou);

  int collective = GETENVINT("COW_HDF5_COLLECTIVE", 0);
  printf("COW_HDF5_COLLECTIVE: %d\n", collective);

  cow_domain *domain = cow_domain_new();
  cow_domain_readsize(domain, finp, wholedset1);
  cow_domain_setguard(domain, 2);
  cow_domain_commit(domain);
  cow_domain_setchunk(domain, 1);
  cow_domain_setcollective(domain, collective);
  cow_domain_setalign(domain, 4*KILOBYTES, 4*MEGABYTES);

  cow_dfield *f = cow_dfield_new(domain, grou);
  cow_dfield_addmember(f, memb);
  cow_dfield_commit(f);
  cow_dfield_read(f, finp);

  double reduc[3]; // min, max, sum
  cow_dfield_reduce(f, take_elem0, reduc);

  cow_histogram *hist = cow_histogram_new();
  cow_histogram_setlower(hist, 0, reduc[0]);
  cow_histogram_setupper(hist, 0, reduc[1]);
  cow_histogram_setnbins(hist, 0, 500);
  cow_histogram_setbinmode(hist, COW_HIST_BINMODE_COUNTS);
  cow_histogram_setdomaincomm(hist, domain);
  cow_histogram_commit(hist);
  cow_histogram_setnickname(hist, wholedset2);
  cow_histogram_populate(hist, f, take_elem0);
  cow_histogram_dumphdf5(hist, fout, "");
  cow_histogram_del(hist);

  cow_dfield_del(f);
  cow_domain_del(domain);
  printf("measured %d cells\n", ncalls);

 done:
#if (COW_MPI)
  MPI_Finalize();
#endif
  return 0;
}
