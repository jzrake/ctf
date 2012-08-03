
#include <stdio.h>
#include <time.h>
#include "cow.h"
#if (COW_MPI)
#include <mpi.h>
#endif

#define KILOBYTES (1<<10)
#define MEGABYTES (1<<20)
#define GETENVINT(a,dflt) (getenv(a) ? atoi(getenv(a)) : dflt)
#define GETENVDBL(a,dflt) (getenv(a) ? atof(getenv(a)) : dflt)

cow_dfield *cow_dfield_new2(cow_domain *domain, char *name)
{
  cow_dfield *f = cow_dfield_new();
  cow_dfield_setdomain(f, domain);
  cow_dfield_setname(f, name);
  return f;
}

int main(int argc, char **argv)
{
  clock_t start = clock();

  int modes = 0;
  int collective = GETENVINT("COW_HDF5_COLLECTIVE", 0);
  int chunk = GETENVINT("COW_HDF5_CHUNK", 1);
  modes |= GETENVINT("COW_NOREOPEN_STDOUT", 0) ? COW_NOREOPEN_STDOUT : 0;
  modes |= GETENVINT("COW_DISABLE_MPI", 0) ? COW_DISABLE_MPI : 0;

  cow_init(argc, argv, modes);

  cow_domain *domain = cow_domain_new();
  cow_dfield *prim = cow_dfield_new2(domain, "prim");

  cow_domain_setndim(domain, 3);
  cow_domain_setguard(domain, 3);
  cow_domain_setsize(domain, 0, 512);
  cow_domain_setsize(domain, 1, 512);
  cow_domain_setsize(domain, 2, 512);
  cow_domain_commit(domain);

  cow_domain_setchunk(domain, chunk);
  cow_domain_setcollective(domain, collective);
  cow_domain_setalign(domain, 4*KILOBYTES, 4*MEGABYTES);

  cow_dfield_addmember(prim, "rho");
  cow_dfield_addmember(prim, "pre");
  cow_dfield_addmember(prim, "vx");
  cow_dfield_addmember(prim, "vy");
  cow_dfield_addmember(prim, "vz");
  cow_dfield_commit(prim);

  double *P = (double*) cow_dfield_getdatabuffer(prim);

  for (int i=0; i<cow_domain_getnumlocalzonesincguard(domain, COW_ALL_DIMS);
       ++i) {
    P[5*i + 0] = 1.0;
    P[5*i + 1] = 2.0;
    P[5*i + 2] = 3.0;
    P[5*i + 3] = 4.0;
    P[5*i + 4] = 5.0;
  }
  cow_dfield_write(prim, "data/cow-io.h5");
  printf("test finished in %3.2f seconds\n",
	 (double)(clock()-start)/CLOCKS_PER_SEC);

  cow_dfield_del(prim);
  cow_domain_del(domain);
  cow_finalize();
  return 0;
}
