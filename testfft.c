
#include <stdio.h>
#include "cow.h"
#if (COW_MPI)
#include <mpi.h>
#endif
#define KILOBYTES (1<<10)
#define MEGABYTES (1<<20)
#define GETENVINT(a,dflt) (getenv(a) ? atoi(getenv(a)) : dflt)
#define GETENVDBL(a,dflt) (getenv(a) ? atof(getenv(a)) : dflt)

void cow_fft_pspecvecfield2(cow_dfield *vel, const char *fout, const char *gname)
{
  cow_histogram *hist = cow_histogram_new();
  cow_histogram_setnbins(hist, 0, 200);
  cow_histogram_setspacing(hist, COW_HIST_SPACING_LOG);
  cow_histogram_setnickname(hist, gname);
  cow_fft_pspecvecfield(vel, hist);
  cow_histogram_dumphdf5(hist, fout, gname);
  cow_histogram_del(hist);
}

cow_dfield *cow_dfield_new2(cow_domain *domain, const char *name)
{
  cow_dfield *f = cow_dfield_new();
  cow_dfield_setdomain(f, domain);
  cow_dfield_setname(f, name);
  return f;
}

int main(int argc, char **argv)
{
  int modes = 0;
  int collective = GETENVINT("COW_HDF5_COLLECTIVE", 0);
  modes |= GETENVINT("COW_NOREOPEN_STDOUT", 0) ? COW_NOREOPEN_STDOUT : 0;
  modes |= GETENVINT("COW_DISABLE_MPI", 0) ? COW_DISABLE_MPI : 0;

  cow_init(argc, argv, modes);
  if (argc == 3) {
    printf("running on input file %s\n", argv[1]);
  }
  else {
    printf("usage: $> srhdhist infile.h5 outfile.h5\n");
    cow_finalize();
    return 0;
  }
  printf("COW_HDF5_COLLECTIVE: %d\n", collective);

  char *finp = argv[1];
  char *fout = argv[2];
  cow_domain *domain = cow_domain_new();
  cow_dfield *vel = cow_dfield_new2(domain, "prim");

  cow_domain_setndim(domain, 3);
  cow_domain_setguard(domain, 2);
  cow_domain_setsize(domain, 0, 16);
  cow_domain_setsize(domain, 1, 16);
  cow_domain_setsize(domain, 2, 16);
  cow_domain_commit(domain);

  cow_domain_setchunk(domain, 1);
  cow_domain_setcollective(domain, collective);
  cow_domain_setalign(domain, 4*KILOBYTES, 4*MEGABYTES);

  cow_dfield_addmember(vel, "vx");
  cow_dfield_addmember(vel, "vy");
  cow_dfield_addmember(vel, "vz");
  cow_dfield_commit(vel);

  cow_dfield_read(vel, finp);
  cow_fft_pspecvecfield2(vel, fout, "pspec");

  cow_dfield *dil = cow_dfield_dup(vel);
  cow_dfield *sol = cow_dfield_dup(vel);
  cow_dfield_setname(dil, "dil");
  cow_dfield_setname(sol, "sol");

  cow_fft_helmholtzdecomp(sol, COW_PROJECT_OUT_DIV);
  cow_fft_helmholtzdecomp(dil, COW_PROJECT_OUT_CURL);
  cow_dfield_write(dil, fout);
  cow_dfield_write(sol, fout);
  cow_dfield_write(vel, fout);

  cow_dfield_del(dil);
  cow_dfield_del(sol);
  cow_dfield_del(vel);
  cow_domain_del(domain);

  cow_finalize();
  return 0;
}
