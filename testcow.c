
#include <stdio.h>
#include "cow.h"
#if (COW_MPI)
#include <mpi.h>
#endif

void stencildiv(double *result, double **args, int **s, void *u)
{
  double *x = args[0];
  double *y = result;
  int si = s[0][0];
  y[0] = x[+si] - x[-si];
}
void pickmember1(double *result, double **args, int **s, void *u)
{
  *result = args[0][0];
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
  cow_dfield *prim = cow_dfield_new2(domain, "primitive");
  cow_dfield *magf = cow_dfield_new2(domain, "magnetic");

  cow_domain_setndim(domain, 1);
  cow_domain_setguard(domain, 3);
  cow_domain_setsize(domain, 0, 10);
  cow_domain_commit(domain);

  cow_dfield_addmember(prim, "vx");
  cow_dfield_addmember(prim, "vy");
  cow_dfield_addmember(prim, "vz");
  cow_dfield_commit(prim);

  cow_dfield_addmember(magf, "Bx");
  cow_dfield_addmember(magf, "By");
  cow_dfield_addmember(magf, "Bz");
  cow_dfield_commit(magf);

  printf("%s\n", cow_dfield_getname(prim));
  for (const char *m = cow_dfield_iteratemembers(magf);
       m != NULL; m = cow_dfield_nextmember(magf)) {
    printf("\t%s\n", m);
  }
  printf("%s\n", cow_dfield_getname(magf));
  for (const char *m = cow_dfield_iteratemembers(magf);
       m != NULL; m = cow_dfield_nextmember(magf)) {
    printf("\t%s\n", m);
  }

  int si = cow_dfield_getstride(prim, 0);
  int ng = cow_domain_getguard(domain);

  double *P = (double*) cow_dfield_getbuffer(prim);
  double *B = (double*) cow_dfield_getbuffer(magf);
  for (int i=ng; i<cow_domain_getnumlocalzonesinterior(domain, 0)+ng; ++i) {
    P[si*i + 0] = 1.0;
    P[si*i + 1] = 2.0;
    P[si*i + 2] = 3.0;
    B[si*i + 0] = 1.0;
    B[si*i + 1] = 2.0;
    B[si*i + 2] = 3.0;
    printf("(%02d) %f %f %f\n", i, P[si*i + 0], P[si*i + 1], P[si*i + 2]);
  }

  cow_dfield_syncguard(prim);
  printf("\n\n ------------------------------------------------- \n\n");
  for (int i=0; i<cow_domain_getnumlocalzonesinterior(domain, 0)+2*ng; ++i) {
    printf("(%02d) %f %f %f\n", i, P[si*i + 0], P[si*i + 1], P[si*i + 2]);
  }

  cow_dfield *divB = cow_dfield_new2(domain, "divB");
  cow_dfield_addmember(divB, "divB");
  cow_dfield_commit(divB);
  cow_dfield_syncguard(magf);
  cow_dfield_syncguard(divB);

  cow_dfield_clearargs(divB);
  cow_dfield_pusharg(divB, magf);
  cow_dfield_settransform(divB, stencildiv);
  cow_dfield_setuserdata(divB, NULL);
  cow_dfield_transformexecute(divB);

  int I0[] = { 3, 0, 0 };
  int I1[] = { 5, 0, 0 };
  double *subarray = (double*) malloc(2 * 3 * sizeof(double));
  cow_dfield_extract(prim, I0, I1, subarray);
  printf("%f %f\n", subarray[0], subarray[3]);
  printf("%f %f\n", subarray[1], subarray[4]);
  printf("%f %f\n", subarray[2], subarray[5]);
  if (cow_domain_getglobalstartindex(domain, 0) == 0) {
    subarray[0] = 10.0;
    cow_dfield_replace(prim, I0, I1, subarray);
  }
  printf("%f %f\n", subarray[0], subarray[1]);
  free(subarray);
  cow_dfield *divB_copy = cow_dfield_dup(divB);
  cow_dfield_setname(divB_copy, "divB_copy");

  double reduction[3];
  cow_dfield_clearargs(prim);
  cow_dfield_settransform(prim, cow_trans_component);
  cow_dfield_setuserdata(prim, prim);
  cow_dfield_setiparam(prim, 0);
  cow_dfield_reduce(prim, reduction);
  printf("(min max sum): %f %f %f\n", reduction[0], reduction[1], reduction[2]);

  cow_domain_setchunk(domain, 1);
  cow_domain_setcollective(domain, 0);
  cow_domain_setalign(domain, 4096, 4*1024*1024);


  cow_dfield_write(divB_copy, "thefile.h5");
  cow_dfield_write(divB, "thefile.h5");
  cow_dfield_write(magf, "thefile.h5");
  cow_dfield_write(prim, "thefile.h5");

  cow_dfield_read(magf, "thefile.h5");

  cow_dfield_del(divB_copy);
  cow_dfield_del(divB);
  cow_dfield_del(prim);
  cow_dfield_del(magf);
  cow_domain_del(domain);
#if (COW_MPI)
  MPI_Finalize();
#endif
  return 0;
}
