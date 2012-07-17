
#include <stdio.h>
#include "cow.h"
#if (COW_MPI)
#include <mpi.h>
#endif
#define KILOBYTES (1<<10)
#define MEGABYTES (1<<20)
#define GETENVINT(a,dflt) (getenv(a) ? atoi(getenv(a)) : dflt)
#define GETENVDBL(a,dflt) (getenv(a) ? atof(getenv(a)) : dflt)

static void divcorner(double *result, double **args, int **s, void *u)
{
#define M(i,j,k) ((i)*s[0][0] + (j)*s[0][1] + (k)*s[0][2])
  double *fx = &args[0][0];
  double *fy = &args[0][1];
  double *fz = &args[0][2];
  *result = ((fx[M(1,0,0)] + fx[M(1,1,0)] + fx[M(1,0,1)] + fx[M(1,1,1)]) -
             (fx[M(0,0,0)] + fx[M(0,1,0)] + fx[M(0,0,1)] + fx[M(0,1,1)])) / 4.0
    +       ((fy[M(0,1,0)] + fy[M(0,1,1)] + fy[M(1,1,0)] + fy[M(1,1,1)]) -
             (fy[M(0,0,0)] + fy[M(0,0,1)] + fy[M(1,0,0)] + fy[M(1,0,1)])) / 4.0
    +       ((fz[M(0,0,1)] + fz[M(1,0,1)] + fz[M(0,1,1)] + fz[M(1,1,1)]) -
             (fz[M(0,0,0)] + fz[M(1,0,0)] + fz[M(0,1,0)] + fz[M(1,1,0)])) / 4.0;
#undef M
}
static void div5(double *result, double **args, int **s, void *u)
{
#define diff5(f,s) ((-f[2*s] + 8*f[s] - 8*f[-s] + f[-2*s]) / 12.0)
  double *f0 = &args[0][0];
  double *f1 = &args[0][1];
  double *f2 = &args[0][2];
  *result = diff5(f0, s[0][0]) + diff5(f1, s[0][1]) + diff5(f2, s[0][2]);
#undef diff5
}
static void curl(double *result, double **args, int **s, void *u)
{
  // http://en.wikipedia.org/wiki/Five-point_stencil
#define diff5(f,s) ((-f[2*s] + 8*f[s] - 8*f[-s] + f[-2*s]) / 12.0)
  double *f0 = &args[0][0];
  double *f1 = &args[0][1];
  double *f2 = &args[0][2];
  result[0] = diff5(f2, s[0][1]) - diff5(f1, s[0][2]);
  result[1] = diff5(f0, s[0][2]) - diff5(f2, s[0][0]);
  result[2] = diff5(f1, s[0][0]) - diff5(f0, s[0][1]);
#undef diff5
}
static void crossprod(double *result, double **args, int **s, void *u)
{
  double a0 = args[0][0];
  double a1 = args[0][1];
  double a2 = args[0][2];
  double b0 = args[1][0];
  double b1 = args[1][1];
  double b2 = args[1][2];
  result[0] = a1 * b2 - a2 * b1;
  result[1] = a2 * b0 - a0 * b2;
  result[2] = a0 * b1 - a1 * b0;
}
static void dotprod(double *result, double **args, int **s, void *u)
{
  double a0 = args[0][0];
  double a1 = args[0][1];
  double a2 = args[0][2];
  double b0 = args[1][0];
  double b1 = args[1][1];
  double b2 = args[1][2];
  *result = a0 * b0 + a1 * b1 + a2 * b2;
}
cow_dfield *cow_vectorfield(cow_domain *domain, const char *name)
{
  cow_dfield *f = cow_dfield_new(domain, name);
  cow_dfield_addmember(f, "fx");
  cow_dfield_addmember(f, "fy");
  cow_dfield_addmember(f, "fz");
  cow_dfield_commit(f);
  return f;
}
cow_dfield *cow_scalarfield(cow_domain *domain, const char *name)
{
  cow_dfield *f = cow_dfield_new(domain, name);
  cow_dfield_addmember(f, "f");
  cow_dfield_commit(f);
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

  if (argc == 3) {
    printf("running on input file %s\n", argv[1]);
  }
  else {
    printf("usage: $> milos infile.h5 outfile.h5\n");
    goto done;
  }
  char *finp = argv[1];
  char *fout = argv[2];
  int collective = GETENVINT("COW_HDF5_COLLECTIVE", 0);
  printf("COW_HDF5_COLLECTIVE: %d\n", collective);

  cow_domain *domain = cow_domain_new();
  cow_domain_readsize(domain, finp, "prim/vx");
  cow_domain_setguard(domain, 2);
  cow_domain_commit(domain);
  cow_domain_setchunk(domain, 1);
  cow_domain_setcollective(domain, collective);
  cow_domain_setalign(domain, 4*KILOBYTES, 4*MEGABYTES);

  cow_dfield *vel = cow_dfield_new(domain, "prim");
  cow_dfield *mag = cow_dfield_new(domain, "prim");
  cow_dfield_addmember(vel, "vx");
  cow_dfield_addmember(vel, "vy");
  cow_dfield_addmember(vel, "vz");
  cow_dfield_addmember(mag, "Bx");
  cow_dfield_addmember(mag, "By");
  cow_dfield_addmember(mag, "Bz");
  cow_dfield_commit(vel);
  cow_dfield_commit(mag);
  cow_dfield_read(vel, finp);
  cow_dfield_read(mag, finp);

  cow_dfield *divB = cow_scalarfield(domain, "divB");
  cow_dfield *divV = cow_scalarfield(domain, "divV");
  cow_dfield *curlB = cow_vectorfield(domain, "curlB");
  cow_dfield *curlV = cow_vectorfield(domain, "curlV");
  cow_dfield *vcrossB = cow_vectorfield(domain, "vcrossB");
  cow_dfield *curlBdotvcrossB = cow_scalarfield(domain, "curlBdotvcrossB");
  cow_dfield *curlBdotB = cow_scalarfield(domain, "curlBdotB");
  cow_dfield *vcrossBcrossB = cow_vectorfield(domain, "vcrossBcrossB");
  cow_dfield *divvcrossBcrossB = cow_scalarfield(domain, "divvcrossBcrossB");

  cow_dfield_transform(divB, &mag, 1, divcorner, NULL);
  cow_dfield_transform(divV, &vel, 1, divcorner, NULL);
  cow_dfield_transform(curlB, &mag, 1, curl, NULL);
  cow_dfield_transform(curlV, &vel, 1, curl, NULL);

  cow_dfield_write(divB, fout);
  cow_dfield_write(divV, fout);
  cow_dfield_write(curlB, fout);
  cow_dfield_write(curlV, fout);


  struct cow_dfield *vcrossBargs[2] = { vel, mag };
  struct cow_dfield *vcrossBcrossBargs[2] = { vel, vcrossB };
  struct cow_dfield *curlBdotBargs[2] = { curlB, mag };
  struct cow_dfield *curlBdotvcrossBargs[2] = { curlB, vcrossB };

  cow_dfield_transform(vcrossB, vcrossBargs, 2, crossprod, NULL);
  cow_dfield_transform(vcrossBcrossB, vcrossBcrossBargs, 2, crossprod, NULL);

  cow_dfield_transform(curlBdotvcrossB, curlBdotvcrossBargs, 2, dotprod, NULL);
  cow_dfield_transform(curlBdotB, curlBdotBargs, 2, dotprod, NULL);
  cow_dfield_transform(divvcrossBcrossB, &vcrossBcrossB, 1, div5, NULL);

  cow_dfield_write(curlBdotvcrossB, fout);
  cow_dfield_write(curlBdotB, fout);
  cow_dfield_write(divvcrossBcrossB, fout);

  cow_dfield_del(vel);
  cow_dfield_del(mag);
  cow_dfield_del(divB);
  cow_dfield_del(divV);
  cow_dfield_del(curlB);
  cow_dfield_del(curlV);
  cow_dfield_del(vcrossB);
  cow_dfield_del(curlBdotvcrossB);
  cow_dfield_del(curlBdotB);
  cow_dfield_del(vcrossBcrossB);
  cow_dfield_del(divvcrossBcrossB);
  cow_domain_del(domain);

 done:
#if (COW_MPI)
  MPI_Finalize();
#endif
  return 0;
}
