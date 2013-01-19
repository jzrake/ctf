
#include <stdio.h>
#include <math.h>
#include "cow.h"
#if (COW_MPI)
#include <mpi.h>
#endif
#define KILOBYTES (1<<10)
#define MEGABYTES (1<<20)
#define PI (4*atan(1))
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
static void take_elem0(double *result, double **args, int **s, void *u)
{
  *result = args[0][0];
}
static void take_mag3(double *result, double **args, int **s, void *u)
{
  double *m = args[0];
  *result = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
}
static void take_sqr3(double *result, double **args, int **s, void *u)
{
  double *m = args[0];
  *result = m[0]*m[0] + m[1]*m[1] + m[2]*m[2];
}
static void magEtrans(double *result, double **args, int **s, void *u)
{
  double *B = args[0];
  *result = (B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) / (8*PI);
}
static void kinEtrans(double *result, double **args, int **s, void *u)
{
  double *rho = args[0];
  double *v = args[1];
  *result = 0.5 * rho[0] * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

cow_dfield *cow_dfield_new2(cow_domain *domain, char *name)
{
  cow_dfield *f = cow_dfield_new();
  cow_dfield_setdomain(f, domain);
  cow_dfield_setname(f, name);
  return f;
}
cow_dfield *cow_vectorfield(cow_domain *domain, char *name)
{
  cow_dfield *f = cow_dfield_new2(domain, name);
  cow_dfield_addmember(f, "0");
  cow_dfield_addmember(f, "1");
  cow_dfield_addmember(f, "2");
  cow_dfield_commit(f);
  return f;
}
cow_dfield *cow_scalarfield(cow_domain *domain, char *name)
{
  cow_dfield *f = cow_dfield_new2(domain, name);
  cow_dfield_addmember(f, "0");
  cow_dfield_commit(f);
  return f;
}

void cow_dfield_transform(cow_dfield *f, cow_dfield **args, int narg,
			  cow_transform op, void *userdata)
{
  cow_dfield_clearargs(f);
  for (int n=0; n<narg; ++n) {
    cow_dfield_pusharg(f, args[n]);
  }
  cow_dfield_settransform(f, op);
  cow_dfield_setuserdata(f, userdata);
  cow_dfield_transformexecute(f);
}
void cow_dfield_reduce2(cow_dfield *f, cow_transform op, double reduc[3])
{
  cow_dfield_clearargs(f);
  cow_dfield_settransform(f, op);
  cow_dfield_setuserdata(f, NULL);
  cow_dfield_reduce(f, reduc);
}


void make_hist(cow_dfield *f, cow_transform op, char *fout, char *m)
{
  char nickname[1024];
  snprintf(nickname, 1024, "%s-hist", m ? m : cow_dfield_getname(f));

  double reduc[3]; // min, max, sum
  cow_dfield_reduce2(f, op, reduc);
  printf("max, min on %s = %e, %e\n", nickname, reduc[0], reduc[1]);

  cow_histogram *hist = cow_histogram_new();
  cow_histogram_setlower(hist, 0, reduc[0]);
  cow_histogram_setupper(hist, 0, reduc[1]);
  cow_histogram_setnbins(hist, 0, 500);
  cow_histogram_setbinmode(hist, COW_HIST_BINMODE_COUNTS);
  cow_histogram_setdomaincomm(hist, cow_dfield_getdomain(f));
  cow_histogram_commit(hist);
  cow_histogram_setnickname(hist, nickname);
  cow_histogram_populate(hist, f, op);
  cow_histogram_seal(hist);
  cow_histogram_dumphdf5(hist, fout, "");
  cow_histogram_del(hist);
}

int main(int argc, char **argv)
{

  int modes = 0;
  int collective = GETENVINT("COW_HDF5_COLLECTIVE", 0);
  int chunk = GETENVINT("COW_HDF5_CHUNK", 1);
  modes |= GETENVINT("COW_NOREOPEN_STDOUT", 0) ? COW_NOREOPEN_STDOUT : 0;
  modes |= GETENVINT("COW_DISABLE_MPI", 0) ? COW_DISABLE_MPI : 0;

  cow_init(argc, argv, modes);
  if (argc == 3) {
    printf("running on input file %s\n", argv[1]);
  }
  else {
    printf("usage: $> mhdstats infile.h5 outfile.h5\n");
    cow_finalize();
    return 0;
  }
  printf("COW_HDF5_COLLECTIVE: %d\n", collective);

  char *finp = argv[1];
  char *fout = argv[2];
  int derivfields = 0;
  int energies = 1;

  cow_domain *domain = cow_domain_new();
  cow_domain_readsize(domain, finp, "prim/vx");
  cow_domain_setguard(domain, 2);
  cow_domain_commit(domain);

  cow_domain_setchunk(domain, chunk);
  cow_domain_setcollective(domain, collective);
  cow_domain_setalign(domain, 4*KILOBYTES, 4*MEGABYTES);

  cow_dfield *vel = cow_dfield_new2(domain, "prim");
  cow_dfield *mag = cow_dfield_new2(domain, "prim");
  cow_dfield *rho = cow_dfield_new2(domain, "prim");
  cow_dfield *pre = cow_dfield_new2(domain, "prim");
  cow_dfield_addmember(vel, "vx");
  cow_dfield_addmember(vel, "vy");
  cow_dfield_addmember(vel, "vz");
  cow_dfield_addmember(mag, "Bx");
  cow_dfield_addmember(mag, "By");
  cow_dfield_addmember(mag, "Bz");
  cow_dfield_addmember(rho, "rho");
  cow_dfield_addmember(pre, "pre"); // really the internal energy here
  cow_dfield_commit(vel);
  cow_dfield_commit(mag);
  cow_dfield_commit(rho);
  cow_dfield_commit(pre);
  cow_dfield_read(vel, finp);
  cow_dfield_read(mag, finp);
  cow_dfield_read(rho, finp);
  cow_dfield_read(pre, finp);

  if (derivfields) {
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
    cow_dfield_transform(divV, &vel, 1, div5, NULL);
    cow_dfield_transform(curlB, &mag, 1, curl, NULL);
    cow_dfield_transform(curlV, &vel, 1, curl, NULL);

    struct cow_dfield *vcrossBargs[2] = { vel, mag };
    struct cow_dfield *vcrossBcrossBargs[2] = { vel, vcrossB };
    struct cow_dfield *curlBdotBargs[2] = { curlB, mag };
    struct cow_dfield *curlBdotvcrossBargs[2] = { curlB, vcrossB };

    cow_dfield_transform(vcrossB, vcrossBargs, 2, crossprod, NULL);
    cow_dfield_transform(vcrossBcrossB, vcrossBcrossBargs, 2, crossprod, NULL);
    cow_dfield_transform(curlBdotvcrossB, curlBdotvcrossBargs, 2, dotprod, NULL);
    cow_dfield_transform(curlBdotB, curlBdotBargs, 2, dotprod, NULL);
    cow_dfield_transform(divvcrossBcrossB, &vcrossBcrossB, 1, div5, NULL);

    make_hist(divB, take_elem0, fout, NULL);
    make_hist(divV, take_elem0, fout, NULL);
    make_hist(mag, take_sqr3, fout, "B2");
    make_hist(curlB, take_mag3, fout, NULL);
    make_hist(curlV, take_mag3, fout, NULL);
    make_hist(curlBdotvcrossB, take_elem0, fout, NULL);
    make_hist(curlBdotB, take_elem0, fout, NULL);
    make_hist(divvcrossBcrossB, take_elem0, fout, NULL);

    cow_dfield_del(divB);
    cow_dfield_del(divV);
    cow_dfield_del(curlB);
    cow_dfield_del(curlV);
    cow_dfield_del(vcrossB);
    cow_dfield_del(curlBdotvcrossB);
    cow_dfield_del(curlBdotB);
    cow_dfield_del(vcrossBcrossB);
    cow_dfield_del(divvcrossBcrossB);
  }

  if (energies) {
    cow_dfield *kinE = cow_scalarfield(domain, "kinE");
    cow_dfield *magE = cow_scalarfield(domain, "magE");
    cow_dfield *intE = cow_scalarfield(domain, "intE");

    struct cow_dfield *kinEargs[2] = { rho, vel };
    cow_dfield_transform(kinE, kinEargs, 2, kinEtrans, NULL);
    cow_dfield_transform(magE, &mag, 1, magEtrans, NULL);
    cow_dfield_transform(intE, &pre, 1, take_elem0, NULL);

    make_hist(kinE, take_elem0, fout, NULL);
    make_hist(magE, take_elem0, fout, NULL);
    make_hist(intE, take_elem0, fout, NULL);

    cow_dfield_del(kinE);
    cow_dfield_del(magE);
    cow_dfield_del(intE);
  }

  cow_dfield_del(rho);
  cow_dfield_del(pre);
  cow_dfield_del(vel);
  cow_dfield_del(mag);
  cow_domain_del(domain);

  cow_finalize();
  return 0;
}
