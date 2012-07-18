
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


void relative_lorentz_factor(cow_dfield *vel, cow_histogram *hist, int N,
                             char mode);

static void div5(double *result, double **args, int **s, void *u)
{
#define diff5(f,s) ((-f[2*s] + 8*f[s] - 8*f[-s] + f[-2*s]) / 12.0)
  double *f0 = &args[0][0];
  double *f1 = &args[0][1];
  double *f2 = &args[0][2];
  *result = diff5(f0, s[0][0]) + diff5(f1, s[0][1]) + diff5(f2, s[0][2]);
#undef diff5
}
static void rot5(double *result, double **args, int **s, void *u)
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
static void take_elm0(double *result, double **args, int **s, void *u)
{
  *result = args[0][0];
}
static void take_mag3(double *result, double **args, int **s, void *u)
{
  double *m = args[0];
  *result = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
}
static void take_lorentzfactor(double *result, double **args, int **s, void *u)
{
  double *m = args[0];
  *result = 1.0 / sqrt(1.0 - (m[0]*m[0] + m[1]*m[1] + m[2]*m[2]));
}

cow_dfield *cow_vectorfield(cow_domain *domain, const char *name)
{
  cow_dfield *f = cow_dfield_new(domain, name);
  cow_dfield_addmember(f, "0");
  cow_dfield_addmember(f, "1");
  cow_dfield_addmember(f, "2");
  cow_dfield_commit(f);
  return f;
}
cow_dfield *cow_scalarfield(cow_domain *domain, const char *name)
{
  cow_dfield *f = cow_dfield_new(domain, name);
  cow_dfield_addmember(f, "0");
  cow_dfield_commit(f);
  return f;
}
void cow_dfield_reduce2(cow_dfield *f, cow_transform op, double reduc[3])
{
  cow_dfield_clearargs(f);
  cow_dfield_settransform(f, op);
  cow_dfield_setuserdata(f, NULL);
  cow_dfield_reduce(f, reduc);
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

void make_hist(cow_dfield *f, cow_transform op, const char *fout, const char *m)
{
  char nickname[1024];
  snprintf(nickname, 1024, "%s-hist", m ? m : cow_dfield_getname(f));

  double reduc[3]; // min, max, sum
  cow_dfield_reduce2(f, op, reduc);

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
  cow_domain_readsize(domain, finp, "prim/vx");
  cow_domain_setguard(domain, 2);
  cow_domain_commit(domain);
  cow_domain_setchunk(domain, 1);
  cow_domain_setcollective(domain, collective);
  cow_domain_setalign(domain, 4*KILOBYTES, 4*MEGABYTES);

  srand(cow_domain_getcartrank(domain));

  cow_dfield *vel = cow_dfield_new(domain, "prim");
  cow_dfield *rho = cow_dfield_new(domain, "prim");
  cow_dfield_addmember(vel, "vx");
  cow_dfield_addmember(vel, "vy");
  cow_dfield_addmember(vel, "vz");
  cow_dfield_addmember(rho, "rho");
  cow_dfield_commit(vel);
  cow_dfield_commit(rho);
  cow_dfield_read(vel, finp);
  cow_dfield_read(rho, finp);

  int dohist = 0;
  int dopair = 1;

  if (dohist) {
    make_hist(vel, take_lorentzfactor, fout, "gamma");
    make_hist(rho, take_elm0, fout, "rho");

    cow_dfield *divV = cow_scalarfield(domain, "divV");
    cow_dfield *rotV = cow_vectorfield(domain, "rotV");
    cow_dfield_transform(divV, &vel, 1, div5, NULL);
    cow_dfield_transform(rotV, &vel, 1, rot5, NULL);

    make_hist(divV, take_elm0, fout, NULL);
    make_hist(rotV, take_mag3, fout, NULL);
    cow_dfield_del(divV);
    cow_dfield_del(rotV);
  }

  if (dopair) {
    cow_histogram *histpro = cow_histogram_new();
    cow_histogram_setlower(histpro, 0, 0.1);
    cow_histogram_setupper(histpro, 0, 1.5);
    cow_histogram_setnbins(histpro, 0, 500);
    cow_histogram_setbinmode(histpro, COW_HIST_BINMODE_AVERAGE);
    cow_histogram_setnickname(histpro, "gamma-rel-drprop-hist");
    cow_histogram_commit(histpro);

    cow_histogram *histlab = cow_histogram_new();
    cow_histogram_setlower(histlab, 0, 0.1);
    cow_histogram_setupper(histlab, 0, 1.5);
    cow_histogram_setnbins(histlab, 0, 500);
    cow_histogram_setbinmode(histlab, COW_HIST_BINMODE_AVERAGE);
    cow_histogram_setnickname(histlab, "gamma-rel-drlab-hist");
    cow_histogram_commit(histlab);

    int nbatch = 5;
    for (int n=0; n<nbatch; ++n) {
      relative_lorentz_factor(vel, histlab, 10000, 'l');
      relative_lorentz_factor(vel, histpro, 10000, 'p');
      printf("batch=%d\n", n);
    }
    cow_histogram_seal(histpro);
    cow_histogram_seal(histlab);
    cow_histogram_dumphdf5(histpro, fout, "");
    cow_histogram_dumphdf5(histlab, fout, "");
    cow_histogram_del(histlab);
    cow_histogram_del(histpro);
  }

  cow_dfield_del(vel);
  cow_dfield_del(rho);
  cow_domain_del(domain);
  cow_finalize();
  return 0;
}



static void lorentz_boost(double u[4], double x[4], double xp[4])
// -----------------------------------------------------------------------------
// Maps the 4-vector x into xp by the boost formed from the 4-velocity u.
//
// see:
// http://en.wikipedia.org/wiki/Lorentz_transformation#Boost_in_any_direction
// -----------------------------------------------------------------------------
{
  double u2 = u[1]*u[1] + u[2]*u[2] + u[3]*u[3];
  double gm = sqrt(1.0 + u2);
  double L[4][4];
  L[0][0] = gm;
  L[0][1] = L[1][0] = -u[1];
  L[0][2] = L[2][0] = -u[2];
  L[0][3] = L[3][0] = -u[3];
  for (int i=1; i<4; ++i) {
    for (int j=1; j<4; ++j) {
      L[i][j] = (gm - 1) * u[i]*u[j] / u2 + (i == j);
    }
  }
  for (int i=0; i<4; ++i) {
    xp[i] = 0.0;
    for (int j=0; j<4; ++j) {
      xp[i] += L[i][j] * x[j];
    }
  }
}

void relative_lorentz_factor(cow_dfield *vel, cow_histogram *hist, int N,
                             char mode)
{
  double *x = (double*) malloc(N * 3 * sizeof(double));
  double *v;

  for (int n=0; n<N; ++n) {
    x[3*n + 0] = (double) rand() / RAND_MAX;
    x[3*n + 1] = (double) rand() / RAND_MAX;
    x[3*n + 2] = (double) rand() / RAND_MAX;
  }

  cow_dfield_setsamplemode(vel, COW_SAMPLE_LINEAR);
  int err = cow_dfield_setsamplecoords(vel, x, N, 3);
  if (err) {
    printf("error on setsamplecoords: %d\n", err);
  }
  free(x);

  int nout1, nout2;
  cow_dfield_sampleexecute(vel);
  cow_dfield_getsamplecoords(vel, &x, &nout1, NULL);
  cow_dfield_getsampleresult(vel, &v, &nout2, NULL);

  for (int n=0; n<nout1/2; ++n) {
    double *x1 = &x[6*n + 0];
    double *x2 = &x[6*n + 3];
    double *v1 = &v[6*n + 0];
    double *v2 = &v[6*n + 3];
    double g1 = 1.0 / sqrt(1.0 - (v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]));
    double g2 = 1.0 / sqrt(1.0 - (v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]));
    double umu1[4] = { g1, g1*v1[0], g1*v1[1], g1*v1[2] };
    double umu2[4] = { g2, g2*v2[0], g2*v2[1], g2*v2[2] };
    double xmu1[4] = { 0.0, x1[0], x1[1], x1[2] };
    double xmu2[4] = { 0.0, x2[0], x2[1], x2[2] };
    double xmu1p[4];
    double xmu2p[4];
    double xrel[4];
    double urel[4];
    lorentz_boost(umu1, xmu1, xmu1p);
    lorentz_boost(umu1, xmu2, xmu2p);
    lorentz_boost(umu1, umu2, urel);
    xrel[0] = xmu2p[0] - xmu1p[0];
    xrel[1] = xmu2p[1] - xmu1p[1];
    xrel[2] = xmu2p[2] - xmu1p[2];
    xrel[3] = xmu2p[3] - xmu1p[3];
    double dx[3] = { x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2] };
    double drlab = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
    double drprop = sqrt(xrel[1]*xrel[1] + xrel[2]*xrel[2] + xrel[3]*xrel[3]);
    double gammarel = urel[0];

    if (mode == 'p') { // proper
      cow_histogram_addsample1(hist, drprop, gammarel);
    }
    else if (mode == 'l') {
      cow_histogram_addsample1(hist, drlab, gammarel);
    }
  }
}
