#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_ROWS (1<<16)
#define MAX_COLS (1<<8)
#define PI (atan(1.0)*4.0)
#define NVAR 3

static double OutputRadius = 10;
static double StopDensity = 1e-8;
static double PressureConstant = 1.0;
static double PolytropeIndex = 2.0;
static double LambdaValue = 1.0;
static double FourPiG = 1.0;


static double OutputColumnCur[MAX_COLS];
static char  *OutputColumnDes[MAX_COLS];
static int    OutputColumnKey[MAX_COLS];
static int    OutputColumnNum = 0;


static void process_command_line(int argc, char **argv);
static void report();


struct {
  int help;
  int quiet;
  int noprint_rows;
} opts = {0,0};

enum { OC_Z,    // non-dimensional radial coordinate
       OC_R,    // dimensional radial coordinate
       OC_W,    // W = T'
       OC_T,    // T : T'' + T'/x + T^n = 0
       OC_M,    // enclosed mass
       OC_RHO,  // density
       OC_PRE,  // pressure
       OC_PHI,  // gravitational potential
       OC_GPH,  // grad phi
};

#define CH_CASE(key, descr) do {					\
    if (strcmp(argv[n], #key) == 0) {					\
      OutputColumnKey[OutputColumnNum] = OC_##key;			\
	OutputColumnDes[OutputColumnNum] = descr;			\
	MATCHED = 1;							\
	++OutputColumnNum;						\
    }									\
  } while (0)								\

void process_command_line(int argc, char **argv)
{
  int n;
  for (n=1; n<argc; ++n) {
    if (argv[n][0] == '-' && strlen(&argv[n][1]) == 1) {
      switch (argv[n][1]) {
      case 'h': opts.help = 1; break;
      case 'q': opts.quiet = 1; break;
      case 'r': opts.noprint_rows = 1; break;
      default:
	fprintf(stderr, "ERROR: unknown option %s\n", argv[n]);
	exit(2);
      }
    }
    else if (strpbrk(argv[n], "=")) {

      if (argv[n][0] == 'r') {
	OutputRadius = atof(&argv[n][2]);
      }
      else if (argv[n][0] == 'n') {
	PolytropeIndex = atof(&argv[n][2]);
      }
      else if (argv[n][0] == 's') {
	StopDensity = atof(&argv[n][2]);
      }
      else {
	fprintf(stderr, "ERROR: unknown parameter: %c\n", argv[n][0]);
	exit(2);
      }
    }
    else {
      int MATCHED = 0;
      CH_CASE(Z,    "non-dimensional radial coordinate");
      CH_CASE(R,    "dimensional radial coordinate");
      CH_CASE(W,    "W = T'");
      CH_CASE(T,    "T : T'' + T'/x + T^n = 0");
      CH_CASE(M,    "enclosed mass");
      CH_CASE(RHO,  "density");
      CH_CASE(PRE,  "pressure");
      CH_CASE(PHI,  "gravitational potential");
      CH_CASE(GPH,  "grad phi");
      if (!MATCHED) {
      	fprintf(stderr, "ERROR: unknown column header %s\n", argv[n]);
      	exit(2);
      }
    }
  }

  if (opts.help) {
    printf("usage: polytrope [<options>] [R RHO PRE PHI GPH]\n");
    printf("calculate solutions to stellar polytropes\n\n");
    printf(" -h: output help message\n");
    printf("\n");
    exit(1);
  }
}



static void rungekutta4(double *X, double dt, double t);
static void dXdt(double *X, double z, double *dXdt);
static void integrate();


void rungekutta4(double *X, double dt, double t)
{
  int n;
  double L1[NVAR], X1[NVAR];
  double L2[NVAR], X2[NVAR];
  double L3[NVAR], X3[NVAR];
  double L4[NVAR], X4[NVAR];
  for (n=0; n<NVAR; ++n) { X1[n] = X[n] +         0.0 * dt; }
  dXdt(X1, t + 0.0 * dt, L1);
  for (n=0; n<NVAR; ++n) { X2[n] = X[n] + L1[n] * 0.5 * dt; }
  dXdt(X2, t + 0.5 * dt, L2);
  for (n=0; n<NVAR; ++n) { X3[n] = X[n] + L2[n] * 0.5 * dt; }
  dXdt(X3, t + 0.5 * dt, L3);
  for (n=0; n<NVAR; ++n) { X4[n] = X[n] + L3[n] * 1.0 * dt; }
  dXdt(X4, t + 1.0 * dt, L4);
  for (n=0; n<NVAR; ++n) {
    X[n] += (L1[n] + 2*L2[n] + 2*L3[n] + L4[n]) * dt/6.0;
  }
}

void dXdt(double *X, double z, double *dXdt)
{
  double K      = PressureConstant;
  double n      = PolytropeIndex;
  double Lambda = LambdaValue;

  double Z = z;
  double R = z * sqrt((1 + n) * K * pow(Lambda, 1.0/n - 1.0) / FourPiG);
  double W = X[0];
  double T = X[1];
  double density = Lambda * pow(T, n);

  dXdt[0] = -W / Z - pow(T, PolytropeIndex);
  dXdt[1] =  W;
  dXdt[2] =  2*PI*density*R;
}


void output(double z, double *X)
{
  int j;

  double K      = PressureConstant;
  double n      = PolytropeIndex;
  double Lambda = LambdaValue;
  double Gamma  = 1.0 + 1.0/n;

  double Z = z;
  double R = z * sqrt((1 + n) * K * pow(Lambda, 1.0/n - 1.0) / FourPiG);
  double W = X[0];
  double T = X[1];
  double M = X[2];

  double V = T - 1; // equation 8
  double B = pow(Lambda, 1.0/n) * K * (1 + n) * V; // equation 7
  double B_prime = B * W / V;

  double density  = Lambda * pow(T, n);
  double pressure = K * pow(density, Gamma);
  double potential = -B;
  double gradphi = -B_prime;

  for (j=0; j<OutputColumnNum; ++j) {
    double colval;
    switch (OutputColumnKey[j]) {
    case OC_Z: colval = Z; break;
    case OC_R: colval = R; break;
    case OC_T: colval = T; break;
    case OC_W: colval = W; break;
    case OC_M: colval = M; break;
    case OC_RHO: colval = density; break;
    case OC_PRE: colval = pressure; break;
    case OC_PHI: colval = potential; break;
    case OC_GPH: colval = gradphi; break;
    default:
      fprintf(stderr, "ERROR: column header not implemented\n");
      exit(2);
    }
    OutputColumnCur[j] = colval;
    if (!opts.noprint_rows) printf("%lf ", colval);
  }
  if (!opts.noprint_rows) printf("\n");
}

void integrate()
{
  double K      = PressureConstant;
  double n      = PolytropeIndex;
  double Lambda = LambdaValue;

  double R, W, T, density;
  double Z = 1e-10;
  double dZ = 1e-4;
  double X[3];
  X[0] = 0.0; // W
  X[1] = 1.0; // T
  X[2] = 0.0; // M

  do {
    R = Z * sqrt( (1 + n) * K * pow(Lambda, 1.0/n - 1.0) / FourPiG);
    W = X[0];
    T = X[1];
    density = LambdaValue * pow(T, n);

    if (OutputColumnNum > 0) {
      output(Z, X);
    }

    rungekutta4(X, dZ, Z);
    Z += dZ;

  } while (R < OutputRadius && density > StopDensity);
}


int main(int argc, char **argv)
{
  int i;
  process_command_line(argc, argv);

  if (!opts.quiet) {
    report();
  }

  integrate();

  if (!opts.quiet) {
    printf("FINAL VALUES:\n");
    for (i=0; i<OutputColumnNum; ++i) printf("%lf ", OutputColumnCur[i]);
    printf("\n");
  }

  return 0;
}



void report()
{
  int i;

  printf("r=%f\n", OutputRadius);

  printf("COLUMN OUTPUT:\n");
  for (i=0; i<OutputColumnNum; ++i) {
    printf("col %d: %s\n", i+1, OutputColumnDes[i]);
  }
}
