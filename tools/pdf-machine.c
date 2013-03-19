
/* -----------------------------------------------------------------------------
 *
 * PROGRAM: pdf-machine
 *
 * AUTHOR: Jonathan Zrake
 *
 * PURPOSE: Process histogram output
 *
 *
 * tests:
 *
 * gsl-randist 12345 100000 gaussian 1 | gsl-histogram -10 10 1000 | pdf-machine
 *
 *
 * -----------------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_ROWS (1<<16)
#define MAX_COLS (1<<8)
#define PI (atan(1.0)*4.0)

/* Be careful with portability of large integer types */
#include <stdint.h>
#include <inttypes.h>
#define COUNT_FORMAT "%"PRIu64 // guarentee 64 byte integer
typedef uint64_t count_t;


static double BinEdges[MAX_ROWS];
static count_t dN[MAX_ROWS][MAX_COLS];
static int NR = 0;
static int NC = 0;

static double OutputColumnTot[MAX_COLS];
static char  *OutputColumnDes[MAX_COLS];
static int    OutputColumnKey[MAX_COLS];
static int    OutputColumnNum = 0;

static void process_command_line(int argc, char **argv);
static void report();
static void read_input();
static void output(int nr);
static double pdf_mean();
static double pdf_variance();
static count_t distribution_total_counts(int nc);
static count_t distribution_total_counts_all();


/* PDF models */
static double pdf_model_normal(double x);
static double pdf_model_lognormal(double x);
static double pdf_model_hopkins2013(double x);


/* Bessel function evaluations: http://www.astro.rug.nl/~gipsy/sub/bessel.c */
double bessj( int, double );
double bessy( int, double );
double bessi( int, double );
double bessk( int, double );


/* cached values of repetitive function calls */
static double cached_pdf_mean;
static double cached_pdf_variance;
static count_t cached_distribution_total_counts_all;

struct {
  int help;
  int quiet;
  int skip_zeros;
  int noprint_bins;
  double hopkins_T;
} opts = {0,0,0,0,0.1};

enum { OC_E0,  // bin left-edge
       OC_E1,  // bin right-edge
       OC_MP,  // bin mid-point
       OC_MC,  // bin mean counts
       OC_TC,  // bin total counts
       OC_P,   // normalized PDF value
       OC_STD, // bin standard deviation over histograms
       OC_GD,  // PDF model: Gaussian
       OC_LN,  // PDF model: log-normal
       OC_PH,  // PDF model: Hopkins (2013)
       OC_R,   // model residual: |P - model|
       OC_LR,  // model residual, log10: |log10(P) - log10(model)|
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
      case 'z': opts.skip_zeros = 1; break;
      case 'r': opts.noprint_bins = 1; break;
      case 'T': opts.hopkins_T = atof(argv[++n]); break;
      default:
	fprintf(stderr, "ERROR: unknown option %s\n", argv[n]);
	exit(2);
      }
    }
    else {
      int MATCHED = 0;
      CH_CASE(E0,  "bin left-edge");
      CH_CASE(E1,  "bin right-edge");
      CH_CASE(MP,  "bin mid-point");
      CH_CASE(MC,  "bin mean counts");
      CH_CASE(TC,  "bin total counts");
      CH_CASE(P,   "normalized PDF value");
      CH_CASE(STD, "standard deviation of counts in bin");
      CH_CASE(GD,  "PDF model: Gaussian");
      CH_CASE(LN,  "PDF model: log-normal");
      CH_CASE(PH,  "PDF model: Hopkins (2013)");
      CH_CASE(R,   "model residual: |P - model|");
      CH_CASE(LR,  "model residual, log10: |log10(P) - log10(model)|");

      if (!MATCHED) {
      	fprintf(stderr, "ERROR: unknown column header %s\n", argv[n]);
      	exit(2);
      }
    }
  }

  if (opts.help) {
    printf("usage: pdf-machine [<options>] [E0 E1 MP MC TC P STD GD LN PH]\n");
    printf("output model PDF's from an ensemble of histograms\n\n");
    printf(" -h: output help message\n");
    printf(" -q: limit output to raw PDF only, good if redirecting output\n");
    printf(" -z: skip zeros in PDF output\n");
    printf(" -r: skip bin output\n");
    printf(" -T 0.1: parameter for Hopkins (2013) model\n");
    printf("\n");
    printf("  example input:\n\n");
    printf("  0.0 0.1 8 4 2 1\n");
    printf("  0.1 0.2 3 8 3 6\n");
    printf("  0.2 0.3 6 5 8 1\n");
    printf("  0.3 0.4 1 2 8 3\n\n");

    printf("  allows up to %d bin edges and %d simultaneous PDF's\n",
	   MAX_ROWS, MAX_COLS);
    printf("  uses %ld MB of RAM\n", (sizeof(dN) + sizeof(BinEdges)) / (1<<20));
    exit(1);
  }
}

int main(int argc, char **argv)
{
  int i;
  process_command_line(argc, argv);
  read_input();

  cached_pdf_mean = pdf_mean();
  cached_pdf_variance = pdf_variance();
  cached_distribution_total_counts_all = distribution_total_counts_all();

  for (i=0; i<MAX_COLS; ++i) OutputColumnTot[i] = 0.0;

  if (!opts.quiet) {
    report();
  }

  if (OutputColumnNum > 0) {
    for (i=0; i<NR; ++i) {
      output(i);
    }
  }

  if (!opts.quiet) {
    printf("COLUMN INTEGRALS:\n");
    for (i=0; i<OutputColumnNum; ++i) printf("%lf ", OutputColumnTot[i]);
    printf("\n");
  }

  return 0;
}

void report()
{
  int i;
  printf("number of bins ... %d\n", NR);
  printf("number of distributions ... %d\n", NC);

  printf("PDF mean ... %f\n", pdf_mean());
  printf("PDF variance (S) ... %f\n", pdf_variance());
  printf("-S/2 (=M if log-normal and unit mean) ... %f\n", -pdf_variance()/2);
  printf("Hopkins (2013) T parameter ... %f\n", opts.hopkins_T);

  printf("total counts over all distributions ... "COUNT_FORMAT"\n",
	 cached_distribution_total_counts_all);

  printf("COLUMN OUTPUT:\n");
  for (i=0; i<OutputColumnNum; ++i) {
    printf("col %d: %s\n", i+1, OutputColumnDes[i]);
  }
}

count_t distribution_total_counts(int nc)
{
  int i;
  count_t N=0;
  for (i=0; i<NR; ++i) {
    N += dN[i][nc];
  }
  return N;
}

count_t distribution_total_counts_all()
{
  int i;
  count_t c=0;
  for (i=0; i<NC; ++i) c += distribution_total_counts(i);
  return c;
}

double pdf_variance()
{
  int nbin, i, j;
  count_t ntot=0;
  double sum=0.0, x, f;
  double mean = pdf_mean();
  for (i=0; i<NR; ++i) {
    nbin = 0;
    x = 0.5 * (BinEdges[i+1] + BinEdges[i]);
    f = (x - mean) * (x - mean);
    for (j=0; j<NC; ++j) {
      ntot += dN[i][j];
      nbin += dN[i][j];
    }
    sum += f * nbin;
  }
  return sum / ntot;
}

double pdf_mean()
{
  int nbin, i, j;
  count_t ntot=0;
  double sum=0.0, x;
  for (i=0; i<NR; ++i) {
    nbin = 0;
    x = 0.5 * (BinEdges[i+1] + BinEdges[i]);
    for (j=0; j<NC; ++j) {
      ntot += dN[i][j];
      nbin += dN[i][j];
    }
    sum += x * nbin;
  }
  return sum / ntot;
}


void output(int nr)
{
  int j, counts=0;
  count_t N = cached_distribution_total_counts_all;
  double x = 0.5 * (BinEdges[nr+1] + BinEdges[nr]);
  double h = 1.0 * (BinEdges[nr+1] - BinEdges[nr]);
  double v = 0.0; // bin count variance
  double P, last_model_pdf=0.0;

  for (j=0; j<NC; ++j) {
    counts += dN[nr][j];
  }
  if (counts == 0 && opts.skip_zeros) return;

  for (j=0; j<NC; ++j) {
    v += pow(dN[nr][j] - (double) counts / NC, 2.0) / NC;
  }

  P = (double)counts / N / h; // pdf value

  for (j=0; j<OutputColumnNum; ++j) {
    double colval;
    switch (OutputColumnKey[j]) {
    case OC_E0  : colval = BinEdges[nr]; break;
    case OC_E1  : colval = BinEdges[nr+1]; break;
    case OC_MP  : colval = x; break;
    case OC_TC  : colval = (double)counts; break;
    case OC_MC  : colval = (double)counts / NC; break;
    case OC_P   : colval = P; break;
    case OC_STD : colval = sqrt(v); break;
    case OC_GD  : colval = last_model_pdf = pdf_model_normal(x); break;
    case OC_LN  : colval = last_model_pdf = pdf_model_lognormal(x); break;
    case OC_PH  : colval = last_model_pdf = pdf_model_hopkins2013(x); break;
    case OC_R   : colval = fabs(last_model_pdf - P); break;
    case OC_LR  : colval = fabs(log10(last_model_pdf) - log10(P)); break;
    default:
      fprintf(stderr, "ERROR: column header not implemented\n");
      exit(2);
    }
    OutputColumnTot[j] += colval * h;
    if (!opts.noprint_bins) printf("%lf ", colval);
  }
  if (!opts.noprint_bins) printf("\n");
}

double pdf_model_normal(double x)
{
  double S = cached_pdf_variance;
  double M = cached_pdf_mean;
  double P = 1.0 / sqrt(2 * PI * S) * exp(-0.5*pow(x - M, 2.0)/S);
  return P;
}

double pdf_model_lognormal(double x)
{
  double S = cached_pdf_variance;
  double M = -S/2;
  double P = 1.0 / sqrt(2 * PI * S) * exp(-0.5*pow(x - M, 2.0)/S);
  return P;
}

double pdf_model_hopkins2013(double x)
{
  double S = cached_pdf_variance;
  double T = opts.hopkins_T;
  double L = S / (2*T*T);
  double u = L / (1 + T) - x / T;
  double P = u < 0.0 ? 0.0 : bessi(1, 2*sqrt(L*u)) * exp(-(L + u)) * sqrt(L/u) / T;

  if (P != P) {
    fprintf(stderr, "ERROR: Bessel function failed, T is too small\n");
    exit(2);
  }

  return P;
}

void read_input()
{
  char c;
  int nr=0, nc, v;
  double x0, x1;
  while (1) {
    c = fgetc(stdin);
    if (feof(stdin)) {
      break;
    }
    if (c == '\n') { nr++; NC = nc - 1; nc = 0; }
    else {
      while (c == ' ') c = fgetc(stdin);
      ungetc(c, stdin);
      if (nc == 0) {
	if (nr+1 >= MAX_ROWS) {
	  fprintf(stderr, "ERROR: exceeds max number of bins %d\n", MAX_ROWS);
	  exit(2);
	}
	fscanf(stdin, "%lf %lf", &x0, &x1);
	BinEdges[nr  ] = x0;
	BinEdges[nr+1] = x1;
      }
      else {
	if (nc-1 >= MAX_COLS) {
	  fprintf(stderr, "ERROR: exceeds max number of PDF's %d\n", MAX_COLS);
	  exit(2);
	}
	fscanf(stdin, "%d", &v);
	dN[nr][nc-1] = v;
      }
      nc++;
    }
  }
  NR = nr;
}









/* bessel.c
*/

/*
#>            bessel.dc2

Function:     BESSEL

Purpose:      Evaluate Bessel function J, Y, I, K of integer order.

Category:     MATH

File:         bessel.c

Author:       M.G.R. Vogelaar

Use:          See bessj.dc2, bessy.dc2, bessi.dc2 or bessk.dc2
                      
Description:  The differential equation 

                       2
                   2  d w       dw      2   2
                  x . --- + x . --- + (x - v ).w = 0
                        2       dx
                      dx
                      
              has two solutions called Bessel functions of the first kind
              Jv(x) and Bessel functions of the second kind Yv(x).
              The routines bessj and bessy return the J and Y for 
              integer v and therefore are called Bessel functions 
              of integer order.
              
              The differential equation 

                       2
                   2  d w       dw      2   2
                  x . --- + x . --- - (x + v ).w = 0
                        2       dx
                      dx
                      
              has two solutions called modified Bessel functions
              Iv(x) and Kv(x).
              The routines bessi and bessk return the I and K for 
              integer v and therefore are called Modified Bessel 
              functions of integer order.
              (Abramowitz & Stegun, Handbook of mathematical 
              functions, ch. 9, pages 358,- and 374,- )
                            
              The implementation is based on the ideas from 
              Numerical Recipes, Press et. al. 
              This routine is NOT callable in FORTRAN.

Updates:      Jun 29, 1998: VOG, Document created.
#<
*/

/*
#> bessel.h
#if !defined(_bessel_h_) 
#define _bessel_h_
extern double bessj( int, double );
extern double bessy( int, double );
extern double bessi( int, double );
extern double bessk( int, double );
#endif
#<
*/


#define setdblank_c(a) *a=0.0
#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10


static double bessj0( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of first kind and order  */
/*          0 at input x                                      */
/*------------------------------------------------------------*/
{
   double ax,z;
   double xx,y,ans,ans1,ans2;

   if ((ax=fabs(x)) < 8.0) {
      y=x*x;
      ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
         +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
      ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
         +y*(59272.64853+y*(267.8532712+y*1.0))));
      ans=ans1/ans2;
   } else {
      z=8.0/ax;
      y=z*z;
      xx=ax-0.785398164;
      ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
         +y*(-0.2073370639e-5+y*0.2093887211e-6)));
      ans2 = -0.1562499995e-1+y*(0.1430488765e-3
         +y*(-0.6911147651e-5+y*(0.7621095161e-6
         -y*0.934935152e-7)));
      ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
   }
   return ans;
}



static double bessj1( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of first kind and order  */
/*          1 at input x                                      */
/*------------------------------------------------------------*/
{
   double ax,z;
   double xx,y,ans,ans1,ans2;

   if ((ax=fabs(x)) < 8.0) {
      y=x*x;
      ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
         +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
      ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
         +y*(99447.43394+y*(376.9991397+y*1.0))));
      ans=ans1/ans2;
   } else {
      z=8.0/ax;
      y=z*z;
      xx=ax-2.356194491;
      ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
         +y*(0.2457520174e-5+y*(-0.240337019e-6))));
      ans2=0.04687499995+y*(-0.2002690873e-3
         +y*(0.8449199096e-5+y*(-0.88228987e-6
         +y*0.105787412e-6)));
      ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
      if (x < 0.0) ans = -ans;
   }
   return ans;
}



/*
#>            bessj.dc2

Function:     bessj

Purpose:      Evaluate Bessel function of first kind of integer order.

Category:     MATH

File:         bessel.c

Author:       M.G.R. Vogelaar

Use:          #include "bessel.h"
              double   result; 
              result = bessj( int n,
                              double x )


              bessj    Return the Bessel function of integer order
                       for input value x.
              n        Integer order of Bessel function.
              x        Double at which the function is evaluated.

                      
Description:  bessj evaluates at x the Bessel function of the first kind 
              and of integer order n. 
              This routine is NOT callable in FORTRAN.

Updates:      Jun 29, 1998: VOG, Document created.
#<
*/


double bessj( int n, double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of first kind and order  */
/*          n at input x                                      */
/* The function can also be called for n = 0 and n = 1.       */
/*------------------------------------------------------------*/
{
   int    j, jsum, m;
   double ax, bj, bjm, bjp, sum, tox, ans;


   if (n < 0)
   {
      double   dblank;
      setdblank_c( &dblank );
      return( dblank );
   }
   ax=fabs(x);
   if (n == 0)
      return( bessj0(ax) );
   if (n == 1)
      return( bessj1(ax) );
      

   if (ax == 0.0)
      return 0.0;
   else if (ax > (double) n) {
      tox=2.0/ax;
      bjm=bessj0(ax);
      bj=bessj1(ax);
      for (j=1;j<n;j++) {
         bjp=j*tox*bj-bjm;
         bjm=bj;
         bj=bjp;
      }
      ans=bj;
   } else {
      tox=2.0/ax;
      m=2*((n+(int) sqrt(ACC*n))/2);
      jsum=0;
      bjp=ans=sum=0.0;
      bj=1.0;
      for (j=m;j>0;j--) {
         bjm=j*tox*bj-bjp;
         bjp=bj;
         bj=bjm;
         if (fabs(bj) > BIGNO) {
            bj *= BIGNI;
            bjp *= BIGNI;
            ans *= BIGNI;
            sum *= BIGNI;
         }
         if (jsum) sum += bj;
         jsum=!jsum;
         if (j == n) ans=bjp;
      }
      sum=2.0*sum-bj;
      ans /= sum;
   }
   return  x < 0.0 && n%2 == 1 ? -ans : ans;
}




static double bessy0( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of second kind and order */
/*          0 at input x.                                     */
/*------------------------------------------------------------*/
{
   double z;
   double xx,y,ans,ans1,ans2;

   if (x < 8.0) {
      y=x*x;
      ans1 = -2957821389.0+y*(7062834065.0+y*(-512359803.6
         +y*(10879881.29+y*(-86327.92757+y*228.4622733))));
      ans2=40076544269.0+y*(745249964.8+y*(7189466.438
         +y*(47447.26470+y*(226.1030244+y*1.0))));
      ans=(ans1/ans2)+0.636619772*bessj0(x)*log(x);
   } else {
      z=8.0/x;
      y=z*z;
      xx=x-0.785398164;
      ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
         +y*(-0.2073370639e-5+y*0.2093887211e-6)));
      ans2 = -0.1562499995e-1+y*(0.1430488765e-3
         +y*(-0.6911147651e-5+y*(0.7621095161e-6
         +y*(-0.934945152e-7))));
      ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
   }
   return ans;
}



static double bessy1( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of second kind and order */
/*          1 at input x.                                     */
/*------------------------------------------------------------*/
{
   double z;
   double xx,y,ans,ans1,ans2;

   if (x < 8.0) {
      y=x*x;
      ans1=x*(-0.4900604943e13+y*(0.1275274390e13
         +y*(-0.5153438139e11+y*(0.7349264551e9
         +y*(-0.4237922726e7+y*0.8511937935e4)))));
      ans2=0.2499580570e14+y*(0.4244419664e12
         +y*(0.3733650367e10+y*(0.2245904002e8
         +y*(0.1020426050e6+y*(0.3549632885e3+y)))));
      ans=(ans1/ans2)+0.636619772*(bessj1(x)*log(x)-1.0/x);
   } else {
      z=8.0/x;
      y=z*z;
      xx=x-2.356194491;
      ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
         +y*(0.2457520174e-5+y*(-0.240337019e-6))));
      ans2=0.04687499995+y*(-0.2002690873e-3
         +y*(0.8449199096e-5+y*(-0.88228987e-6
         +y*0.105787412e-6)));
      ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
   }
   return ans;
}



/*
#>            bessy.dc2

Function:     bessy

Purpose:      Evaluate Bessel function second kind and of integer order.

Category:     MATH

File:         bessel.c

Author:       M.G.R. Vogelaar

Use:          #include "bessel.h"
              double   result; 
              result = bessy( int n,
                              double x )


              bessy    Return the Bessel function of second kind and
                       of integer order, for input value x.
              n        Integer order of Bessel function.
              x        Double at which the function is evaluated.

                      
Description:  bessy evaluates at x the Bessel function of the second kind
              and of integer order n.
              This routine is NOT callable in FORTRAN.

Updates:      Jun 29, 1998: VOG, Document created.
#<
*/


double bessy( int n, double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of second kind and order */
/*          n for input x. (n >= 0)                           */
/* Note that for x == 0 the functions bessy and bessk are not */
/* defined and a blank is returned.                           */
/*------------------------------------------------------------*/
{
   int j;
   double by,bym,byp,tox;


   if (n < 0 || x == 0.0)
   {
      double   dblank;
      setdblank_c( &dblank );
      return( dblank );
   }
   if (n == 0)
      return( bessy0(x) );
   if (n == 1)
      return( bessy1(x) );

   tox=2.0/x;
   by=bessy1(x);
   bym=bessy0(x);
   for (j=1;j<n;j++) {
      byp=j*tox*by-bym;
      bym=by;
      by=byp;
   }
   return by;
}





static double bessi0( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=0.  */
/*------------------------------------------------------------*/
{
   double ax,ans;
   double y;


   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
         +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
   } else {
      y=3.75/ax;
      ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
         +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
         +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
         +y*0.392377e-2))))))));
   }
   return ans;
}




static double bessi1( double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=1.  */
/*------------------------------------------------------------*/
{
   double ax,ans;
   double y;


   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
         +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
   } else {
      y=3.75/ax;
      ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
         -y*0.420059e-2));
      ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
         +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
      ans *= (exp(ax)/sqrt(ax));
   }
   return x < 0.0 ? -ans : ans;
}




/*
#>            bessi.dc2

Function:     bessi

Purpose:      Evaluate Modified Bessel function of integer order.

Category:     MATH

File:         bessel.c

Author:       M.G.R. Vogelaar

Use:          #include "bessel.h"
              double   result; 
              result = bessi( int n,
                              double x )


              bessi    Return the Modified  Bessel function Iv(x) of 
                       integer order for input value x.
              n        Integer order of Bessel function.
              x        Double at which the function is evaluated.

                      
Description:  bessy evaluates at x the Modified Bessel function of 
              integer order n.
              This routine is NOT callable in FORTRAN.

Updates:      Jun 29, 1998: VOG, Document created.
#<
*/



double bessi( int n, double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) for n >= 0*/
/*------------------------------------------------------------*/
{
   int j;
   double bi,bim,bip,tox,ans;


   if (n < 0)
   {
      double   dblank;
      setdblank_c( &dblank );
      return( dblank );
   }
   if (n == 0)
      return( bessi0(x) );
   if (n == 1)
      return( bessi1(x) );


   if (x == 0.0)
      return 0.0;
   else {
      tox=2.0/fabs(x);
      bip=ans=0.0;
      bi=1.0;
      for (j=2*(n+(int) sqrt(ACC*n));j>0;j--) {
         bim=bip+j*tox*bi;
         bip=bi;
         bi=bim;
         if (fabs(bi) > BIGNO) {
            ans *= BIGNI;
            bi *= BIGNI;
            bip *= BIGNI;
         }
         if (j == n) ans=bip;
      }
      ans *= bessi0(x)/bi;
      return  x < 0.0 && n%2 == 1 ? -ans : ans;
   }
}




static double bessk0( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n=0.  */
/*------------------------------------------------------------*/
{
   double y,ans;

   if (x <= 2.0) {
      y=x*x/4.0;
      ans=(-log(x/2.0)*bessi0(x))+(-0.57721566+y*(0.42278420
         +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
         +y*(0.10750e-3+y*0.74e-5))))));
   } else {
      y=2.0/x;
      ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
         +y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
         +y*(-0.251540e-2+y*0.53208e-3))))));
   }
   return ans;
}




static double bessk1( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n=1.  */
/*------------------------------------------------------------*/
{
   double y,ans;

   if (x <= 2.0) {
      y=x*x/4.0;
      ans=(log(x/2.0)*bessi1(x))+(1.0/x)*(1.0+y*(0.15443144
         +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
         +y*(-0.110404e-2+y*(-0.4686e-4)))))));
   } else {
      y=2.0/x;
      ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
         +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
         +y*(0.325614e-2+y*(-0.68245e-3)))))));
   }
   return ans;
}




/*
#>            bessk.dc2

Function:     bessk

Purpose:      Evaluate Modified Bessel function Kv(x) of integer order.

Category:     MATH

File:         bessel.c

Author:       M.G.R. Vogelaar

Use:          #include "bessel.h"
              double   result; 
              result = bessk( int n,
                              double x )


              bessk    Return the Modified Bessel function Kv(x) of 
                       integer order for input value x.
              n        Integer order of Bessel function.
              x        Double at which the function is evaluated.

                      
Description:  bessk evaluates at x the Modified Bessel function Kv(x) of 
              integer order n.
              This routine is NOT callable in FORTRAN.

Updates:      Jun 29, 1998: VOG, Document created.
#<
*/



double bessk( int n, double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n >= 0*/
/* Note that for x == 0 the functions bessy and bessk are not */
/* defined and a blank is returned.                           */
/*------------------------------------------------------------*/
{
   int j;
   double bk,bkm,bkp,tox;


   if (n < 0 || x == 0.0)
   {
      double   dblank;
      setdblank_c( &dblank );
      return( dblank );
   }
   if (n == 0)
      return( bessk0(x) );
   if (n == 1)
      return( bessk1(x) );

   tox=2.0/x;
   bkm=bessk0(x);
   bk=bessk1(x);
   for (j=1;j<n;j++) {
      bkp=bkm+j*tox*bk;
      bkm=bk;
      bk=bkp;
   }
   return bk;
}

#undef setdblank_c
#undef ACC
#undef BIGNO
#undef BIGNI
