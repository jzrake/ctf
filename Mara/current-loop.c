#include <math.h>
#include "current-loop.h"
#define PI (4.0*atan(1.0))


/* private functions */
static double integrate_RK4(double k, double (*f)(double k, double t));
static double elliptic_K_integrand(double k2, double t);
static double elliptic_E_integrand(double k2, double t);
static double K_ell(double k2);
static double E_ell(double k2);


double integrate_RK4(double k, double (*f)(double k, double t))
{
  int nstep = 10;
  int niter = 0;
  double dt = (0.5 * PI) / nstep;

  double t = 0.0;
  double y = 0.0;

  while (niter++ < nstep) {

    double k1 = f(k, t + 0.0 * dt);
    double k2 = f(k, t + 0.5 * dt);
    double k3 = f(k, t + 0.5 * dt);
    double k4 = f(k, t + 1.0 * dt);

    y += dt * (k1 + 2*k2 + 2*k3 + k4) / 6.0;
    t += dt;
  }

  return y;
}

double elliptic_K_integrand(double k2, double t)
{
  if (fabs(k2 - 1.0) < 1e-12) { /* little fudge to avoid inf's */
    k2 -= 1e-12;
  }
  return 1.0 / sqrt(1.0 - k2*sin(t)*sin(t));
}

double elliptic_E_integrand(double k2, double t)
{
  if (fabs(k2 - 1.0) < 1e-12) { /* little fudge to avoid inf's */
    k2 -= 1e-12;
  }
  return 1.0 * sqrt(1.0 - k2*sin(t)*sin(t));
}

double K_ell(double k2)
{
  //  printf("K got k2=%f\n", k2);
  return integrate_RK4(k2, elliptic_K_integrand);
}

double E_ell(double k2)
{
  //  printf("E got k2=%f\n", k2);
  return integrate_RK4(k2, elliptic_E_integrand);
}

double vector_potential_phi_exact(double I, double a, double r, double t)
/*
 * Return the phi-component of vector potential due to a current loop at the
 * origin.
 *
 *
 * I: current, in SI is actually mu_0 I / (4 pi)
 * a: radius of current loop
 * r: distance from origin (loop center)
 * t: theta (polar angle, from loop normal)
 */
{
  double num = 4*a*r*sin(t);
  double den = r*r + a*a + 2*a*r*sin(t);
  double k2 = num / den;

  if (fabs(num) < 1e-12) {
    return 0.0;
  }
  else if (fabs(den) < 1e-12) {
    return 0.0;
  }
  else {
    return -((1.0 / sqrt(r*r + a*a + 2*a*r*sin(t))) *
	     (4*I*a/k2) * ((k2 - 2) * K_ell(k2) + 2*E_ell(k2)));
  }
}

double vector_potential_phi_far(double I, double a, double r, double t)
{
  if (t < 1e-12) {
    return 0.0;
  }
  else {
    return I*PI*a*a * sin(t) / (r*r);
  }
}

void vector_potential_cartesian(double I, double a, double x[3], double A[3])
{
  double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

  if (fabs(r) < 1e-12) {
    A[0] = 0.0;
    A[1] = 0.0;
    A[2] = 0.0;
    return;
  }

  double t = acos(x[2] / r);
  double Aphi = vector_potential_phi_exact(I, a, r, t);
  A[0] = -Aphi * x[1]/r;
  A[1] =  Aphi * x[0]/r;
  A[2] =  0.0;
}

void current_loop_magnetic_field(double I, double a, double x[3],
				 double dx[3], double B[3])
/*
 * Return the magnetic field vector B of the current loop at x in cartesian
 * coordinates. dx is the numerical step used to finite differenc the vector
 * potential, which will respect Mara's "corner" type divergence stencil.
 *
 */
{
  double XxL0[3] = { x[0] - 0.5*dx[0], x[1] - 0.5*dx[1], x[2] - 0.5*dx[2] };
  double XxR0[3] = { x[0] + 0.5*dx[0], x[1] - 0.5*dx[1], x[2] - 0.5*dx[2] };
  double XxL1[3] = { x[0] - 0.5*dx[0], x[1] - 0.5*dx[1], x[2] + 0.5*dx[2] };
  double XxR1[3] = { x[0] + 0.5*dx[0], x[1] - 0.5*dx[1], x[2] + 0.5*dx[2] };
  double XxL2[3] = { x[0] - 0.5*dx[0], x[1] + 0.5*dx[1], x[2] - 0.5*dx[2] };
  double XxR2[3] = { x[0] + 0.5*dx[0], x[1] + 0.5*dx[1], x[2] - 0.5*dx[2] };
  double XxL3[3] = { x[0] - 0.5*dx[0], x[1] + 0.5*dx[1], x[2] + 0.5*dx[2] };
  double XxR3[3] = { x[0] + 0.5*dx[0], x[1] + 0.5*dx[1], x[2] + 0.5*dx[2] };

  double XyL0[3] = { x[0] - 0.5*dx[0], x[1] - 0.5*dx[1], x[2] - 0.5*dx[2] };
  double XyR0[3] = { x[0] - 0.5*dx[0], x[1] + 0.5*dx[1], x[2] - 0.5*dx[2] };
  double XyL1[3] = { x[0] + 0.5*dx[0], x[1] - 0.5*dx[1], x[2] - 0.5*dx[2] };
  double XyR1[3] = { x[0] + 0.5*dx[0], x[1] + 0.5*dx[1], x[2] - 0.5*dx[2] };
  double XyL2[3] = { x[0] - 0.5*dx[0], x[1] - 0.5*dx[1], x[2] + 0.5*dx[2] };
  double XyR2[3] = { x[0] - 0.5*dx[0], x[1] + 0.5*dx[1], x[2] + 0.5*dx[2] };
  double XyL3[3] = { x[0] + 0.5*dx[0], x[1] - 0.5*dx[1], x[2] + 0.5*dx[2] };
  double XyR3[3] = { x[0] + 0.5*dx[0], x[1] + 0.5*dx[1], x[2] + 0.5*dx[2] };

  double XzL0[3] = { x[0] - 0.5*dx[0], x[1] - 0.5*dx[1], x[2] - 0.5*dx[2] };
  double XzR0[3] = { x[0] - 0.5*dx[0], x[1] - 0.5*dx[1], x[2] + 0.5*dx[2] };
  double XzL1[3] = { x[0] - 0.5*dx[0], x[1] + 0.5*dx[1], x[2] - 0.5*dx[2] };
  double XzR1[3] = { x[0] - 0.5*dx[0], x[1] + 0.5*dx[1], x[2] + 0.5*dx[2] };
  double XzL2[3] = { x[0] + 0.5*dx[0], x[1] - 0.5*dx[1], x[2] - 0.5*dx[2] };
  double XzR2[3] = { x[0] + 0.5*dx[0], x[1] - 0.5*dx[1], x[2] + 0.5*dx[2] };
  double XzL3[3] = { x[0] + 0.5*dx[0], x[1] + 0.5*dx[1], x[2] - 0.5*dx[2] };
  double XzR3[3] = { x[0] + 0.5*dx[0], x[1] + 0.5*dx[1], x[2] + 0.5*dx[2] };

  double AxL0[3], AyL0[3], AzL0[3];
  double AxR0[3], AyR0[3], AzR0[3];
  double AxL1[3], AyL1[3], AzL1[3];
  double AxR1[3], AyR1[3], AzR1[3];
  double AxL2[3], AyL2[3], AzL2[3];
  double AxR2[3], AyR2[3], AzR2[3];
  double AxL3[3], AyL3[3], AzL3[3];
  double AxR3[3], AyR3[3], AzR3[3];

  vector_potential_cartesian(I, a, XxL0, AxL0);
  vector_potential_cartesian(I, a, XxR0, AxR0);
  vector_potential_cartesian(I, a, XxL1, AxL1);
  vector_potential_cartesian(I, a, XxR1, AxR1);
  vector_potential_cartesian(I, a, XxL2, AxL2);
  vector_potential_cartesian(I, a, XxR2, AxR2);
  vector_potential_cartesian(I, a, XxL3, AxL3);
  vector_potential_cartesian(I, a, XxR3, AxR3);

  vector_potential_cartesian(I, a, XyL0, AyL0);
  vector_potential_cartesian(I, a, XyR0, AyR0);
  vector_potential_cartesian(I, a, XyL1, AyL1);
  vector_potential_cartesian(I, a, XyR1, AyR1);
  vector_potential_cartesian(I, a, XyL2, AyL2);
  vector_potential_cartesian(I, a, XyR2, AyR2);
  vector_potential_cartesian(I, a, XyL3, AyL3);
  vector_potential_cartesian(I, a, XyR3, AyR3);

  vector_potential_cartesian(I, a, XzL0, AzL0);
  vector_potential_cartesian(I, a, XzR0, AzR0);
  vector_potential_cartesian(I, a, XzL1, AzL1);
  vector_potential_cartesian(I, a, XzR1, AzR1);
  vector_potential_cartesian(I, a, XzL2, AzL2);
  vector_potential_cartesian(I, a, XzR2, AzR2);
  vector_potential_cartesian(I, a, XzL3, AzL3);
  vector_potential_cartesian(I, a, XzR3, AzR3);

  double dA[3][3]; /* dA[i][j] = dA_i / dx_j */
  int i;

  for (i=0; i<3; ++i) {
    dA[i][0] = 0.25 * ((AxR0[i] + AxR1[i] + AxR2[i] + AxR3[i]) -
		       (AxL0[i] + AxL1[i] + AxL2[i] + AxL3[i])) / dx[0];
    dA[i][1] = 0.25 * ((AyR0[i] + AyR1[i] + AyR2[i] + AyR3[i]) -
		       (AyL0[i] + AyL1[i] + AyL2[i] + AyL3[i])) / dx[1];
    dA[i][2] = 0.25 * ((AzR0[i] + AzR1[i] + AzR2[i] + AzR3[i]) -
		       (AzL0[i] + AzL1[i] + AzL2[i] + AzL3[i])) / dx[2];
  }

  B[0] = dA[2][1] - dA[1][2];
  B[1] = dA[0][2] - dA[2][0];
  B[2] = dA[1][0] - dA[0][1];
}

#ifdef __MAIN__
#include <stdio.h>

int main()
{

  double k2 = 0.1;
  double K = K_ell(k2);
  double E = E_ell(k2);

  /* wolfram: K(0.1) */
  double Ktrue = 1.6124413487202193982299163630853741545268463794958952;

  printf("K(%f) = %12.10f\n", k2, K);
  printf("E(%f) = %12.10f\n", k2, E);
  printf("error in K = %e\n", K - Ktrue);

  double A_exact = vector_potential_phi_exact(1.5, 2.0, 100.0, 0.5);
  double A_far   = vector_potential_phi_far  (1.5, 2.0, 100.0, 0.5);

  printf("far   : %10.8e\n", A_far);
  printf("exact : %10.8e\n", A_exact);

  {
    double B[3];
    double x[3] = {-0.09375, -0.03125, 0.03125 };
    double dx[3] = { 1.0/16, 1.0/16, 1.0/16 };
    current_loop_magnetic_field(1.0, 0.125, x, dx, B);
    printf("B = %f %f %f\n\n", B[0], B[1], B[2]);
  }
  {
    double B[3];
    double x[3] = {-0.03125, -0.03125, 0.96875};
    double dx[3] = { 1.0/16, 1.0/16, 1.0/16 };
    current_loop_magnetic_field(1.0, 2.0, x, dx, B);
    printf("B = %f %f %f\n\n", B[0], B[1], B[2]);
  }
  {
    double B[3];
    double x[3] = {-0.03125, -0.03125, -0.03125};
    double dx[3] = { 1.0/16, 1.0/16, 1.0/16 };
    current_loop_magnetic_field(1.0, 2.0, x, dx, B);
    printf("B = %f %f %f\n\n", B[0], B[1], B[2]);
  }

  return 0;
}

#endif
