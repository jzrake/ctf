

#include <vector>
#include "hydro.hpp"
#include "sampling.hpp"


enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
enum { rho, pre, vx, vy, vz };             // Primitive


static void _streamline(const double *r0, double s1, double ds,
			int fx, int fy, int fz,
			std::vector<double> &points);


std::vector<double> Mara_streamline_velocity(const double *r0, double s1, double ds)
{
  std::vector<double> points;
  _streamline(r0, s1, ds, vx, vy, vz, points);
  return points;
}

std::vector<double> Mara_streamline_magnetic(const double *r0, double s1, double ds)
{
  std::vector<double> points;
  _streamline(r0, s1, ds, Bx, By, Bz, points);
  return points;
}




void _streamline(const double *r0, double s1, double ds,
		 int fx, int fy, int fz,
		 std::vector<double> &points)
{
  const PhysicalDomain &domain = *HydroModule::Mara->domain;
  const int Nq = domain.get_Nq();

  double *P1 = new double[Nq];

  double r[3];
  double r1[3], r2[3], r3[3], r4[3];
  double k1[3], k2[3], k3[3], k4[3];

  double s = 0.0;
  std::memcpy(r, r0, 3*sizeof(double));

  while (s < s1) {

    Mara_prim_at_point(r, P1);
    {
      double v1[3] = { P1[fx], P1[fy], P1[fz] };
      double v = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
      if (fabs(v) < 1e-14) v = 1e-14;

      k1[0] = ds * v1[0]/v;
      k1[1] = ds * v1[1]/v;
      k1[2] = ds * v1[2]/v;

      r1[0] = r[0] + 0.5*k1[0];
      r1[1] = r[1] + 0.5*k1[1];
      r1[2] = r[2] + 0.5*k1[2];
    }
    Mara_prim_at_point(r1, P1);
    {
      double v1[3] = { P1[fx], P1[fy], P1[fz] };
      double v = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
      if (fabs(v) < 1e-14) v = 1e-14;

      k2[0] = 0.5*ds * v1[0]/v;
      k2[1] = 0.5*ds * v1[1]/v;
      k2[2] = 0.5*ds * v1[2]/v;

      r2[0] = r[0] + 0.5*k2[0];
      r2[1] = r[1] + 0.5*k2[1];
      r2[2] = r[2] + 0.5*k2[2];
    }
    Mara_prim_at_point(r2, P1);
    {
      double v1[3] = { P1[fx], P1[fy], P1[fz] };
      double v = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
      if (fabs(v) < 1e-14) v = 1e-14;

      k3[0] = 0.5*ds * v1[0]/v;
      k3[1] = 0.5*ds * v1[1]/v;
      k3[2] = 0.5*ds * v1[2]/v;

      r3[0] = r[0] + k3[0];
      r3[1] = r[1] + k3[1];
      r3[2] = r[2] + k3[2];
    }
    Mara_prim_at_point(r3, P1);
    {
      double v1[3] = { P1[fx], P1[fy], P1[fz] };
      double v = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
      if (fabs(v) < 1e-14) v = 1e-14;

      k4[0] = ds * v1[0]/v;
      k4[1] = ds * v1[1]/v;
      k4[2] = ds * v1[2]/v;

      r4[0] = r[0] + k4[0];
      r4[1] = r[1] + k4[1];
      r4[2] = r[2] + k4[2];
    }

    r[0] += (1./6.)*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
    r[1] += (1./6.)*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    r[2] += (1./6.)*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);

    points.push_back(r[0]);
    points.push_back(r[1]);
    points.push_back(r[2]);

    s += ds;
  }

  delete [] P1;
}


