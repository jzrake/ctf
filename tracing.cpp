

#include <vector>
#include "hydro.hpp"
#include "sampling.hpp"


enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
enum { rho, pre, vx, vy, vz };             // Primitive



std::vector<double> Mara_streamline_velocity(const double *r0, double s1, double ds)
{
  const PhysicalDomain &domain = *HydroModule::Mara->domain;
  const int Nq = domain.get_Nq();

  double *P1 = new double[Nq];
  double *r1 = new double[3];
  double *v1 = new double[3];

  double s = 0.0;
  std::vector<double> points;
  std::memcpy(r1, r0, 3*sizeof(double));


  while (s < s1) {
    Mara_prim_at_point(r1, P1);

    v1[0] = P1[vx];
    v1[1] = P1[vy];
    v1[2] = P1[vz];

    const double v = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);

    double vhat[3] = { v1[0]/v, v1[1]/v, v1[2]/v };

    r1[0] += ds * vhat[0];
    r1[1] += ds * vhat[1];
    r1[2] += ds * vhat[2];

    points.push_back(r1[0]);
    points.push_back(r1[1]);
    points.push_back(r1[2]);

    s += ds;
  }

  delete [] P1;
  delete [] r1;
  delete [] v1;

  return points;
}


