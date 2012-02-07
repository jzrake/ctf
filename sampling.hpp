

#ifndef __SamplingModule_HEADER__
#define __SamplingModule_HEADER__

#include <vector>

void Mara_prim_at_point(const double *r0, double *P1);
void Mara_prim_at_point_many(const double *Rin, double *Rlist, double *Plist,
                             int Nsamp);
std::vector<double> Mara_streamline_velocity(const double *r0, double s1, double ds,
					     double(*f)(double *P));
std::vector<double> Mara_streamline_magnetic(const double *r0, double s1, double ds,
					     double(*f)(double *P));

double Mara_streamline_scalars_velocity(double *P);
double Mara_streamline_scalars_magnetic(double *P);

#endif // __SamplingModule_HEADER__
