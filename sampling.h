

#ifndef __SamplingModule_HEADER__
#define __SamplingModule_HEADER__

extern "C" {
void Mara_prim_at_point(const double *r0, double *P1);
void Mara_prim_at_point_many(const double *Rin, double *Rlist, double *Plist,
                             int Nsamp);
}

#endif // __SamplingModule_HEADER__
