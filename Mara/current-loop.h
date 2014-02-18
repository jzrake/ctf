
#ifdef __cplusplus
extern "C" {
#endif

#ifndef __MaraCurrentLoop_HEADER__
#define __MaraCurrentLoop_HEADER__

void current_loop_magnetic_field(double I, double a, double x[3],
				 double dx[3], double B[3]);

#endif // __MaraCurrentLoop_HEADER__

#ifdef __cplusplus
}
#endif
