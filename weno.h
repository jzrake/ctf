
#ifdef __cplusplus
extern "C" {
#endif

#ifndef __MaraWenoLibrary_HEADER__
#define __MaraWenoLibrary_HEADER__

enum WenoOperation { PLM_C2L, PLM_C2R,
		     WENO5_FD_C2R, WENO5_FD_C2L,
		     WENO5_FV_C2R, WENO5_FV_C2L,
		     WENO5_FV_C2A, WENO5_FV_A2C };
double reconstruct(const double *v, enum WenoOperation type);
void reconstruct_set_plm_theta(double theta);

#endif // __MaraWenoLibrary_HEADER__

#ifdef __cplusplus
}
#endif
