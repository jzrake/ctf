
#ifdef __cplusplus
extern "C" {
#endif

#ifndef __MaraWenoLibrary_HEADER__
#define __MaraWenoLibrary_HEADER__

  enum ReconstructOperation { PLM_C2L, PLM_C2R,
			      WENO5_FD_C2R, WENO5_FD_C2L,
			      WENO5_FV_C2R, WENO5_FV_C2L,
			      WENO5_FV_C2A, WENO5_FV_A2C };

  enum SmoothnessIndicator { OriginalJiangShu96,
			     ImprovedBorges08,
			     ImprovedShenZha10 };

  double reconstruct(const double *v, enum ReconstructOperation type);
  void reconstruct_set_smoothness_indicator(enum SmoothnessIndicator IS);
  void reconstruct_set_plm_theta(double theta);
  void reconstruct_set_shenzha10_A(double A);

#endif // __MaraWenoLibrary_HEADER__

#ifdef __cplusplus
}
#endif
