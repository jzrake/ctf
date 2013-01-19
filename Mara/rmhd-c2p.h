

#ifndef __RmhdConsToPrim_HEADER__
#define __RmhdConsToPrim_HEADER__

#ifdef __cplusplus
extern "C" {
#endif

  enum RmhdConsToPrimError {
    RMHD_C2P_SUCCESS,
    RMHD_C2P_CONS_CONTAINS_NAN,
    RMHD_C2P_CONS_NEGATIVE_DENSITY,
    RMHD_C2P_CONS_NEGATIVE_ENERGY,
    RMHD_C2P_PRIM_CONTAINS_NAN,
    RMHD_C2P_PRIM_NEGATIVE_PRESSURE,
    RMHD_C2P_PRIM_NEGATIVE_RESTMASS,
    RMHD_C2P_PRIM_SUPERLUMINAL,
    RMHD_C2P_MAXITER
  } ;

  void rmhd_c2p_set_gamma(double adiabatic_gamma);
  void rmhd_c2p_new_state(const double *U);
  void rmhd_c2p_estimate_from_cons();
  void rmhd_c2p_set_starting_prim(const double *P);
  void rmhd_c2p_get_starting_prim(double *P);
  int rmhd_c2p_reconstruct_prim(double Z, double W, double *P);
  int rmhd_c2p_solve_anton2dzw(double *P);
  int rmhd_c2p_solve_noble1dw(double *P);
  int rmhd_c2p_get_iterations();
  int rmhd_c2p_check_cons(const double *U);
  int rmhd_c2p_check_prim(const double *P);
  char *rmhd_c2p_get_error(int error);

  void rmhd_c2p_eos_set_eos(const void *eos_);
  void rmhd_c2p_eos_new_state(const double *U);
  void rmhd_c2p_eos_estimate_from_cons();
  void rmhd_c2p_eos_set_starting_prim(const double *P);
  int rmhd_c2p_eos_solve_noble2dzt(double *P);
  int rmhd_c2p_eos_solve_duffell3d(double *P);
  int rmhd_c2p_eos_get_iterations();

#ifdef __cplusplus
}
#endif



#endif // __RmhdConsToPrim_HEADER__
