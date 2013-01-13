
#ifdef __cplusplus
extern "C" {
#endif

#ifndef __MaraMpiWrappers_HEADER__
#define __MaraMpiWrappers_HEADER__

int Mara_mpi_active();
int Mara_mpi_get_rank();
int Mara_mpi_get_size();
void Mara_mpi_barrier();
double Mara_mpi_dbl_min(double myval);
double Mara_mpi_dbl_max(double myval);
double Mara_mpi_dbl_sum(double myval);
int Mara_mpi_int_prod(int myval);
int Mara_mpi_int_sum(int myval);

#endif // __MaraMpiWrappers_HEADER__

#ifdef __cplusplus
}
#endif
