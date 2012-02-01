
/*------------------------------------------------------------------------------
 * FILE: matrix.h
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP
 *
 * DESCRIPTION: Header file for a simple matrix module
 *
 * REFERENCES:
 *------------------------------------------------------------------------------
 */


#ifdef __cplusplus
extern "C" {
#endif

#ifndef __MatrixModule_HEADER__
#define __MatrixModule_HEADER__

double matrix_determinant();
int matrix_inverse(const double *A, double *B, int rank);
int matrix_matrix_product(const double *A, const double *B, double *C,
			  int ni, int nj, int nk);
int matrix_vector_product(const double *A, const double *x, double *b,
			  int ni, int nj);
int vector_matrix_product(const double *A, const double *x, double *b,
			  int ni, int nj);

#endif // __MatrixModule_HEADER__

#ifdef __cplusplus
}
#endif
