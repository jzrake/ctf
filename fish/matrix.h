
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

#ifndef MATRIX_HEADER_INCLUDED
#define MATRIX_HEADER_INCLUDED

int matrix_matrix_product(const double *A, const double *B, double *C, int ni,
			  int nj, int nk);
int matrix_vector_product(const double *A, const double *x, double *b,
			  int ni, int nj);
int vector_matrix_product(const double *A, const double *x, double *b,
			  int ni, int nj);
double matrix_inverse(const double *A, double *B, int rank);

#endif // MATRIX_HEADER_INCLUDED
