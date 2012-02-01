
/*------------------------------------------------------------------------------
 * FILE: matmult.c
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP
 *
 * DESCRIPTION: Calculates the matrix-matrix and matrix-vector product
 *
 * REFERENCES:
 *
 *------------------------------------------------------------------------------
 */

#include "matrix.h"

int matrix_matrix_product(const double *A, const double *B, double *C,
			  int ni, int nj, int nk)
{
  /* ---------------------------------------------------------------------------
   * Multiplies matrices C = A*B according to
   *
   * C_ij = Sum_k{ A_ik * B_kj }
   *
   * A is (ni x nk)
   * B is (nk x nj)
   * C is (ni x nj)
   * ---------------------------------------------------------------------------
   */
  int i,j,k;
  for (i=0; i<ni; ++i) {
    for (j=0; j<nj; ++j) {
      C[i*nj+j] = 0.0;
      for (k=0; k<nk; ++k) {
	C[i*nj+j] += A[i*nk+k] * B[k*nj+j];
      }
    }
  }
  return 0;
}

int matrix_vector_product(const double *A, const double *x, double *b,
			  int ni, int nj)
{
  /* ---------------------------------------------------------------------------
   * Multiplies matrix to vector b = A*x according to
   *
   * b_i = Sum_j{ A_ij * x_j }
   *
   * A is (ni x nj)
   * x is (nj x 1)
   * b is (ni x 1)
   * ---------------------------------------------------------------------------
   */
  int i,j;
  for (i=0; i<ni; ++i) {
    b[i] = 0.0;
    for (j=0; j<nj; ++j) {
      b[i] += A[i*nj+j] * x[j];
    }
  }
  return 0;
}

int vector_matrix_product(const double *A, const double *x, double *b,
			  int ni, int nj)
{
  /* ---------------------------------------------------------------------------
   * Multiplies matrix to vector b = x*A according to
   *
   * b_j = Sum_i{ x_i * A_ij }
   *
   * A is (ni x nj)
   * x is (1  x nj)
   * b is (1  x ni)
   * ---------------------------------------------------------------------------
   */
  int i,j;
  for (j=0; j<nj; ++j) {
    b[j] = 0.0;
    for (i=0; i<ni; ++i) {
      b[j] += x[i] * A[i*nj+j];
    }
  }
  return 0;
}
