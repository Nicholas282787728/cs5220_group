#include <nmmintrin.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*
 * Dimensions for a "kernel" multiply.  We use define statements in
 * order to make sure these are treated as compile-time constants
 * (which the optimizer likes)
 */
#define M 4
#define N 4
#define P 4

/*
 * The ktimer driver expects these variables to be set to whatever
 * the dimensions of a kernel multiply are.  It uses them both for
 * space allocation and for flop rate computations.
 */
int DIM_M=M;
int DIM_N=N;
int DIM_P=P;

/*
 * Block matrix multiply kernel.
 * Inputs:
 *    A: 8-by-8 matrix in row major format.
 *    B: 8-by-8 matrix in column major format.
 * Outputs:
 *    C: 8-by-8 matrix in row major format.
 */
void kdgemm(const double * restrict A,
            const double * restrict B,
            double * restrict C) {
  
#pragma unroll(4)
  for (int i = 0; i < P; i++) {
    
    __m128d c0;
    __m128d c1;
    __m128d c2;
    __m128d c3;
    __m128d c4;
    __m128d c5;
    __m128d c6;
    __m128d c7;
    __m128d c8;
    __m128d cval;
    
    __m128d ca0 = _mm_load_pd(A + i*4+0);
    __m128d ca1 = _mm_load_pd(A + i*4+2);
    
    cval = _mm_load_pd(C + i*4+0);
    
    c0 = _mm_mul_pd(ca0, _mm_load_pd(B + 0*4+0));
    c1 = _mm_mul_pd(ca1, _mm_load_pd(B + 0*4+2));
    c4 = _mm_add_pd(c0, c1);
    
    c0 = _mm_mul_pd(ca0, _mm_load_pd(B + 1*4+0));
    c1 = _mm_mul_pd(ca1, _mm_load_pd(B + 1*4+2));
    c5 = _mm_add_pd(c0, c1);
    c7 = _mm_hadd_pd(c4, c5);
    cval = _mm_add_pd(cval, c7);
    
    _mm_store_pd(C + i*4+0, cval);
    
    cval = _mm_load_pd(C + i*4+2);
    
    c0 = _mm_mul_pd(ca0, _mm_load_pd(B + 2*4+0));
    c1 = _mm_mul_pd(ca1, _mm_load_pd(B + 2*4+2));
    c4 = _mm_add_pd(c0, c1);
    
    c0 = _mm_mul_pd(ca0, _mm_load_pd(B + 3*4+0));
    c1 = _mm_mul_pd(ca1, _mm_load_pd(B + 3*4+2));
    c5 = _mm_add_pd(c0, c1);
    c7 = _mm_hadd_pd(c4, c5);

    cval = _mm_add_pd(cval, c7);
    
    _mm_store_pd(C + i*4+2, cval);

    
  }
  
}

/*
 * Conversion routines that take a matrix block in column-major form
 * and put it into whatever form the kdgemm routine likes.
 */

void to_kdgemm_A(int ldA, const double* restrict A, double * restrict Ak) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      Ak[i*ldA + j] = A[i + j*M];
    }
  }
}

void to_kdgemm_B(int ldB, const double* restrict B, double * restrict Bk) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      Bk[i*P + j] = B[i*ldB + j];
    }
  }
}

void from_kdgemm_C(int ldC, const double* restrict Ck, double * restrict C) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      C[i*ldC + j] = Ck[i + j*M];
    }
  }
}
