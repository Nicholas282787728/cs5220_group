#include <nmmintrin.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifndef BLOCK_SIZE
#define BLOCK_SIZE ((int) 8)
#endif


const char* dgemm_desc = "Nick's dgemm";

// Tranpose of matrix
void trans(const double *A, double *At, const int m) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      At[i*m + j] = A[i + j*m];
    }
  }
}

// Attempting to use memcpy to copy something
// from an unaligned array into an aligned array
// gives weird compiler errors
void mcopy(const double *A, double *At, const int m) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      At[i*m + j] = A[i*m + j];
    }
  }
}

/*
 * Block matrix multiply kernel.
 * Inputs:
 *    A: 8-by-8 matrix in row major format.
 *    B: 8-by-8 matrix in column major format.
 * Outputs:
 *    C: 8-by-8 matrix in row major format.
 */
void kdgemm(double * restrict C, const double * restrict A,  const double * restrict B) {
  
  for (int i = 0; i < 8; i++) {
    
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
    
    __m128d ca0 = _mm_load_pd(A + i*8+0);
    __m128d ca1 = _mm_load_pd(A + i*8+2);
    __m128d ca2 = _mm_load_pd(A + i*8+4);
    __m128d ca3 = _mm_load_pd(A + i*8+6);
    
    cval = _mm_load_pd(C + i*8+0);
    
    c0 = _mm_mul_pd(ca0, _mm_load_pd(B + 0*8+0));
    c1 = _mm_mul_pd(ca1, _mm_load_pd(B + 0*8+2));
    c2 = _mm_mul_pd(ca2, _mm_load_pd(B + 0*8+4));
    c3 = _mm_mul_pd(ca3, _mm_load_pd(B + 0*8+6));
    c4 = _mm_add_pd(c0, c1);
    c5 = _mm_add_pd(c2, c3);
    c6 = _mm_add_pd(c4, c5);
    
    c0 = _mm_mul_pd(ca0, _mm_load_pd(B + 1*8+0));
    c1 = _mm_mul_pd(ca1, _mm_load_pd(B + 1*8+2));
    c2 = _mm_mul_pd(ca2, _mm_load_pd(B + 1*8+4));
    c3 = _mm_mul_pd(ca3, _mm_load_pd(B + 1*8+6));
    c4 = _mm_add_pd(c0, c1);
    c5 = _mm_add_pd(c2, c3);
    c8 = _mm_add_pd(c4, c5);
    c7 = _mm_hadd_pd(c6, c8);
    cval = _mm_add_pd(cval, c7);
    
    _mm_store_pd(C + i*8+0, cval);
    
    cval = _mm_load_pd(C + i*8+2);
    
    c0 = _mm_mul_pd(ca0, _mm_load_pd(B + 2*8+0));
    c1 = _mm_mul_pd(ca1, _mm_load_pd(B + 2*8+2));
    c2 = _mm_mul_pd(ca2, _mm_load_pd(B + 2*8+4));
    c3 = _mm_mul_pd(ca3, _mm_load_pd(B + 2*8+6));
    c4 = _mm_add_pd(c0, c1);
    c5 = _mm_add_pd(c2, c3);
    c6 = _mm_add_pd(c4, c5);
    
    c0 = _mm_mul_pd(ca0, _mm_load_pd(B + 3*8+0));
    c1 = _mm_mul_pd(ca1, _mm_load_pd(B + 3*8+2));
    c2 = _mm_mul_pd(ca2, _mm_load_pd(B + 3*8+4));
    c3 = _mm_mul_pd(ca3, _mm_load_pd(B + 3*8+6));
    c4 = _mm_add_pd(c0, c1);
    c5 = _mm_add_pd(c2, c3);
    c8 = _mm_add_pd(c4, c5);
    c7 = _mm_hadd_pd(c6, c8);
    cval = _mm_add_pd(cval, c7);
    
    _mm_store_pd(C + i*8+2, cval);
    
    cval = _mm_load_pd(C + i*8+4);
    
    c0 = _mm_mul_pd(ca0, _mm_load_pd(B + 4*8+0));
    c1 = _mm_mul_pd(ca1, _mm_load_pd(B + 4*8+2));
    c2 = _mm_mul_pd(ca2, _mm_load_pd(B + 4*8+4));
    c3 = _mm_mul_pd(ca3, _mm_load_pd(B + 4*8+6));
    c4 = _mm_add_pd(c0, c1);
    c5 = _mm_add_pd(c2, c3);
    c6 = _mm_add_pd(c4, c5);
    
    c0 = _mm_mul_pd(ca0, _mm_load_pd(B + 5*8+0));
    c1 = _mm_mul_pd(ca1, _mm_load_pd(B + 5*8+2));
    c2 = _mm_mul_pd(ca2, _mm_load_pd(B + 5*8+4));
    c3 = _mm_mul_pd(ca3, _mm_load_pd(B + 5*8+6));
    c4 = _mm_add_pd(c0, c1);
    c5 = _mm_add_pd(c2, c3);
    c8 = _mm_add_pd(c4, c5);
    c7 = _mm_hadd_pd(c6, c8);
    cval = _mm_add_pd(cval, c7);
    
    _mm_store_pd(C + i*8+4, cval);
    
    cval = _mm_load_pd(C + i*8+6);
    
    c0 = _mm_mul_pd(ca0, _mm_load_pd(B + 6*8+0));
    c1 = _mm_mul_pd(ca1, _mm_load_pd(B + 6*8+2));
    c2 = _mm_mul_pd(ca2, _mm_load_pd(B + 6*8+4));
    c3 = _mm_mul_pd(ca3, _mm_load_pd(B + 6*8+6));
    c4 = _mm_add_pd(c0, c1);
    c5 = _mm_add_pd(c2, c3);
    c6 = _mm_add_pd(c4, c5);
    
    c0 = _mm_mul_pd(ca0, _mm_load_pd(B + 7*8+0));
    c1 = _mm_mul_pd(ca1, _mm_load_pd(B + 7*8+2));
    c2 = _mm_mul_pd(ca2, _mm_load_pd(B + 7*8+4));
    c3 = _mm_mul_pd(ca3, _mm_load_pd(B + 7*8+6));
    c4 = _mm_add_pd(c0, c1);
    c5 = _mm_add_pd(c2, c3);
    c8 = _mm_add_pd(c4, c5);
    c7 = _mm_hadd_pd(c6, c8);
    cval = _mm_add_pd(cval, c7);
    
    
    _mm_store_pd(C + i*8+6, cval);
    
  }
  
}

void square_dgemm(const int M, const double *A, const double *B, double *C) {
  
  const int M2 = 8 * 8;
  
  static double transA[64] __attribute__((aligned(16))) = {0};
  static double alignedB[64] __attribute__((aligned(16))) = {0};
  static double transC[64] __attribute__((aligned(16))) = {0};
  
  trans(A, transA, 8);
  mcopy(B, alignedB, 8);
  trans(C, transC, 8);
  
  kdgemm(transC, transA, B);
  

  trans(transC, C, 8);
  
}
