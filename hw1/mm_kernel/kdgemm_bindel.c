#include <nmmintrin.h>

/*
 * Dimensions for a "kernel" multiply.  We use define statements in
 * order to make sure these are treated as compile-time constants
 * (which the optimizer likes)
 */
#define M 2
#define N 2
#define P 5000

#ifdef USE_SHUFPD
#  define swap_sse_doubles(a) _mm_shuffle_pd(a, a, 1)
#else
#  define swap_sse_doubles(a) (__m128d) _mm_shuffle_epi32((__m128i) a, 0x4e)
#endif

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
 *    A: 2-by-P matrix in column major format.
 *    B: P-by-2 matrix in row major format.
 * Outputs:
 *    C: 2-by-2 matrix with element order [c11, c22, c12, c21]
 *       (diagonals stored first, then off-diagonals)
 */
void kdgemm(const double * restrict A,
            const double * restrict B,
            double * restrict C) {
  // Load diagonal and off-diagonals
  __m128d cd = _mm_load_pd(C+0);
  __m128d co = _mm_load_pd(C+2);
  
  /*
   * Do block dot product.  Each iteration adds the result of a two-by-two
   * matrix multiply into the accumulated 2-by-2 product matrix, which is
   * stored in the registers cd (diagonal part) and co (off-diagonal part).
   */
  for (int k = 0; k < P; k += 2) {
    
    __m128d a0 = _mm_load_pd(A+2*k+0);
    __m128d b0 = _mm_load_pd(B+2*k+0);
    __m128d td0 = _mm_mul_pd(a0, b0);
    __m128d bs0 = swap_sse_doubles(b0);
    __m128d to0 = _mm_mul_pd(a0, bs0);
    
    __m128d a1 = _mm_load_pd(A+2*k+2);
    __m128d b1 = _mm_load_pd(B+2*k+2);
    __m128d td1 = _mm_mul_pd(a1, b1);
    __m128d bs1 = swap_sse_doubles(b1);
    __m128d to1 = _mm_mul_pd(a1, bs1);
    
    __m128d td_sum = _mm_add_pd(td0, td1);
    __m128d to_sum = _mm_add_pd(to0, to1);
    
    cd = _mm_add_pd(cd, td_sum);
    co = _mm_add_pd(co, to_sum);
  }
  
  // Write back sum
  _mm_store_pd(C+0, cd);
  _mm_store_pd(C+2, co);
}

/*
 * Conversion routines that take a matrix block in column-major form
 * and put it into whatever form the kdgemm routine likes.
 */

void to_kdgemm_A(int ldA, const double* restrict A, double * restrict Ak) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      Ak[i*P + j] = A[i*P + j];
    }
  }
}

void to_kdgemm_B(int ldB, const double* restrict B, double * restrict Bk) {
  for (int i = 0; i < P; i++) {
    for (int j = 0; j < N; j++) {
      Bk[i*M + j] = B[i + j*P];
    }
  }
}

void from_kdgemm_C(int ldC, const double* restrict Ck, double * restrict C) {
  C[0] = Ck[0];
  C[1] = Ck[3];
  C[2] = Ck[2];
  C[3] = Ck[1];
}
