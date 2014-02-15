#include <nmmintrin.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

const char* dgemm_desc = "Simple blocked dgemm.";


/*
 ** On the Nehalem architecture, shufpd and multiplication use the same port.
 ** 32-bit integer shuffle is a different matter.  If we want to try to make
 ** it as easy as possible for the compiler to schedule multiplies along
 ** with adds, it therefore makes sense to abuse the integer shuffle
 ** instruction.  See also
 **   http://locklessinc.com/articles/interval_arithmetic/
 **/
#ifdef USE_SHUFPD
#  define swap_sse_doubles(a) _mm_shuffle_pd(a, a, 1)
#else
#  define swap_sse_doubles(a) (__m128d) _mm_shuffle_epi32((__m128i) a, 0x4e)
#endif

#ifndef BLOCK_SIZE
#define BLOCK_SIZE ((int) 2)
#endif

/*
   A is M-by-K
   B is K-by-N
   C is M-by-N

   lda is the leading dimension of the matrix (the M of square_dgemm).
   */
void kdgemm2P2(const int P, double * restrict C,
    const double * restrict A,
    const double * restrict B)
{
  // This is really implicit in using the aligned ops...
  //__assume_aligned(A, 16);
  //__assume_aligned(B, 16);
  //__assume_aligned(C, 16);
  __builtin_assume_aligned(A, 16);
  __builtin_assume_aligned(B, 16);
  __builtin_assume_aligned(C, 16);

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


// Tranpose of matrix
void trans(const double *A, double *At, const int m) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      At[i*m + j] = A[i + j*m];
    }
  }
}

// Print the square matrix (row major)
void printMatrixC(const double *A, const int n) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      printf("%f ", A[j*n + i]);
    }
    printf("\n");
  }
}

// Print the square matrix (column major) 
void printMatrixR(const double *A, const int n) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      printf("%f ", A[i + j*n]);
    }
    printf("\n");
  }
}

// n is the size of A, m is the number of blocks
void bijectA(const double *A, double *At, const int n, const int m, const int odd) {
  if (odd == 0) {
    for (int j = 0; j<n; j++) {
      for (int i = 0; i < m; i++) {
        At[i*n*2 + 2*j]     = A[j*n + i*2];
        At[i*n*2 + 2*j + 1] = A[j*n + i*2 + 1];
      }
    }
  }
  else {
    // First loop over the blocks
    for (int j = 0; j< n-1; j++) {
      // Then for each block, go to each column, except the last one
      for (int i = 0; i < m; i++) {
        At[i*n*2 + 2*j] = A[j*(n-1) + i*2];
        //printf("Stored (%d, %d) in (%d, %d)\n", i*2, j*(n-1), i*n*2, 2*j);

        // For the case where of the last row being padded with zeros
        if (i != m-1) {
          At[i*n*2 + 2*j +1] = A[j*(n-1) + i*2 +1];
          //printf("Stored %d in %d\n", i*2+1+ j*(n-1), i*n*2+ 2*j+1);
        }
        else {
          At[i*n*2 + 2*j + 1] = 0;
          //printf("Stored %d in %d\n", 0, i*n*2+ 2*j+1);
        }
        // Last columns padded with zero
        At[i*n*2 + 2*(n-1)] = 0;
        At[i*n*2 + 2*(n-1) + 1] = 0;

      }

    }
  }
}

// n is the size of B, m is the number of blocks
void bijectB(const double *B, double *Bt, const int n, const int m, const int odd) {
  // Slightly different style from the bijectA function, but same thing.
  int n2 = n*n;
  if (odd == 0) {
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        Bt[n2 + i*n*2 + 2*j]     = B[j + (i*2)*n];
        Bt[n2 + i*n*2 + 2*j + 1] = B[j + (i*2 + 1)*n];
      }
    }
  }

  else {
    for (int i = 0; i< m; i++) {
      for (int j = 0; j < n-1; j++) {
        Bt[n2 + i*n*2 + 2*j] = B[j + i*2*(n-1)];

        if (i != m-1) {
          Bt[n2 + i*n*2 + 2*j + 1] = B[j + (i*2 + 1)*(n-1)];
        }
        else {
          Bt[n2 + i*n*2 + 2*j + 1] = 0;
        }
      }

      // Last columns padded with zero
      Bt[n2 + i*n*2 + 2*(n-1)] = 0;
      Bt[n2 + i*n*2 + 2*(n-1) + 1] = 0;
    }
  }
}

void square_dgemm(const int M, const double *A, const double *B, double *C)
{
  int bi, bj, bk, temp, n;
  int odd = 0; 

  // Check if the matrices is even or not.
  if (M % 2 == 1) {
    odd = 1;
    n = M+1;
  }
  else {
    n = M;
  }

  // Used for just setting up the arrays
  temp = n*n;

  // Setup
  const int n_blocks = M / BLOCK_SIZE + (M%BLOCK_SIZE? 1 : 0);
  double tempC[4] __attribute__((aligned(16))) = {0};
  double* trans __attribute__((aligned(16))) = malloc(2*temp*sizeof(double));
  bijectA(A, trans, n, n_blocks, odd);
  bijectB(B, trans, n, n_blocks, odd);

  for (bi = 0; bi < n_blocks; ++bi) {
    const int i = bi * BLOCK_SIZE;
    for (bj = 0; bj < n_blocks; ++bj) {
      const int j = bj * BLOCK_SIZE;

      kdgemm2P2(n, tempC, &trans[bi*n*2], &trans[temp + bj*n*2]);
      C[i+j*M] += tempC[0];

      if (odd == 1) {
        // Case 1: at bottom-most row of an odd row, don't store bottom ones
        if (bi == n_blocks - 1 && bj != n_blocks -1 ) { 
          C[i+(j+1)*M] += tempC[2];
        }
        // Case 2: At right-most, need to not store the top
        else if (bj == n_blocks - 1 && bi != n_blocks - 1) {
          C[i+j*M+1] += tempC[3];
        }
        // Case 3: 
        else if (bj != n_blocks - 1 && bi != n_blocks -1) {
          C[i+j*M+1] += tempC[3];
          C[i+(j+1)*M] += tempC[2];
          C[i+(j+1)*M+1] += tempC[1];
        }

      } 
      else {
        C[i+j*M+1] += tempC[3];
        C[i+(j+1)*M] += tempC[2];
        C[i+(j+1)*M+1] += tempC[1];
      }

      memset(tempC, 0, 4*sizeof(double));
    }
  }

  free(trans);
}
