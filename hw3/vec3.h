#ifndef VEC3_H
#define VEC3_H

#include "nmmintrin.h"

/*@T
 * \subsection{Vector operations}
 * 
 * A few inline functions save us from typing the same annoying bits
 * of code repeatedly.
 *
 *@c*/

inline void vec3_set(float* result, float x, float y, float z)
{
    result[0] = x;
    result[1] = y;
    result[2] = z;
}

inline void vec3_copy(float* result, float* v)
{
    vec3_set(result, v[0], v[1], v[2]);
}

inline void vec3_diff(float* result, float* a, float* b)
{
  __m128 ar = _mm_load_ps(a);
  __m128 br = _mm_load_ps(b);
  __m128 cr = _mm_sub_ps(ar, br);
  _mm_store_ps(result, cr);
}

inline void vec3_scale(float* result, float alpha, float* v)
{
    __m128 alphar = _mm_set1_ps(alpha);
    __m128 vr = _mm_load_ps(v);
    __m128 resultr = _mm_mul_ps(vr, alphar);
    _mm_store_ps(result, resultr);
}

inline float vec3_dist2(float* a, float* b)
{
    __m128 ar = _mm_load_ps(a);
    __m128 br = _mm_load_ps(b);
    __m128 cr = _mm_sub_ps(ar, br);
    return _mm_cvtss_f32(_mm_dp_ps(cr, cr, 0x7F));
}

inline float vec3_len2(float* a)
{
    __m128 ar = _mm_load_ps(a);
    return _mm_cvtss_f32(_mm_dp_ps(ar, ar, 0x7F));
}

inline void vec3_saxpy(float* result, float alpha, float* v)
{
    __m128 resultr = _mm_load_ps(result);
    __m128 alphar = _mm_set1_ps(alpha);
    __m128 vr = _mm_load_ps(v);
    resultr = _mm_add_ps(resultr, _mm_mul_ps(vr, alphar));
    _mm_store_ps(result, resultr);
}

inline void vec3_scalev(float* result, float alpha)
{
    __m128 resultr = _mm_load_ps(result);
    __m128 alphar = _mm_set1_ps(alpha);
    resultr = _mm_mul_ps(resultr, alphar);
    _mm_store_ps(result, resultr);
}

/*@q*/
#endif /* VEC3_H */
