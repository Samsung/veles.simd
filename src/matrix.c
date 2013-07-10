/*! @file matrix.c
 *  @brief New file description.
 *  @author Markovtsev Vadim <v.markovtsev@samsung.com>
 *  @version 1.0
 *
 *  @section Notes
 *  This code partially conforms to <a href="http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml">Google C++ Style Guide</a>.
 *
 *  @section Copyright
 *  Copyright 2013 Samsung R&D Institute Russia
 */

#define LIBSIMD_IMPLEMENTATION
#include "inc/simd/matrix.h"
#include <assert.h>
#include "inc/simd/memory.h"
#ifdef __AVX__
#include <immintrin.h>
#elif defined(__ARM_NEON__)
#include <arm_neon.h>
#endif

static void matrix_multiply_novec(const float *m1, const float *m2,
                                  size_t w1, size_t h1, size_t w2,
                                  size_t h2 UNUSED, float *res) {
  for (int i = 0; i < (int)w2; i++) {
    for (int j = 0; j < (int)h1; j++) {
      float sum = 0;
      for (int k = 0; k < (int)w1; k++) {
        sum += m1[j * w1 + k] * m2[k * w2 + i];
      }
      res[j * w2 + i] = sum;
    }
  }
}

#ifdef __ARM_NEON__
static void matrix_multiply_neon(const float *m1, const float *m2,
                                 size_t w1, size_t h1, size_t w2, size_t h2,
                                 float *res) {

}
#endif

#ifdef __AVX__
static void matrix_multiply_avx(const float *m1, const float *m2,
                                size_t w1, size_t h1, size_t w2,
                                size_t h2 UNUSED, float *res) {
  assert(w1 % 8 == 0);
  assert(align_complement_f32(m1) == 0);
  float col2[w1] __attribute__((aligned(32)));
  for (int i = 0; i < (int)w2; i++) {
    for (int k = 0; k < (int)w1; k++) {
      col2[k] = m2[k * w2 + i];
    }
    for (int j = 0; j < (int)h1; j++) {
      __m256 sum = _mm256_setzero_ps();
      for (int k = 0; k < (int)w1; k += 8) {
        __m256 v1 = _mm256_load_ps(m1 + j * w1 + k);
        __m256 v2 = _mm256_load_ps(col2 + k);
        __m256 dp = _mm256_dp_ps(v1, v2, 0xFF);
        sum = _mm256_add_ps(sum, dp);
      }
      res[j * w2 + i] = sum[0] + sum[4];
    }
  }
}
#endif

void matrix_multiply(int simd, const float *m1, const float *m2,
                     size_t w1, size_t h1, size_t w2, size_t h2,
                     float *res) {
  assert(w1 == h2);
  assert(m1);
  assert(m2);
  assert(res);
  assert(w1 > 0);
  assert(h1 > 0);
  assert(w2 > 0);
  if (!simd) {
    matrix_multiply_novec(m1, m2, w1, h1, w2, h2, res);
  } else {
#ifdef __ARM_NEON__
    matrix_multiply_neon(m1, m2, w1, h1, w2, h2, res);
#elif defined(__AVX__)
    matrix_multiply_avx(m1, m2, w1, h1, w2, h2, res);
#else
#error SIMD version of this function is not implemented for the selected target
#endif
  }
}

