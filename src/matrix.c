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

static void matrix_add_novec(const float *m1, const float *m2,
                      size_t w, size_t h, float *res) {
  int size = (int)w * (int)h;
  for (int i = 0; i < size; i++) {
    res[i] = m1[i] + m2[i];
  }
}

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

static void matrix_multiply_transposed_novec(const float *m1, const float *m2,
                                             size_t w1, size_t h1,
                                             size_t w2 UNUSED, size_t h2,
                                             float *res) {
  for (int j = 0; j < (int)h1; j++) {
    for (int i = 0; i < (int)h2; i++) {
      float sum = 0;
      for (int k = 0; k < (int)w1; k++) {
        sum += m1[j * w1 + k] * m2[i * w1 + k];
      }
      res[j * h2 + i] = sum;
    }
  }
}

#ifdef __ARM_NEON__
// TODO(e.sanches): Implement ARM NEON version
static void matrix_add_neon(const float *m1, const float *m2,
                            size_t w, size_t h, float *res) {
  matrix_add_novec(m1, m2, w, h, res);
}

static void matrix_multiply_neon(const float *m1, const float *m2,
                                 size_t w1, size_t h1, size_t w2,
                                 size_t h2 UNUSED, float *res) {
  float col2[w1] __attribute__((aligned(64)));
  for (int i = 0; i < (int)w2; i++) {
    for (int k = 0; k < (int)w1; k++) {
      col2[k] = m2[k * w2 + i];
    }
    for (int j = 0; j < (int)h1; j++) {
      float32x4_t sum = vdupq_n_f32(0.f);
      for (int k = 0; k < (int)w1 - 7; k += 8) {
        float32x4_t v1 = vld1q_f32(m1 + j * w1 + k);
        float32x4_t v2 = vld1q_f32(col2 + k);
        sum = vmlaq_f32(sum, v1, v2);
        v1 = vld1q_f32(m1 + j * w1 + k + 4);
        v2 = vld1q_f32(col2 + k + 4);
        sum = vmlaq_f32(sum, v1, v2);
      }
      float32x2_t sum2 = vpadd_f32(vget_high_f32(sum),
                                   vget_low_f32(sum));
      float rsum = vget_lane_f32(sum2, 0) + vget_lane_f32(sum2, 1);
      for (int k = (w1 & ~0x7); k < (int)w1; k++) {
        rsum += m1[j * w1 + k] * col2[k];
      }
      res[j * w2 + i] = rsum;
    }
  }
}

static void matrix_multiply_transposed_neon(const float *m1, const float *m2,
                                            size_t w1, size_t h1,
                                            size_t w2 UNUSED,  size_t h2,
                                            float *res) {
  for (int j = 0; j < (int)h1; j++) {
    for (int i = 0; i < (int)h2; i++) {
      float32x4_t sum = vdupq_n_f32(0.f);
      for (int k = 0; k < (int)w1 - 7; k += 8) {
        float32x4_t v1 = vld1q_f32(m1 + j * w1 + k);
        float32x4_t v2 = vld1q_f32(m2 + i * w1 + k);
        sum = vmlaq_f32(sum, v1, v2);
        v1 = vld1q_f32(m1 + j * w1 + k + 4);
        v2 = vld1q_f32(m2 + i * w1 + k + 4);
        sum = vmlaq_f32(sum, v1, v2);
      }
      float32x2_t sum2 = vpadd_f32(vget_high_f32(sum),
                                   vget_low_f32(sum));
      float rsum = vget_lane_f32(sum2, 0) + vget_lane_f32(sum2, 1);
      for (int k = (w1 & ~0x7); k < (int)w1; k++) {
        rsum += m1[j * w1 + k] * m2[i * w1 + k];
      }
      res[j * h2 + i] = rsum;
    }
  }
}
#endif

#ifdef __AVX__
// TODO(e.sanches): Implement AVX version
static void matrix_add_avx(const float *m1, const float *m2,
                           size_t w, size_t h, float *res) {
  matrix_add_novec(m1, m2, w, h, res);
}

static void matrix_multiply_avx(const float *m1, const float *m2,
                                size_t w1, size_t h1, size_t w2,
                                size_t h2 UNUSED, float *res) {
  assert(align_complement_f32(m1) == 0);
  float col2[w1] __attribute__((aligned(64)));
  for (int i = 0; i < (int)w2; i++) {
    for (int k = 0; k < (int)w1; k++) {
      col2[k] = m2[k * w2 + i];
    }
    for (int j = 0; j < (int)h1; j++) {
      __m256 sum = _mm256_setzero_ps();
      for (int k = 0; k < (int)w1 - 7; k += 8) {
        __m256 v1 = _mm256_loadu_ps(m1 + j * w1 + k);
        __m256 v2 = _mm256_load_ps(col2 + k);
        __m256 dp = _mm256_mul_ps(v1, v2);
        sum = _mm256_add_ps(sum, dp);
      }
      sum = _mm256_hadd_ps(sum, sum);
      sum = _mm256_hadd_ps(sum, sum);
      float rsum = sum[0] + sum[4];
      for (int k = (w1 & ~0x7); k < (int)w1; k++) {
        rsum += m1[j * w1 + k] * col2[k];
      }
      res[j * w2 + i] = rsum;
    }
  }
}

static void matrix_multiply_transposed_avx(const float *m1, const float *m2,
                                           size_t w1, size_t h1,
                                           size_t w2 UNUSED,  size_t h2,
                                           float *res) {
  assert(align_complement_f32(m1) == 0);
  assert(align_complement_f32(m2) == 0);
  for (int j = 0; j < (int)h1; j++) {
    for (int i = 0; i < (int)h2; i++) {
      __m256 sum = _mm256_setzero_ps();
      for (int k = 0; k < (int)w1 - 7; k += 8) {
        __m256 v1 = _mm256_loadu_ps(m1 + j * w1 + k);
        __m256 v2 = _mm256_loadu_ps(m2 + i * w1 + k);
        __m256 dp = _mm256_mul_ps(v1, v2);
        sum = _mm256_add_ps(sum, dp);
      }
      sum = _mm256_hadd_ps(sum, sum);
      sum = _mm256_hadd_ps(sum, sum);
      float rsum = sum[0] + sum[4];
      for (int k = (w1 & ~0x7); k < (int)w1; k++) {
        rsum += m1[j * w1 + k] * m2[i * w1 + k];
      }
      res[j * h2 + i] = rsum;
    }
  }
}
#endif

void matrix_add(int simd, const float *m1, const float *m2,
                size_t w, size_t h, float *res) {
  assert(m1);
  assert(m2);
  assert(res);
  assert(w > 0);
  assert(h > 0);
  if (simd) {
#ifdef __ARM_NEON__
    matrix_add_neon(m1, m2, w, h, res);
  } else {
#elif defined(__AVX__)
    matrix_add_avx(m1, m2, w, h, res);
  } else {
#else
  } {
#endif
    matrix_add_novec(m1, m2, w, h, res);
  }
}

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
  if (simd) {
#ifdef __ARM_NEON__
    matrix_multiply_neon(m1, m2, w1, h1, w2, h2, res);
  } else {
#elif defined(__AVX__)
    matrix_multiply_avx(m1, m2, w1, h1, w2, h2, res);
  } else {
#else
  } {
#endif    
    matrix_multiply_novec(m1, m2, w1, h1, w2, h2, res);
  }
}

void matrix_multiply_transposed(int simd, const float *m1, const float *m2,
                                size_t w1, size_t h1, size_t w2, size_t h2,
                                float *res) {
  assert(w1 == w2);
  assert(m1);
  assert(m2);
  assert(res);
  assert(w1 > 0);
  assert(h1 > 0);
  assert(h2 > 0);
  if (simd) {
#ifdef __ARM_NEON__
    matrix_multiply_transposed_neon(m1, m2, w1, h1, w2, h2, res);
  } else {
#elif defined(__AVX__)
    matrix_multiply_transposed_avx(m1, m2, w1, h1, w2, h2, res);
  } else {
#else
  } {
#endif
    matrix_multiply_transposed_novec(m1, m2, w1, h1, w2, h2, res);
  }
}

