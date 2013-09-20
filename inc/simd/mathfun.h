/*! @file mathfun.h
 *  @brief SIMD implementation of sin, cos, sincos, exp and log
 *  @author Ernesto Sanches <ernestosanches@gmail.com>
 *  @version 1.0
 *
 *  @section Notes
 *  This code partially conforms to <a href="http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml">Google C++ Style Guide</a>.
 *
 *  @section Copyright
 *  Copyright 2013 Samsung R&D Institute Russia
 */

#ifndef INC_SIMD_MATHFUN_H_
#define INC_SIMD_MATHFUN_H_

#include <math.h>
#include <stddef.h>
#include <simd/common.h>
#include <simd/attributes.h>

SIMD_API_BEGIN

typedef float (*PsvStdFunc)(float);

inline void func_psv_novec(PsvStdFunc func, const float *src,
                           size_t length, float *res) {
  for(size_t i = 0; i < length; ++i) {
    res[i] = func(src[i]);
  }
}

inline void sin_psv_novec(const float *src, size_t length, float *res) {
  func_psv_novec(sinf, src, length, res);
}

inline void cos_psv_novec(const float *src, size_t length, float *res) {
  func_psv_novec(cosf, src, length, res);
}

inline void log_psv_novec(const float *src, size_t length, float *res) {
  func_psv_novec(logf, src, length, res);
}

inline void exp_psv_novec(const float *src, size_t length, float *res) {
  func_psv_novec(expf, src, length, res);
}

#ifdef __ARM_NEON__
#include <simd/neon_mathfun.h>  // NO_LINT

typedef float32x4_t (*PsvNeonFunc)(float32x4_t);

inline void func_psv_neon(PsvNeonFunc neon_func, PsvStdFunc std_func,
                          const float *src, size_t length, float *res
                          ) {
  int ilength = (int)length;
  int i = 0;
  for (; i < ilength - 3; i += 4) {
    float32x4_t vec = vld1q_f32(src + i);
    float32x4_t vec_r = neon_func(vec);
    vst1q_f32(res + i, vec_r);
  }
  for (; i < ilength; i++) {
    res[i] = std_func(src[i]);
  }
}

inline void sin_psv_neon(const float *src, size_t length, float *res) {
  func_psv_neon(sin_ps, sinf, src, length, res);
}

inline void cos_psv_neon(const float *src, size_t length, float *res) {
  func_psv_neon(cos_ps, cosf, src, length, res);
}

inline void log_psv_neon(const float *src, size_t length, float *res) {
  func_psv_neon(log_ps, logf, src, length, res);
}

inline void exp_psv_neon(const float *src, size_t length, float *res) {
  func_psv_neon(exp_ps, expf, src, length, res);
}

#endif

#ifdef __AVX__
#include <simd/avx_mathfun.h>  // NO_LINT

typedef __m256 (*PsvAvxFunc)(__m256);

inline void func_psv_avx(PsvAvxFunc avx_func, PsvStdFunc std_func,
                         const float *src, size_t length, float *res) {
  int ilength = (int)length;
  int i = 0;
  for (; i < ilength - 7; i += 8) {
    __m256 vec = _mm256_loadu_ps(src + i);
    __m256 vec_r = avx_func(vec);
    _mm256_storeu_ps(res + i, vec_r);
  }
  for (; i < ilength; i++) {
    res[i] = std_func(src[i]);
  }
}

inline void sin_psv_avx(const float *src, size_t length, float *res) {
  func_psv_avx(sin256_ps, sinf, src, length, res);
}

inline void cos_psv_avx(const float *src, size_t length, float *res) {
  func_psv_avx(cos256_ps, cosf, src, length, res);
}

inline void log_psv_avx(const float *src, size_t length, float *res) {
  func_psv_avx(log256_ps, logf, src, length, res);
}

inline void exp_psv_avx(const float *src, size_t length, float *res) {
  func_psv_avx(exp256_ps, expf, src, length, res);
}

#endif

INLINE NOTNULL(2, 4) void sin_psv(int simd, const float *src, size_t length,
                                  float *res) {
  if (simd) {
#ifdef __ARM_NEON__
    sin_psv_neon(src, length, res);
  } else {
#elif defined(__AVX__)
    sin_psv_avx(src, length, res);
  } else {
#else
  } {
#endif
    sin_psv_novec(src, length, res);
  }
}

INLINE NOTNULL(2, 4) void cos_psv(int simd, const float *src, size_t length,
                                  float *res) {
  if (simd) {
#ifdef __ARM_NEON__
    cos_psv_neon(src, length, res);
  } else {
#elif defined(__AVX__)
    cos_psv_avx(src, length, res);
  } else {
#else
  } {
#endif
    cos_psv_novec(src, length, res);
  }
}

INLINE NOTNULL(2, 4) void log_psv(int simd, const float *src, size_t length,
                                  float *res) {
  if (simd) {
#ifdef __ARM_NEON__
    log_psv_neon(src, length, res);
  } else {
#elif defined(__AVX__)
    log_psv_avx(src, length, res);
  } else {
#else
  } {
#endif
    log_psv_novec(src, length, res);
  }
}

INLINE NOTNULL(2, 4) void exp_psv(int simd, const float *src, size_t length,
                                  float *res) {
  if (simd) {
#ifdef __ARM_NEON__
    exp_psv_neon(src, length, res);
  } else {
#elif defined(__AVX__)
    exp_psv_avx(src, length, res);
  } else {
#else
  } {
#endif
    exp_psv_novec(src, length, res);
  }
}

//////////////////////////

SIMD_API_END

#endif  // INC_SIMD_MATHFUN_H_
