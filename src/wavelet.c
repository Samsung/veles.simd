/*! @file wavelet.c
 *  @brief Wavelet low level functions.
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
#include "inc/simd/wavelet.h"
#include <assert.h>
#include <string.h>
#include "inc/simd/arithmetic-inl.h"
#define WAVELET_INTERNAL_USE
#include "src/coiflets.h"
#include "src/daubechies.h"
#include "src/symlets.h"
#include "inc/simd/memory.h"

#define max(a,b) \
   ({ \
     __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; \
   })

INLINE void check_length(size_t length) {
  assert(length > 0);
  assert(length % 2 == 0);
}

INLINE size_t aligned_length(size_t length) {
  --length;
  length |= length >> 1;
  length |= length >> 2;
  length |= length >> 4;
  length |= length >> 8;
  length |= length >> 16;
  return length + 1;
}

INLINE NOTNULL(2, 4) void wavelet_prepare_array_memcpy(
    int order, const float *src, size_t length, float *res) {
  check_length(length);

  if (res != src) {
    memcpy(res, src, length * sizeof(float));
  }

  if (length > 8) {
    size_t alength = aligned_length(length);
    int copySize = (alength - 8) * sizeof(float);
    memcpy(res + alength,          src + 2, copySize);
    if (order > 4) {
      memcpy(res + alength * 2 -  8, src + 4, copySize);
      memcpy(res + alength * 3 - 16, src + 6, copySize);
    }
  }
}

int wavelet_validate_order(WaveletType type, int order) {
  size_t uorder = (size_t)order;
  switch (type) {
    case WAVELET_TYPE_DAUBECHIES:
      return (uorder <= sizeof(kDaubechiesF[0]) / sizeof(kDaubechiesF[0][0])) &&
             (uorder % 2 == 0);
    case WAVELET_TYPE_COIFLET:
      return (uorder <= sizeof(kCoifletsF[0]) / sizeof(kCoifletsF[0][0])) &&
             (uorder % 6 == 0);
    case WAVELET_TYPE_SYMLET:
      return (uorder <= sizeof(kSymletsF[0]) / sizeof(kSymletsF[0][0])) &&
             (uorder % 2 == 0);
    default:
      return 0;
  }
}

float *wavelet_prepare_array(int order
#ifndef __AVX__
                            UNUSED
#endif
                            , const float *src, size_t length
#ifndef __AVX__
                             UNUSED
#endif
) {
  check_length(length);
#ifndef __AVX__
  float *res = mallocf(length);
  memcpy(res, src, length * sizeof(src[0]));
#else
  size_t alength = aligned_length(length);
  float *res = mallocf(alength * (order > 4? 4 : 2));
  wavelet_prepare_array_memcpy(order, src, length, res);
#endif
  return res;
}

float *wavelet_allocate_destination(int order
#ifndef __AVX__
                            UNUSED
#endif
                            , size_t sourceLength) {
  check_length(sourceLength);
  assert(sourceLength % 4 == 0);

#ifndef __AVX__
  float *res = mallocf(sourceLength / 2);
#else
  size_t alength = aligned_length(sourceLength);
  float *res = mallocf(alength * (order > 4? 2 : 1));
#endif
  return res;
}

void wavelet_recycle_source(int order
#ifndef __AVX__
                            UNUSED
#endif
                            , float *src, size_t length,
                            float **desthihi, float **desthilo,
                            float **destlohi, float **destlolo) {
  if (length == 0 || length % 4 != 0) {
    *desthihi = NULL;
    *desthilo = NULL;
    *destlohi = NULL;
    *destlolo = NULL;
    return;
  }

#ifndef __AVX__
  int lq = length / 4;
#else
  int lq = aligned_length(length) / (order > 4? 1 : 2);
  if (length < 16) {
    lq = length / 4;
  }
#endif
  *desthihi = src;
  *desthilo = src + lq;
  *destlohi = src + lq * 2;
  *destlolo = src + lq * 3;
}

INLINE void check_wavelet_order(WaveletType type, size_t order) {
  switch (type) {
    case WAVELET_TYPE_DAUBECHIES:
      assert(order <= sizeof(kDaubechiesF[0]) / sizeof(kDaubechiesF[0][0]) &&
             order % 2 == 0 &&
             "Supported Daubechies orders are 2..76 (even numbers only)");
      break;
    case WAVELET_TYPE_COIFLET:
      assert(order <= sizeof(kCoifletsF[0]) / sizeof(kCoifletsF[0][0]) &&
             order % 6 == 0 &&
             "Supported Coiflet orders are 6, 12, 18, 24 and 30");
      break;
    case WAVELET_TYPE_SYMLET:
      assert(order <= sizeof(kSymletsF[0]) / sizeof(kSymletsF[0][0]) &&
             order % 2 == 0 &&
             "Supported Daubechies orders are 2..76 (even numbers only)");
      break;
  }
}

INLINE NOTNULL(3, 4) void initialize_highpass_lowpass(
    WaveletType type, int order, float *restrict highpass,
    float *restrict lowpass) {
  assert(order >= 2);
  check_wavelet_order(type, order);
  for (int i = 0; i < order; i++) {
    float val = 0.f;
    switch (type) {
      case WAVELET_TYPE_DAUBECHIES:
        val =  kDaubechiesF[order / 2 - 1][i];
        break;
      case WAVELET_TYPE_COIFLET:
        assert(order >= 6);
        val = kCoifletsF[order / 6 - 1][i];
        break;
      case WAVELET_TYPE_SYMLET:
        val = kSymletsF[order / 2 - 1][i];
        break;
    }
    lowpass[i] = val;
    highpass[order - i - 1] = (i & 1) ? val : -val;
  }
}

INLINE NOTNULL(4, 5) void stationary_initialize_highpass_lowpass(
    WaveletType type, int size, int level, float *restrict highpass,
    float *restrict lowpass) {
  int stride = 1 << (level - 1);
  if (stride == 1) {
    initialize_highpass_lowpass(type, size, highpass, lowpass);
    return;
  }
  assert(size % stride == 0);
  int order = size / stride;
  assert(order >= 2);
  check_wavelet_order(type, order);
  for (int i = 0; i < size; i++) {
    if (i % stride != 0) {
      lowpass[i] = highpass[size - i] = 0;
      continue;
    }
    int ri = i / stride;
    float val = 0.f;
    switch (type) {
      case WAVELET_TYPE_DAUBECHIES:
        val = kDaubechiesF[order / 2 - 1][ri];
        break;
      case WAVELET_TYPE_COIFLET:
        assert(order >= 6);
        val = kCoifletsF[order / 6 - 1][ri];
        break;
      case WAVELET_TYPE_SYMLET:
        val = kSymletsF[order / 2 - 1][ri];
        break;
    }
    lowpass[i] = val;
    highpass[size - i - stride] = (ri & 1) ? val : -val;
  }
}

INLINE NOTNULL(3, 5) void initialize_extension(ExtensionType ext,
                                               size_t extLength,
                                               const float *restrict src,
                                               size_t length,
                                               float *restrict result) {
  for (int i = 0; i < (int)extLength; i++) {
    switch (ext) {
      case EXTENSION_TYPE_PERIODIC:
        result[i] = src[i % length];
        break;
      case EXTENSION_TYPE_MIRROR:
        result[i] = src[length - 1 - (i % length)];
        break;
      case EXTENSION_TYPE_CONSTANT:
        result[i] = src[length - 1];
        break;
      case EXTENSION_TYPE_ZERO:
        result[i] = 0;
        break;
    }
  }
}

void wavelet_apply_na(WaveletType type, int order, ExtensionType ext,
                      const float *__restrict src, size_t length,
                      float *__restrict desthi, float *__restrict destlo) {
  check_length(length);
  assert(src && desthi && destlo);

  int ilength = (int)length;
  float highpassC[order], lowpassC[order];
  initialize_highpass_lowpass(type, order, highpassC, lowpassC);
  float src_ext[order];
  initialize_extension(ext, order, src, length, src_ext);

  if (ilength != order) {
    for (int i = 0, di = 0; i < ilength; i += 2, di++) {
      float reshi = 0.f, reslo = 0.f;
      for (int j = 0; j < order; j++) {
        int index = i + j;
        float srcval = index < ilength? src[index] : src_ext[index - ilength];
        reshi += highpassC[j] * srcval;
        reslo += lowpassC[j] * srcval;
      }
      desthi[di] = reshi;
      destlo[di] = reslo;
    }
  } else {
    if (order == 8) {
      // Give a chance to the compiler to optimize these two loops
      for (int i = 0, di = 0; i < 8; i += 2, di++) {
        float reshi = 0.f, reslo = 0.f;
        for (int j = 0; j < 8; j++) {
          int index = i + j;
          float srcval = index < 8? src[index] : src_ext[index - 8];
          reshi += highpassC[j] * srcval;
          reslo += lowpassC[j] * srcval;
        }
        desthi[di] = reshi;
        destlo[di] = reslo;
      }
    } else {
      for (int i = 0, di = 0; i < order; i += 2, di++) {
        float reshi = 0.f, reslo = 0.f;
        for (int j = 0; j < order; j++) {
          int index = i + j;
          float srcval = index < order? src[index] : src_ext[index - order];
          reshi += highpassC[j] * srcval;
          reslo += lowpassC[j] * srcval;
        }
        desthi[di] = reshi;
        destlo[di] = reslo;
      }
    }
  }
}

void stationary_wavelet_apply_na(WaveletType type, int order, int level,
                                 ExtensionType ext,
                                 const float *__restrict src, size_t length,
                                 float *__restrict desthi,
                                 float *__restrict destlo) {
  assert(length > 0);
  assert(src && desthi && destlo);

  int size = order * (1 << (level - 1));
  int ilength = (int)length;
  float highpassC[size], lowpassC[size];
  stationary_initialize_highpass_lowpass(type, size, level, highpassC,
                                         lowpassC);
  float src_ext[size];
  initialize_extension(ext, size, src, length, src_ext);

  int stride = 1 << (level - 1);
  if (ilength != size) {
    for (int i = 0; i < ilength; i++) {
      float reshi = 0.f, reslo = 0.f;
      for (int j = 0; j < size; j += stride) {
        int index = i + j;
        float srcval = index < ilength? src[index] : src_ext[index - ilength];
        reshi += highpassC[j] * srcval;
        reslo += lowpassC[j] * srcval;
      }
      desthi[i] = reshi;
      destlo[i] = reslo;
    }
  } else {
    if (size == 8) {
      // Give a chance to the compiler to optimize these two loops
      for (int i = 0; i < 8; i++) {
        float reshi = 0.f, reslo = 0.f;
        for (int j = 0; j < 8; j += stride) {
          int index = i + j;
          float srcval = index < 8? src[index] : src_ext[index - 8];
          reshi += highpassC[j] * srcval;
          reslo += lowpassC[j] * srcval;
        }
        desthi[i] = reshi;
        destlo[i] = reslo;
      }
    } else {
      for (int i = 0; i < size; i++) {
        float reshi = 0.f, reslo = 0.f;
        for (int j = 0; j < size; j += stride) {
          int index = i + j;
          float srcval = index < size? src[index] : src_ext[index - size];
          reshi += highpassC[j] * srcval;
          reslo += lowpassC[j] * srcval;
        }
        desthi[i] = reshi;
        destlo[i] = reslo;
      }
    }
  }
}

#ifdef __AVX__
#define ALIGN_ORDER(order) (order % 8 == 0? order : (order / 8 + 1) * 8)
#else
#define ALIGN_ORDER(order) (order % 4 == 0? order : (order / 4 + 1) * 4)
#define align_complement_f32(x) 0
#endif

#define DECLARE_PASSC(order) \
  float highpassC[ALIGN_ORDER(order)] __attribute__ ((aligned (64))), \
        lowpassC[ALIGN_ORDER(order)] __attribute__ ((aligned (64)))

static void wavelet_apply2(WaveletType type, ExtensionType ext,
                           const float *__restrict src, size_t length,
                           float *__restrict desthi,
                           float *__restrict destlo) {
#ifdef SIMD
  check_length(length);
  assert(src && desthi && destlo);

  if (align_complement_f32(src) != 0 || length < 8) {
    wavelet_apply_na(type, 2, ext, src, length, desthi, destlo);
    return;
  }

  int ilength = (int)length;
  DECLARE_PASSC(2);
  initialize_highpass_lowpass(type, 2, highpassC, lowpassC);
  float src_ext[2];
  initialize_extension(ext, 2, src, length, src_ext);

#ifdef __AVX__
  int simd_end = ilength - 14;
  highpassC[2] = highpassC[0];
  highpassC[3] = highpassC[1];
  lowpassC[2] = lowpassC[0];
  lowpassC[3] = lowpassC[1];
  memcpy(&highpassC[4], &highpassC[0], sizeof(highpassC) / 2);
  memcpy(&lowpassC[4], &lowpassC[0], sizeof(lowpassC) / 2);

  const __m256 hpvec = _mm256_load_ps(highpassC);
  const __m256 lpvec = _mm256_load_ps(lowpassC);
  for (int i = 0; i < simd_end; i += 16) {
    __m256 srcvec1 = _mm256_load_ps(src + i);
    __m256 srcvec2 = _mm256_load_ps(src + i + 8);
    __m256 vechi1 = _mm256_mul_ps(srcvec1, hpvec);
    __m256 vechi2 = _mm256_mul_ps(srcvec2, hpvec);
    __m256 vechi1p = _mm256_permute2f128_ps(vechi1, vechi2, 32);
    __m256 vechi2p = _mm256_permute2f128_ps(vechi1, vechi2, 49);
    __m256 vechi = _mm256_hadd_ps(vechi1p, vechi2p);
    __m256 veclo1 = _mm256_mul_ps(srcvec1, lpvec);
    __m256 veclo2 = _mm256_mul_ps(srcvec2, lpvec);
    __m256 veclo1p = _mm256_permute2f128_ps(veclo1, veclo2, 32);
    __m256 veclo2p = _mm256_permute2f128_ps(veclo1, veclo2, 49);
    __m256 veclo = _mm256_hadd_ps(veclo1p, veclo2p);
    _mm256_store_ps(desthi + i / 2, vechi);
    _mm256_store_ps(destlo + i / 2, veclo);
  }
#elif defined(__ARM_NEON__)
  int simd_end = ilength - 6;
  highpassC[2] = highpassC[0];
  highpassC[3] = highpassC[1];
  lowpassC[2] = lowpassC[0];
  lowpassC[3] = lowpassC[1];
  const float32x4_t hivec = vld1q_f32(highpassC);
  const float32x4_t lovec = vld1q_f32(lowpassC);
  for (int i = 0; i < simd_end; i += 8) {
    float32x4_t srcvec1 = vld1q_f32(src + i);
    float32x4_t srcvec2 = vld1q_f32(src + i + 4);

    float32x4_t vechiadd1 = vmulq_f32(srcvec1, hivec);
    float32x4_t vecloadd1 = vmulq_f32(srcvec1, lovec);
    float32x4_t vechiadd2 = vmulq_f32(srcvec2, hivec);
    float32x4_t vecloadd2 = vmulq_f32(srcvec2, lovec);

    float32x2_t vecreshi1 = vpadd_f32(vget_low_f32(vechiadd1),
                                      vget_high_f32(vechiadd1));
    float32x2_t vecreslo1 = vpadd_f32(vget_low_f32(vecloadd1),
                                      vget_high_f32(vecloadd1));
    float32x2_t vecreshi2 = vpadd_f32(vget_low_f32(vechiadd2),
                                      vget_high_f32(vechiadd2));
    float32x2_t vecreslo2 = vpadd_f32(vget_low_f32(vecloadd2),
                                      vget_high_f32(vecloadd2));

    vst1_f32(desthi + i / 2, vecreshi1);
    vst1_f32(desthi + i / 2 + 2, vecreshi2);
    vst1_f32(destlo + i / 2, vecreslo1);
    vst1_f32(destlo + i / 2 + 2, vecreslo2);
  }
#else
#error This SIMD variant is not supported.
#endif  // #elif defined(__ARM_NEON__)
  // Finish with the extended end
  for (int i = simd_end, di = simd_end / 2; i < ilength; i += 2, di++) {
    float reshi = 0.f, reslo = 0.f;
    for (int j = 0; j < 2; j++) {
      int index = i + j;
      float srcval = index < ilength? src[index] : src_ext[index - ilength];
      reshi += highpassC[j] * srcval;
      reslo += lowpassC[j] * srcval;
    }
    desthi[di] = reshi;
    destlo[di] = reslo;
  }
#else  // #ifdef SIMD
  wavelet_apply_na(type, 2, ext, src, length, desthi, destlo);
#endif
}

static void stationary_wavelet_apply2(WaveletType type, int level,
                                      ExtensionType ext,
                                      const float *__restrict src,
                                      size_t length,
                                      float *__restrict desthi,
                                      float *__restrict destlo) {
  int stride = 1 << (level - 1);
#ifdef SIMD
  assert(length > 0);
  assert(src && desthi && destlo);

  if (stride >= 4 ||
#ifdef __AVX__
      length < 12
#elif defined(__ARM_NEON__)
      length < 8
#else
      0
#endif
  ) {
    stationary_wavelet_apply_na(type, 2 / stride, level, ext, src, length,
                                desthi, destlo);
    return;
  }

  int ilength = (int)length;
  DECLARE_PASSC(2);
  stationary_initialize_highpass_lowpass(type, 2, level, highpassC, lowpassC);
  float src_ext[2];
  initialize_extension(ext, 2, src, length, src_ext);

#ifdef __AVX__
  int simd_end = ilength - 8;
  highpassC[2] = highpassC[0];
  highpassC[3] = highpassC[1];
  lowpassC[2] = lowpassC[0];
  lowpassC[3] = lowpassC[1];
  memcpy(&highpassC[4], &highpassC[0], sizeof(highpassC) / 2);
  memcpy(&lowpassC[4], &lowpassC[0], sizeof(lowpassC) / 2);

  const __m256 hpvec = _mm256_load_ps(highpassC);
  const __m256 lpvec = _mm256_load_ps(lowpassC);
  for (int i = 0; i < simd_end; i += 8) {
    __m256 srcvec1 = _mm256_loadu_ps(src + i);
    __m256 srcvec2 = _mm256_loadu_ps(src + i + 1);
    __m256 vechi1 = _mm256_mul_ps(srcvec1, hpvec);
    __m256 vechi2 = _mm256_mul_ps(srcvec2, hpvec);
    __m256 vechi = _mm256_hadd_ps(vechi1, vechi2);
    vechi = _mm256_permute_ps(vechi, 216);
    __m256 veclo1 = _mm256_mul_ps(srcvec1, lpvec);
    __m256 veclo2 = _mm256_mul_ps(srcvec2, lpvec);
    __m256 veclo = _mm256_hadd_ps(veclo1, veclo2);
    veclo = _mm256_permute_ps(veclo, 216);
    _mm256_store_ps(desthi + i, vechi);
    _mm256_store_ps(destlo + i, veclo);
  }
#elif defined(__ARM_NEON__)
  int simd_end = ilength - 4;
  highpassC[2] = highpassC[0];
  highpassC[3] = highpassC[1];
  lowpassC[2] = lowpassC[0];
  lowpassC[3] = lowpassC[1];
  const float32x4_t hivec = vld1q_f32(highpassC);
  const float32x4_t lovec = vld1q_f32(lowpassC);
  for (int i = 0; i < simd_end; i += 4) {
    float32x4_t srcvec1 = vld1q_f32(src + i);
    float32x4_t srcvec2 = vld1q_f32(src + i + 1);

    float32x4_t vechiadd1 = vmulq_f32(srcvec1, hivec);
    float32x4_t vecloadd1 = vmulq_f32(srcvec1, lovec);
    float32x4_t vechiadd2 = vmulq_f32(srcvec2, hivec);
    float32x4_t vecloadd2 = vmulq_f32(srcvec2, lovec);

    float32x2_t vecreshi1 = vpadd_f32(vget_low_f32(vechiadd1),
                                      vget_high_f32(vechiadd1));
    float32x2_t vecreslo1 = vpadd_f32(vget_low_f32(vecloadd1),
                                      vget_high_f32(vecloadd1));
    float32x2_t vecreshi2 = vpadd_f32(vget_low_f32(vechiadd2),
                                      vget_high_f32(vechiadd2));
    float32x2_t vecreslo2 = vpadd_f32(vget_low_f32(vecloadd2),
                                      vget_high_f32(vecloadd2));
    float32x2x2_t rhi = vtrn_f32(vecreshi1, vecreshi2);
    float32x2x2_t rlo = vtrn_f32(vecreslo1, vecreslo2);
    vst1_f32(desthi + i, rhi.val[0]);
    vst1_f32(desthi + i + 2, rhi.val[1]);
    vst1_f32(destlo + i, rlo.val[0]);
    vst1_f32(destlo + i + 2, rlo.val[1]);
  }
#else
#error This SIMD variant is not supported.
#endif  // #elif defined(__ARM_NEON__)
  // Finish with the extended end
  for (int i = simd_end; i < ilength; i++) {
    float reshi = 0.f, reslo = 0.f;
    for (int j = 0; j < 2; j++) {
      int index = i + j;
      float srcval = index < ilength? src[index] : src_ext[index - ilength];
      reshi += highpassC[j] * srcval;
      reslo += lowpassC[j] * srcval;
    }
    desthi[i] = reshi;
    destlo[i] = reslo;
  }
#else  // #ifdef SIMD
  stationary_wavelet_apply_na(type, 2 / stride, level, ext, src, length,
                              desthi, destlo);
#endif
}

static void wavelet_apply4(WaveletType type, ExtensionType ext,
                           const float *__restrict src, size_t length,
                           float *__restrict desthi,
                           float *__restrict destlo) {
#ifdef SIMD
  check_length(length);
  assert(src && desthi && destlo);

  if (align_complement_f32(src) != 0 || length < 8) {
    wavelet_apply_na(type, 4, ext, src, length, desthi, destlo);
    return;
  }

  int ilength = (int)length;
  DECLARE_PASSC(4);
  initialize_highpass_lowpass(type, 4, highpassC, lowpassC);
  float src_ext[4];
  initialize_extension(ext, 4, src, length, src_ext);

#ifdef __AVX__
  memcpy(&highpassC[4], &highpassC[0], sizeof(highpassC) / 2);
  memcpy(&lowpassC[4], &lowpassC[0], sizeof(lowpassC) / 2);

  const __m256 hpvec = _mm256_load_ps(highpassC);
  const __m256 lpvec = _mm256_load_ps(lowpassC);
  size_t alength = aligned_length(length);
  for (int i = 0, di = 0; i < ilength - 6; i += 8, di += 4) {
    __m256 srcvec1 = _mm256_load_ps(src + i);
    __m256 srcvec2 = _mm256_load_ps(src + i + alength);
    __m256 vechi1 = _mm256_dp_ps(srcvec1, hpvec, 0xFF);
    __m256 veclo1 = _mm256_dp_ps(srcvec1, lpvec, 0xFF);
    __m256 vechi2 = _mm256_dp_ps(srcvec2, hpvec, 0xFF);
    __m256 veclo2 = _mm256_dp_ps(srcvec2, lpvec, 0xFF);
    desthi[di] = vechi1[0];
    destlo[di] = veclo1[0];
    desthi[di + 1] = vechi2[0];
    destlo[di + 1] = veclo2[0];
    desthi[di + 2] = vechi1[4];
    destlo[di + 2] = veclo1[4];
    desthi[di + 3] = vechi2[4];
    destlo[di + 3] = veclo2[4];
  }
#elif defined(__ARM_NEON__)
  const float32x4_t hivec = vld1q_f32(highpassC);
  const float32x4_t lovec = vld1q_f32(lowpassC);
  for (int i = 0, di = 0; i < ilength - 6; i += 4, di += 2) {
    float32x4_t srcvec1 = vld1q_f32(src + i);
    float32x4_t srcvec2 = vld1q_f32(src + i + 2);

    float32x4_t vechiadd1 = vmulq_f32(srcvec1, hivec);
    float32x4_t vecloadd1 = vmulq_f32(srcvec1, lovec);
    float32x4_t vechiadd2 = vmulq_f32(srcvec2, hivec);
    float32x4_t vecloadd2 = vmulq_f32(srcvec2, lovec);

    float32x2_t vechipair1 = vadd_f32(vget_high_f32(vechiadd1),
                                      vget_low_f32(vechiadd1));
    float32x2_t veclopair1 = vadd_f32(vget_high_f32(vecloadd1),
                                      vget_low_f32(vecloadd1));
    float32x2_t vechipair2 = vadd_f32(vget_high_f32(vechiadd2),
                                      vget_low_f32(vechiadd2));
    float32x2_t veclopair2 = vadd_f32(vget_high_f32(vecloadd2),
                                      vget_low_f32(vecloadd2));

    float32x2_t vecres1 = vpadd_f32(vechipair1, veclopair1);
    float32x2_t vecres2 = vpadd_f32(vechipair2, veclopair2);

    desthi[di] = vget_lane_f32(vecres1, 0);
    destlo[di] = vget_lane_f32(vecres1, 1);
    desthi[di + 1] = vget_lane_f32(vecres2, 0);
    destlo[di + 1] = vget_lane_f32(vecres2, 1);
  }
#else
#error This SIMD variant is not supported.
#endif  // #elif defined(__ARM_NEON__)
  // Finish with the extended end
  for (int i = ilength - 6, di = (ilength - 6) / 2;
       i < ilength; i += 2, di++) {
    float reshi = 0.f, reslo = 0.f;
    for (int j = 0; j < 4; j++) {
      int index = i + j;
      float srcval = index < ilength? src[index] : src_ext[index - ilength];
      reshi += highpassC[j] * srcval;
      reslo += lowpassC[j] * srcval;
    }
    desthi[di] = reshi;
    destlo[di] = reslo;
  }
#ifdef __AVX__
  if (length >= 16) {
    wavelet_prepare_array_memcpy(4, desthi, length / 2, desthi);
    wavelet_prepare_array_memcpy(4, destlo, length / 2, destlo);
  }
#endif  // #ifdef __AVX__
#else  // #ifdef SIMD
  wavelet_apply_na(type, 4, ext, src, length, desthi, destlo);
#endif
}

static void stationary_wavelet_apply4(WaveletType type, int level,
                                      ExtensionType ext,
                                      const float *__restrict src,
                                      size_t length,
                                      float *__restrict desthi,
                                      float *__restrict destlo) {
  int stride = 1 << (level - 1);
#ifdef SIMD
  assert(length > 0);
  assert(src && desthi && destlo);

  if (stride >= 4 ||
#ifdef __AVX__
      length < 12
#elif defined(__ARM_NEON__)
      length < 8
#else
      0
#endif
  ) {
    stationary_wavelet_apply_na(type, 4 / stride, level, ext, src, length,
                                desthi, destlo);
    return;
  }

  int ilength = (int)length;
  DECLARE_PASSC(4);
  stationary_initialize_highpass_lowpass(type, 4, level, highpassC, lowpassC);
  float src_ext[4];
  initialize_extension(ext, 4, src, length, src_ext);

#ifdef __AVX__
  int simd_end = ilength - 10;
  memcpy(&highpassC[4], &highpassC[0], sizeof(highpassC) / 2);
  memcpy(&lowpassC[4], &lowpassC[0], sizeof(lowpassC) / 2);

  const __m256 hpvec = _mm256_load_ps(highpassC);
  const __m256 lpvec = _mm256_load_ps(lowpassC);
  for (int i = 0; i < simd_end; i += 8) {
    __m256 srcvec1 = _mm256_loadu_ps(src + i);
    __m256 srcvec2 = _mm256_loadu_ps(src + i + 1);
    __m256 srcvec3 = _mm256_loadu_ps(src + i + 2);
    __m256 srcvec4 = _mm256_loadu_ps(src + i + 3);
    __m256 vechi1 = _mm256_dp_ps(srcvec1, hpvec, 0xFF);
    __m256 veclo1 = _mm256_dp_ps(srcvec1, lpvec, 0xFF);
    __m256 vechi2 = _mm256_dp_ps(srcvec2, hpvec, 0xFF);
    __m256 veclo2 = _mm256_dp_ps(srcvec2, lpvec, 0xFF);
    __m256 vechi3 = _mm256_dp_ps(srcvec3, hpvec, 0xFF);
    __m256 veclo3 = _mm256_dp_ps(srcvec3, lpvec, 0xFF);
    __m256 vechi4 = _mm256_dp_ps(srcvec4, hpvec, 0xFF);
    __m256 veclo4 = _mm256_dp_ps(srcvec4, lpvec, 0xFF);

    desthi[i] = vechi1[0];
    destlo[i] = veclo1[0];
    desthi[i + 1] = vechi2[0];
    destlo[i + 1] = veclo2[0];
    desthi[i + 2] = vechi3[0];
    destlo[i + 2] = veclo3[0];
    desthi[i + 3] = vechi4[0];
    destlo[i + 3] = veclo4[0];

    desthi[i + 4] = vechi1[4];
    destlo[i + 4] = veclo1[4];
    desthi[i + 5] = vechi2[4];
    destlo[i + 5] = veclo2[4];
    desthi[i + 6] = vechi3[4];
    destlo[i + 6] = veclo3[4];
    desthi[i + 7] = vechi4[4];
    destlo[i + 7] = veclo4[4];
  }
#elif defined(__ARM_NEON__)
  int simd_end = ilength - 6;
  const float32x4_t hivec = vld1q_f32(highpassC);
  const float32x4_t lovec = vld1q_f32(lowpassC);
  for (int i = 0; i < simd_end; i += 4) {
    float32x4_t srcvec1 = vld1q_f32(src + i);
    float32x4_t srcvec2 = vld1q_f32(src + i + 1);
    float32x4_t srcvec3 = vld1q_f32(src + i + 2);
    float32x4_t srcvec4 = vld1q_f32(src + i + 3);

    float32x4_t vechiadd1 = vmulq_f32(srcvec1, hivec);
    float32x4_t vecloadd1 = vmulq_f32(srcvec1, lovec);
    float32x4_t vechiadd2 = vmulq_f32(srcvec2, hivec);
    float32x4_t vecloadd2 = vmulq_f32(srcvec2, lovec);
    float32x4_t vechiadd3 = vmulq_f32(srcvec3, hivec);
    float32x4_t vecloadd3 = vmulq_f32(srcvec3, lovec);
    float32x4_t vechiadd4 = vmulq_f32(srcvec4, hivec);
    float32x4_t vecloadd4 = vmulq_f32(srcvec4, lovec);

    float32x2_t vechipair1 = vadd_f32(vget_high_f32(vechiadd1),
                                      vget_low_f32(vechiadd1));
    float32x2_t veclopair1 = vadd_f32(vget_high_f32(vecloadd1),
                                      vget_low_f32(vecloadd1));
    float32x2_t vechipair2 = vadd_f32(vget_high_f32(vechiadd2),
                                      vget_low_f32(vechiadd2));
    float32x2_t veclopair2 = vadd_f32(vget_high_f32(vecloadd2),
                                      vget_low_f32(vecloadd2));
    float32x2_t vechipair3 = vadd_f32(vget_high_f32(vechiadd3),
                                      vget_low_f32(vechiadd3));
    float32x2_t veclopair3 = vadd_f32(vget_high_f32(vecloadd3),
                                      vget_low_f32(vecloadd3));
    float32x2_t vechipair4 = vadd_f32(vget_high_f32(vechiadd4),
                                      vget_low_f32(vechiadd4));
    float32x2_t veclopair4 = vadd_f32(vget_high_f32(vecloadd4),
                                      vget_low_f32(vecloadd4));

    float32x2_t vecres1 = vpadd_f32(vechipair1, veclopair1);
    float32x2_t vecres2 = vpadd_f32(vechipair2, veclopair2);
    float32x2_t vecres3 = vpadd_f32(vechipair3, veclopair3);
    float32x2_t vecres4 = vpadd_f32(vechipair4, veclopair4);

    desthi[i] = vget_lane_f32(vecres1, 0);
    destlo[i] = vget_lane_f32(vecres1, 1);
    desthi[i + 1] = vget_lane_f32(vecres2, 0);
    destlo[i + 1] = vget_lane_f32(vecres2, 1);
    desthi[i + 2] = vget_lane_f32(vecres3, 0);
    destlo[i + 2] = vget_lane_f32(vecres3, 1);
    desthi[i + 3] = vget_lane_f32(vecres4, 0);
    destlo[i + 3] = vget_lane_f32(vecres4, 1);
  }
#else
#error This SIMD variant is not supported.
#endif  // #elif defined(__ARM_NEON__)
  // Finish with the extended end
  for (int i = simd_end; i < ilength; i++) {
    float reshi = 0.f, reslo = 0.f;
    for (int j = 0; j < 4; j++) {
      int index = i + j;
      float srcval = index < ilength? src[index] : src_ext[index - ilength];
      reshi += highpassC[j] * srcval;
      reslo += lowpassC[j] * srcval;
    }
    desthi[i] = reshi;
    destlo[i] = reslo;
  }
#else  // #ifdef SIMD
  stationary_wavelet_apply_na(type, 4 / stride, level, ext, src, length,
                              desthi, destlo);
#endif
}

static void wavelet_apply6(WaveletType type, ExtensionType ext,
                           const float *__restrict src, size_t length,
                           float *__restrict desthi,
                           float *__restrict destlo) {
#ifdef SIMD
  check_length(length);
  assert(src && desthi && destlo);

  if (align_complement_f32(src) != 0 || length < 8) {
    wavelet_apply_na(type, 6, ext, src, length, desthi, destlo);
    return;
  }

  int ilength = (int)length;
  DECLARE_PASSC(6);
  initialize_highpass_lowpass(type, 6, highpassC, lowpassC);
  float src_ext[6];
  initialize_extension(ext, 6, src, length, src_ext);

#ifdef __AVX__
  int simd_end = ilength - 6;
  highpassC[6] = 0.f;
  highpassC[7] = 0.f;
  lowpassC[6] = 0.f;
  lowpassC[7] = 0.f;
  const __m256 hpvec = _mm256_load_ps(highpassC);
  const __m256 lpvec = _mm256_load_ps(lowpassC);
  size_t alength = aligned_length(length);
  for (int i = 0, di = 0; i < simd_end; i += 2, di++) {
    int ex = (i / 2) % 4;
    int offset = i + (ex > 0? ex * (alength - 10) + 8 : 0);
    __m256 srcvec = _mm256_load_ps(src + offset);
    __m256 vechi = _mm256_dp_ps(srcvec, hpvec, 0xFF);
    __m256 veclo = _mm256_dp_ps(srcvec, lpvec, 0xFF);
    float reshi = vechi[0] + vechi[4];
    float reslo = veclo[0] + veclo[4];
    desthi[di] = reshi;
    destlo[di] = reslo;
  }
#elif defined(__ARM_NEON__)
  int simd_end = ilength - 4;
  const float32x4_t hivec1 = vld1q_f32(highpassC);
  const float32x2_t hivec2 = vld1_f32(&highpassC[4]);
  const float32x4_t lovec1 = vld1q_f32(lowpassC);
  const float32x2_t lovec2 = vld1_f32(&lowpassC[4]);
  for (int i = 0, di = 0; i < simd_end; i += 2, di++) {
    float32x4_t srcvec1 = vld1q_f32(src + i);
    float32x2_t srcvec2 = vld1_f32(src + i + 4);

    float32x4_t vechiadd1 = vmulq_f32(srcvec1, hivec1);
    float32x4_t vecloadd1 = vmulq_f32(srcvec1, lovec1);
    float32x2_t vechiadd2 = vmul_f32(srcvec2, hivec2);
    float32x2_t vecloadd2 = vmul_f32(srcvec2, lovec2);

    float32x2_t vechipair = vadd_f32(vget_high_f32(vechiadd1),
                                     vget_low_f32(vechiadd1));
    float32x2_t veclopair = vadd_f32(vget_high_f32(vecloadd1),
                                     vget_low_f32(vecloadd1));

    vechipair = vadd_f32(vechipair, vechiadd2);
    veclopair = vadd_f32(veclopair, vecloadd2);
    float32x2_t vecres = vpadd_f32(vechipair, veclopair);

    desthi[di] = vget_lane_f32(vecres, 0);
    destlo[di] = vget_lane_f32(vecres, 1);
  }
#else
#error This SIMD variant is not supported.
#endif  // #elif defined(__ARM_NEON__)
  // Finish with the extended end
  for (int i = simd_end, di = simd_end / 2; i < ilength; i += 2, di++) {
    float reshi = 0.f, reslo = 0.f;
    for (int j = 0; j < 6; j++) {
      int index = i + j;
      float srcval = index < ilength? src[index] : src_ext[index - ilength];
      reshi += highpassC[j] * srcval;
      reslo += lowpassC[j] * srcval;
    }
    desthi[di] = reshi;
    destlo[di] = reslo;
  }
#ifdef __AVX__
  if (length >= 16) {
    wavelet_prepare_array_memcpy(6, desthi, length / 2, desthi);
    wavelet_prepare_array_memcpy(6, destlo, length / 2, destlo);
  }
#endif  // #ifdef __AVX__
#else  // #ifdef SIMD
  wavelet_apply_na(type, 6, ext, src, length, desthi, destlo);
#endif
}

static void stationary_wavelet_apply6(WaveletType type, int level,
                                      ExtensionType ext,
                                      const float *__restrict src,
                                      size_t length,
                                      float *__restrict desthi,
                                      float *__restrict destlo) {
  int stride = 1 << (level - 1);
#ifdef SIMD
  assert(length > 0);
  assert(src && desthi && destlo);

  if (stride >= 4 || length < 8) {
    stationary_wavelet_apply_na(type, 6 / stride, level, ext, src, length,
                                desthi, destlo);
    return;
  }

  int ilength = (int)length;
  DECLARE_PASSC(6);
  stationary_initialize_highpass_lowpass(type, 6, level, highpassC, lowpassC);
  float src_ext[6];
  initialize_extension(ext, 6, src, length, src_ext);

#ifdef __AVX__
  int simd_end = ilength - 8;
  highpassC[6] = 0.f;
  highpassC[7] = 0.f;
  lowpassC[6] = 0.f;
  lowpassC[7] = 0.f;
  const __m256 hpvec = _mm256_load_ps(highpassC);
  const __m256 lpvec = _mm256_load_ps(lowpassC);
  for (int i = 0; i < simd_end; i += 2) {
    __m256 srcvec1 = _mm256_loadu_ps(src + i);
    __m256 srcvec2 = _mm256_loadu_ps(src + i + 1);
    __m256 vechi1 = _mm256_dp_ps(srcvec1, hpvec, 0xFF);
    __m256 veclo1 = _mm256_dp_ps(srcvec1, lpvec, 0xFF);
    __m256 vechi2 = _mm256_dp_ps(srcvec2, hpvec, 0xFF);
    __m256 veclo2 = _mm256_dp_ps(srcvec2, lpvec, 0xFF);
    float reshi1 = vechi1[0] + vechi1[4];
    float reslo1 = veclo1[0] + veclo1[4];
    float reshi2 = vechi2[0] + vechi2[4];
    float reslo2 = veclo2[0] + veclo2[4];
    desthi[i] = reshi1;
    destlo[i] = reslo1;
    desthi[i + 1] = reshi2;
    destlo[i + 1] = reslo2;
  }
#elif defined(__ARM_NEON__)
  int simd_end = ilength - 8;
  const float32x4_t hivec1 = vld1q_f32(highpassC);
  const float32x2_t hivec2 = vld1_f32(&highpassC[4]);
  const float32x4_t lovec1 = vld1q_f32(lowpassC);
  const float32x2_t lovec2 = vld1_f32(&lowpassC[4]);
  for (int i = 0; i < simd_end; i += 2) {
    float32x4_t srcvec1_1 = vld1q_f32(src + i);
    float32x2_t srcvec1_2 = vld1_f32(src + i + 4);
    float32x4_t srcvec2_1 = vld1q_f32(src + i + 1);
    float32x2_t srcvec2_2 = vld1_f32(src + i + 5);

    float32x4_t vechiadd1_1 = vmulq_f32(srcvec1_1, hivec1);
    float32x4_t vecloadd1_1 = vmulq_f32(srcvec1_1, lovec1);
    float32x2_t vechiadd1_2 = vmul_f32(srcvec1_2, hivec2);
    float32x2_t vecloadd1_2 = vmul_f32(srcvec1_2, lovec2);
    float32x4_t vechiadd2_1 = vmulq_f32(srcvec2_1, hivec1);
    float32x4_t vecloadd2_1 = vmulq_f32(srcvec2_1, lovec1);
    float32x2_t vechiadd2_2 = vmul_f32(srcvec2_2, hivec2);
    float32x2_t vecloadd2_2 = vmul_f32(srcvec2_2, lovec2);

    float32x2_t vechipair1 = vadd_f32(vget_high_f32(vechiadd1_1),
                                      vget_low_f32(vechiadd1_1));
    float32x2_t veclopair1 = vadd_f32(vget_high_f32(vecloadd1_1),
                                      vget_low_f32(vecloadd1_1));
    float32x2_t vechipair2 = vadd_f32(vget_high_f32(vechiadd2_1),
                                      vget_low_f32(vechiadd2_1));
    float32x2_t veclopair2 = vadd_f32(vget_high_f32(vecloadd2_1),
                                      vget_low_f32(vecloadd2_1));

    vechipair1 = vadd_f32(vechipair1, vechiadd1_2);
    veclopair1 = vadd_f32(veclopair1, vecloadd1_2);
    vechipair2 = vadd_f32(vechipair2, vechiadd2_2);
    veclopair2 = vadd_f32(veclopair2, vecloadd2_2);
    float32x2_t vecres1 = vpadd_f32(vechipair1, veclopair1);
    float32x2_t vecres2 = vpadd_f32(vechipair2, veclopair2);

    desthi[i] = vget_lane_f32(vecres1, 0);
    destlo[i] = vget_lane_f32(vecres1, 1);
    desthi[i + 1] = vget_lane_f32(vecres2, 0);
    destlo[i + 1] = vget_lane_f32(vecres2, 1);
  }
#else
#error This SIMD variant is not supported.
#endif  // #elif defined(__ARM_NEON__)
  // Finish with the extended end
  for (int i = simd_end;  i < ilength; i++) {
    float reshi = 0.f, reslo = 0.f;
    for (int j = 0; j < 6; j++) {
      int index = i + j;
      float srcval = index < ilength? src[index] : src_ext[index - ilength];
      reshi += highpassC[j] * srcval;
      reslo += lowpassC[j] * srcval;
    }
    desthi[i] = reshi;
    destlo[i] = reslo;
  }
#else  // #ifdef SIMD
  stationary_wavelet_apply_na(type, 6 / stride, level, ext, src, length,
                              desthi, destlo);
#endif
}

static void wavelet_apply8(WaveletType type, ExtensionType ext,
                           const float *__restrict src, size_t length,
                           float *__restrict desthi,
                           float *__restrict destlo) {
#ifdef SIMD
  check_length(length);
  assert(src && desthi && destlo);
  if (align_complement_f32(src) != 0 || length < 8) {
    wavelet_apply_na(type, 8, ext, src, length, desthi, destlo);
    return;
  }

  int ilength = (int)length;
  DECLARE_PASSC(8);
  initialize_highpass_lowpass(type, 8, highpassC, lowpassC);
  float src_ext[8];
  initialize_extension(ext, 8, src, length, src_ext);

#ifdef __AVX__
  const __m256 hpvec = _mm256_load_ps(highpassC);
  const __m256 lpvec = _mm256_load_ps(lowpassC);
  size_t alength = aligned_length(length);
  for (int i = 0, di = 0; i < ilength - 6; i += 2, di++) {
    int ex = (i / 2) % 4;
    int offset = i + (ex > 0? ex * (alength - 10) + 8 : 0);
    __m256 srcvec = _mm256_load_ps(src + offset);
    __m256 vechi = _mm256_dp_ps(srcvec, hpvec, 0xFF);
    __m256 veclo = _mm256_dp_ps(srcvec, lpvec, 0xFF);
    float reshi = vechi[0] + vechi[4];
    float reslo = veclo[0] + veclo[4];
    desthi[di] = reshi;
    destlo[di] = reslo;
  }
#elif defined(__ARM_NEON__)
  const float32x4_t hivec1 = vld1q_f32(highpassC);
  const float32x4_t hivec2 = vld1q_f32(&highpassC[4]);
  const float32x4_t lovec1 = vld1q_f32(lowpassC);
  const float32x4_t lovec2 = vld1q_f32(&lowpassC[4]);
  for (int i = 0, di = 0; i < ilength - 6; i += 2, di++) {
    float32x4_t srcvec1 = vld1q_f32(src + i);
    float32x4_t srcvec2 = vld1q_f32(src + i + 4);

    float32x4_t vechiadd = vmulq_f32(srcvec1, hivec1);
    float32x4_t vecloadd = vmulq_f32(srcvec1, lovec1);
    vechiadd = vmlaq_f32(vechiadd, srcvec2, hivec2);
    vecloadd = vmlaq_f32(vecloadd, srcvec2, lovec2);

    float32x2_t vechipair = vadd_f32(vget_high_f32(vechiadd),
                                     vget_low_f32(vechiadd));
    float32x2_t veclopair = vadd_f32(vget_high_f32(vecloadd),
                                     vget_low_f32(vecloadd));

    float32x2_t vecres = vpadd_f32(vechipair, veclopair);

    desthi[di] = vget_lane_f32(vecres, 0);
    destlo[di] = vget_lane_f32(vecres, 1);
  }
#else
#error This SIMD variant is not supported.
#endif  // #elif defined(__ARM_NEON__)
  // Finish with the extended end
  for (int i = ilength - 6, di = (ilength - 6) / 2;
       i < ilength; i += 2, di++) {
    float reshi = 0.f, reslo = 0.f;
    for (int j = 0; j < 8; j++) {
      int index = i + j;
      float srcval = index < ilength? src[index] : src_ext[index - ilength];
      reshi += highpassC[j] * srcval;
      reslo += lowpassC[j] * srcval;
    }
    desthi[di] = reshi;
    destlo[di] = reslo;
  }
#ifdef __AVX__
  if (length >= 16) {
    wavelet_prepare_array_memcpy(8, desthi, length / 2, desthi);
    wavelet_prepare_array_memcpy(8, destlo, length / 2, destlo);
  }
#endif  // #ifdef __AVX__
#else  // #ifdef SIMD
  wavelet_apply_na(type, 8, ext, src, length, desthi, destlo);
#endif
}

static void stationary_wavelet_apply8(WaveletType type, int level,
                                      ExtensionType ext,
                                      const float *__restrict src,
                                      size_t length,
                                      float *__restrict desthi,
                                      float *__restrict destlo) {
  int stride = 1 << (level - 1);
#ifdef SIMD
  assert(length > 0);
  assert(src && desthi && destlo);

  if (stride >= 4 || length < 8) {
    stationary_wavelet_apply_na(type, 8 / stride, level, ext, src, length,
                                desthi, destlo);
    return;
  }

  int ilength = (int)length;
  DECLARE_PASSC(8);
  stationary_initialize_highpass_lowpass(type, 8, level, highpassC, lowpassC);
  float src_ext[8];
  initialize_extension(ext, 8, src, length, src_ext);

#ifdef __AVX__
  int simd_end = ilength - 8;
  const __m256 hpvec = _mm256_load_ps(highpassC);
  const __m256 lpvec = _mm256_load_ps(lowpassC);
  for (int i = 0; i < simd_end; i += 2) {
    __m256 srcvec1 = _mm256_loadu_ps(src + i);
    __m256 srcvec2 = _mm256_loadu_ps(src + i + 1);
    __m256 vechi1 = _mm256_dp_ps(srcvec1, hpvec, 0xFF);
    __m256 veclo1 = _mm256_dp_ps(srcvec1, lpvec, 0xFF);
    __m256 vechi2 = _mm256_dp_ps(srcvec2, hpvec, 0xFF);
    __m256 veclo2 = _mm256_dp_ps(srcvec2, lpvec, 0xFF);
    float reshi1 = vechi1[0] + vechi1[4];
    float reslo1 = veclo1[0] + veclo1[4];
    float reshi2 = vechi2[0] + vechi2[4];
    float reslo2 = veclo2[0] + veclo2[4];
    desthi[i] = reshi1;
    destlo[i] = reslo1;
    desthi[i + 1] = reshi2;
    destlo[i + 1] = reslo2;
  }
#elif defined(__ARM_NEON__)
  int simd_end = ilength - 7;
  const float32x4_t hivec1 = vld1q_f32(highpassC);
  const float32x4_t hivec2 = vld1q_f32(&highpassC[4]);
  const float32x4_t lovec1 = vld1q_f32(lowpassC);
  const float32x4_t lovec2 = vld1q_f32(&lowpassC[4]);
  for (int i = 0; i < simd_end; i++) {
    float32x4_t srcvec1 = vld1q_f32(src + i);
    float32x4_t srcvec2 = vld1q_f32(src + i + 4);

    float32x4_t vechiadd = vmulq_f32(srcvec1, hivec1);
    float32x4_t vecloadd = vmulq_f32(srcvec1, lovec1);
    vechiadd = vmlaq_f32(vechiadd, srcvec2, hivec2);
    vecloadd = vmlaq_f32(vecloadd, srcvec2, lovec2);

    float32x2_t vechipair = vadd_f32(vget_high_f32(vechiadd),
                                     vget_low_f32(vechiadd));
    float32x2_t veclopair = vadd_f32(vget_high_f32(vecloadd),
                                     vget_low_f32(vecloadd));

    float32x2_t vecres = vpadd_f32(vechipair, veclopair);

    desthi[i] = vget_lane_f32(vecres, 0);
    destlo[i] = vget_lane_f32(vecres, 1);
  }
#else
#error This SIMD variant is not supported.
#endif  // #elif defined(__ARM_NEON__)
  // Finish with the extended end
  for (int i = simd_end; i < ilength; i++) {
    float reshi = 0.f, reslo = 0.f;
    for (int j = 0; j < 8; j++) {
      int index = i + j;
      float srcval = index < ilength? src[index] : src_ext[index - ilength];
      reshi += highpassC[j] * srcval;
      reslo += lowpassC[j] * srcval;
    }
    desthi[i] = reshi;
    destlo[i] = reslo;
  }
#else  // #ifdef SIMD
  stationary_wavelet_apply_na(type, 8 / stride, level, ext, src, length,
                              desthi, destlo);
#endif
}

static void wavelet_apply12(WaveletType type, ExtensionType ext,
                            const float *__restrict src, size_t length,
                            float *__restrict desthi,
                            float *__restrict destlo) {
#ifdef SIMD
  check_length(length);
  assert(src && desthi && destlo);

  if (align_complement_f32(src) != 0 ||
#ifdef __AVX__
      length < 16
#elif defined(__ARM_NEON__)
      length < 12
#endif
  ) {
    wavelet_apply_na(type, 12, ext, src, length, desthi, destlo);
    return;
  }

  int ilength = (int)length;
  DECLARE_PASSC(12);
  initialize_highpass_lowpass(type, 12, highpassC, lowpassC);
  float src_ext[12];
  initialize_extension(ext, 12, src, length, src_ext);

#ifdef __AVX__
  int simd_end = ilength - 14;
  highpassC[12] = 0.f;
  highpassC[13] = 0.f;
  highpassC[14] = 0.f;
  highpassC[15] = 0.f;
  lowpassC[12] = 0.f;
  lowpassC[13] = 0.f;
  lowpassC[14] = 0.f;
  lowpassC[15] = 0.f;
  const __m256 hpvec1 = _mm256_load_ps(highpassC);
  const __m256 lpvec1 = _mm256_load_ps(lowpassC);
  const __m256 hpvec2 = _mm256_load_ps(&highpassC[8]);
  const __m256 lpvec2 = _mm256_load_ps(&lowpassC[8]);
  size_t alength = aligned_length(length);
  for (int i = 0, di = 0; i < simd_end; i += 2, di++) {
    int ex = (i / 2) % 4;
    int offset = i + (ex > 0? ex * (alength - 10) + 8 : 0);
    __m256 srcvec1 = _mm256_load_ps(src + offset);
    __m256 srcvec2 = _mm256_load_ps(src + offset + 8);
    __m256 vechi1 = _mm256_dp_ps(srcvec1, hpvec1, 0xFF);
    __m256 veclo1 = _mm256_dp_ps(srcvec1, lpvec1, 0xFF);
    __m256 vechi2 = _mm256_dp_ps(srcvec2, hpvec2, 0xFF);
    __m256 veclo2 = _mm256_dp_ps(srcvec2, lpvec2, 0xFF);
    __m256 vechi = _mm256_add_ps(vechi1, vechi2);
    __m256 veclo = _mm256_add_ps(veclo1, veclo2);
    float reshi = vechi[0] + vechi[4];
    float reslo = veclo[0] + veclo[4];
    desthi[di] = reshi;
    destlo[di] = reslo;
  }
#elif defined(__ARM_NEON__)
  int simd_end = ilength - 10;
  const float32x4_t hivec1 = vld1q_f32(highpassC);
  const float32x4_t hivec2 = vld1q_f32(&highpassC[4]);
  const float32x4_t hivec3 = vld1q_f32(&highpassC[8]);
  const float32x4_t lovec1 = vld1q_f32(lowpassC);
  const float32x4_t lovec2 = vld1q_f32(&lowpassC[4]);
  const float32x4_t lovec3 = vld1q_f32(&lowpassC[8]);
  for (int i = 0, di = 0; i < simd_end; i += 2, di++) {
    float32x4_t srcvec1 = vld1q_f32(src + i);
    float32x4_t srcvec2 = vld1q_f32(src + i + 4);
    float32x4_t srcvec3 = vld1q_f32(src + i + 8);

    float32x4_t vechiadd = vmulq_f32(srcvec1, hivec1);
    float32x4_t vecloadd = vmulq_f32(srcvec1, lovec1);
    vechiadd = vmlaq_f32(vechiadd, srcvec2, hivec2);
    vecloadd = vmlaq_f32(vecloadd, srcvec2, lovec2);
    vechiadd = vmlaq_f32(vechiadd, srcvec3, hivec3);
    vecloadd = vmlaq_f32(vecloadd, srcvec3, lovec3);

    float32x2_t vechipair = vadd_f32(vget_high_f32(vechiadd),
                                     vget_low_f32(vechiadd));
    float32x2_t veclopair = vadd_f32(vget_high_f32(vecloadd),
                                     vget_low_f32(vecloadd));

    float32x2_t vecres = vpadd_f32(vechipair, veclopair);

    desthi[di] = vget_lane_f32(vecres, 0);
    destlo[di] = vget_lane_f32(vecres, 1);
  }
#else
#error This SIMD variant is not supported.
#endif  // #elif defined(__ARM_NEON__)
  // Finish with the extended end
  for (int i = simd_end, di = simd_end / 2; i < ilength; i += 2, di++) {
    float reshi = 0.f, reslo = 0.f;
    for (int j = 0; j < 12; j++) {
      int index = i + j;
      float srcval = index < ilength? src[index] : src_ext[index - ilength];
      reshi += highpassC[j] * srcval;
      reslo += lowpassC[j] * srcval;
    }
    desthi[di] = reshi;
    destlo[di] = reslo;
  }
#ifdef __AVX__
  if (length >= 24) {
    wavelet_prepare_array_memcpy(12, desthi, length / 2, desthi);
    wavelet_prepare_array_memcpy(12, destlo, length / 2, destlo);
  }
#endif  // #ifdef __AVX__
#else  // #ifdef SIMD
  wavelet_apply_na(type, 12, ext, src, length, desthi, destlo);
#endif
}

static void stationary_wavelet_apply12(WaveletType type, int level,
                                       ExtensionType ext,
                                       const float *__restrict src,
                                       size_t length,
                                       float *__restrict desthi,
                                       float *__restrict destlo) {
  int stride = 1 << (level - 1);
#ifdef SIMD
  assert(length > 0);
  assert(src && desthi && destlo);

  if (stride >= 4 ||
#ifdef __AVX__
      length < 16
#elif defined(__ARM_NEON__)
      length < 12
#else
      0
#endif
  ) {
    stationary_wavelet_apply_na(type, 12 / stride, level, ext, src, length,
                                desthi, destlo);
    return;
  }

  int ilength = (int)length;
  DECLARE_PASSC(12);
  stationary_initialize_highpass_lowpass(type, 12, level, highpassC, lowpassC);
  float src_ext[12];
  initialize_extension(ext, 12, src, length, src_ext);

#ifdef __AVX__
  int simd_end = ilength - 16;
  highpassC[12] = 0.f;
  highpassC[13] = 0.f;
  highpassC[14] = 0.f;
  highpassC[15] = 0.f;
  lowpassC[12] = 0.f;
  lowpassC[13] = 0.f;
  lowpassC[14] = 0.f;
  lowpassC[15] = 0.f;
  const __m256 hpvec1 = _mm256_load_ps(highpassC);
  const __m256 lpvec1 = _mm256_load_ps(lowpassC);
  const __m256 hpvec2 = _mm256_load_ps(&highpassC[8]);
  const __m256 lpvec2 = _mm256_load_ps(&lowpassC[8]);
  for (int i = 0; i < simd_end; i += 2) {
    __m256 srcvec1_1 = _mm256_loadu_ps(src + i);
    __m256 srcvec1_2 = _mm256_loadu_ps(src + i + 8);
    __m256 srcvec2_1 = _mm256_loadu_ps(src + i + 1);
    __m256 srcvec2_2 = _mm256_loadu_ps(src + i + 9);

    __m256 vechi1_1 = _mm256_dp_ps(srcvec1_1, hpvec1, 0xFF);
    __m256 veclo1_1 = _mm256_dp_ps(srcvec1_1, lpvec1, 0xFF);
    __m256 vechi1_2 = _mm256_dp_ps(srcvec1_2, hpvec2, 0xFF);
    __m256 veclo1_2 = _mm256_dp_ps(srcvec1_2, lpvec2, 0xFF);
    __m256 vechi2_1 = _mm256_dp_ps(srcvec2_1, hpvec1, 0xFF);
    __m256 veclo2_1 = _mm256_dp_ps(srcvec2_1, lpvec1, 0xFF);
    __m256 vechi2_2 = _mm256_dp_ps(srcvec2_2, hpvec2, 0xFF);
    __m256 veclo2_2 = _mm256_dp_ps(srcvec2_2, lpvec2, 0xFF);

    __m256 vechi1 = _mm256_add_ps(vechi1_1, vechi1_2);
    __m256 veclo1 = _mm256_add_ps(veclo1_1, veclo1_2);
    __m256 vechi2 = _mm256_add_ps(vechi2_1, vechi2_2);
    __m256 veclo2 = _mm256_add_ps(veclo2_1, veclo2_2);
    float reshi1 = vechi1[0] + vechi1[4];
    float reslo1 = veclo1[0] + veclo1[4];
    float reshi2 = vechi2[0] + vechi2[4];
    float reslo2 = veclo2[0] + veclo2[4];
    desthi[i] = reshi1;
    destlo[i] = reslo1;
    desthi[i + 1] = reshi2;
    destlo[i + 1] = reslo2;
  }
#elif defined(__ARM_NEON__)
  int simd_end = ilength - 11;
  const float32x4_t hivec1 = vld1q_f32(highpassC);
  const float32x4_t hivec2 = vld1q_f32(&highpassC[4]);
  const float32x4_t hivec3 = vld1q_f32(&highpassC[8]);
  const float32x4_t lovec1 = vld1q_f32(lowpassC);
  const float32x4_t lovec2 = vld1q_f32(&lowpassC[4]);
  const float32x4_t lovec3 = vld1q_f32(&lowpassC[8]);
  for (int i = 0; i < simd_end; i++) {
    float32x4_t srcvec1 = vld1q_f32(src + i);
    float32x4_t srcvec2 = vld1q_f32(src + i + 4);
    float32x4_t srcvec3 = vld1q_f32(src + i + 8);

    float32x4_t vechiadd = vmulq_f32(srcvec1, hivec1);
    float32x4_t vecloadd = vmulq_f32(srcvec1, lovec1);
    vechiadd = vmlaq_f32(vechiadd, srcvec2, hivec2);
    vecloadd = vmlaq_f32(vecloadd, srcvec2, lovec2);
    vechiadd = vmlaq_f32(vechiadd, srcvec3, hivec3);
    vecloadd = vmlaq_f32(vecloadd, srcvec3, lovec3);

    float32x2_t vechipair = vadd_f32(vget_high_f32(vechiadd),
                                     vget_low_f32(vechiadd));
    float32x2_t veclopair = vadd_f32(vget_high_f32(vecloadd),
                                     vget_low_f32(vecloadd));

    float32x2_t vecres = vpadd_f32(vechipair, veclopair);

    desthi[i] = vget_lane_f32(vecres, 0);
    destlo[i] = vget_lane_f32(vecres, 1);
  }
#else
#error This SIMD variant is not supported.
#endif  // #elif defined(__ARM_NEON__)
  // Finish with the extended end
  for (int i = simd_end; i < ilength; i++) {
    float reshi = 0.f, reslo = 0.f;
    for (int j = 0; j < 12; j++) {
      int index = i + j;
      float srcval = index < ilength? src[index] : src_ext[index - ilength];
      reshi += highpassC[j] * srcval;
      reslo += lowpassC[j] * srcval;
    }
    desthi[i] = reshi;
    destlo[i] = reslo;
  }
#else  // #ifdef SIMD
  stationary_wavelet_apply_na(type, 12 / stride, level, ext, src, length,
                              desthi, destlo);
#endif
}

static void wavelet_apply16(WaveletType type, ExtensionType ext,
                            const float *__restrict src, size_t length,
                            float *__restrict desthi,
                            float *__restrict destlo) {
#ifdef SIMD
  check_length(length);
  assert(src && desthi && destlo);

  if (align_complement_f32(src) != 0 || length < 16) {
    wavelet_apply_na(type, 16, ext, src, length, desthi, destlo);
    return;
  }

  int ilength = (int)length;
  DECLARE_PASSC(16);
  initialize_highpass_lowpass(type, 16, highpassC, lowpassC);
  float src_ext[16];
  initialize_extension(ext, 16, src, length, src_ext);

#ifdef __AVX__
  const __m256 hpvec1 = _mm256_load_ps(highpassC);
  const __m256 lpvec1 = _mm256_load_ps(lowpassC);
  const __m256 hpvec2 = _mm256_load_ps(&highpassC[8]);
  const __m256 lpvec2 = _mm256_load_ps(&lowpassC[8]);
  size_t alength = aligned_length(length);
  for (int i = 0, di = 0; i < ilength - 14; i += 2, di++) {
    int ex = (i / 2) % 4;
    int offset = i + (ex > 0? ex * (alength - 10) + 8 : 0);
    __m256 srcvec1 = _mm256_load_ps(src + offset);
    __m256 srcvec2 = _mm256_load_ps(src + offset + 8);
    __m256 vechi1 = _mm256_dp_ps(srcvec1, hpvec1, 0xFF);
    __m256 veclo1 = _mm256_dp_ps(srcvec1, lpvec1, 0xFF);
    __m256 vechi2 = _mm256_dp_ps(srcvec2, hpvec2, 0xFF);
    __m256 veclo2 = _mm256_dp_ps(srcvec2, lpvec2, 0xFF);
    __m256 vechi = _mm256_add_ps(vechi1, vechi2);
    __m256 veclo = _mm256_add_ps(veclo1, veclo2);
    float reshi = vechi[0] + vechi[4];
    float reslo = veclo[0] + veclo[4];
    desthi[di] = reshi;
    destlo[di] = reslo;
  }
#elif defined(__ARM_NEON__)
  const float32x4_t hivec1 = vld1q_f32(highpassC);
  const float32x4_t hivec2 = vld1q_f32(&highpassC[4]);
  const float32x4_t hivec3 = vld1q_f32(&highpassC[8]);
  const float32x4_t hivec4 = vld1q_f32(&highpassC[12]);
  const float32x4_t lovec1 = vld1q_f32(lowpassC);
  const float32x4_t lovec2 = vld1q_f32(&lowpassC[4]);
  const float32x4_t lovec3 = vld1q_f32(&lowpassC[8]);
  const float32x4_t lovec4 = vld1q_f32(&lowpassC[12]);
  for (int i = 0, di = 0; i < ilength - 15; i += 2, di++) {
    float32x4_t srcvec1 = vld1q_f32(src + i);
    float32x4_t srcvec2 = vld1q_f32(src + i + 4);
    float32x4_t srcvec3 = vld1q_f32(src + i + 8);
    float32x4_t srcvec4 = vld1q_f32(src + i + 12);

    float32x4_t vechiadd = vmulq_f32(srcvec1, hivec1);
    float32x4_t vecloadd = vmulq_f32(srcvec1, lovec1);
    vechiadd = vmlaq_f32(vechiadd, srcvec2, hivec2);
    vecloadd = vmlaq_f32(vecloadd, srcvec2, lovec2);
    vechiadd = vmlaq_f32(vechiadd, srcvec3, hivec3);
    vecloadd = vmlaq_f32(vecloadd, srcvec3, lovec3);
    vechiadd = vmlaq_f32(vechiadd, srcvec4, hivec4);
    vecloadd = vmlaq_f32(vecloadd, srcvec4, lovec4);

    float32x2_t vechipair = vadd_f32(vget_high_f32(vechiadd),
                                     vget_low_f32(vechiadd));
    float32x2_t veclopair = vadd_f32(vget_high_f32(vecloadd),
                                     vget_low_f32(vecloadd));

    float32x2_t vecres = vpadd_f32(vechipair, veclopair);

    float reshi = vget_lane_f32(vecres, 0);
    float reslo = vget_lane_f32(vecres, 1);

    desthi[di] = reshi;
    destlo[di] = reslo;
  }
#else
#error This SIMD variant is not supported.
#endif  // #elif defined(__ARM_NEON__)
  // Finish with the extended end
  for (int i = ilength - 14, di = (ilength - 14) / 2;
       i < ilength; i += 2, di++) {
    float reshi = 0.f, reslo = 0.f;
    for (int j = 0; j < 16; j++) {
      int index = i + j;
      float srcval = index < ilength? src[index] : src_ext[index - ilength];
      reshi += highpassC[j] * srcval;
      reslo += lowpassC[j] * srcval;
    }
    desthi[di] = reshi;
    destlo[di] = reslo;
  }
#ifdef __AVX__
  if (length >= 32) {
    wavelet_prepare_array_memcpy(16, desthi, length / 2, desthi);
    wavelet_prepare_array_memcpy(16, destlo, length / 2, destlo);
  }
#endif  // #ifdef __AVX__
#else  // #ifdef SIMD
  wavelet_apply_na(type, 16, ext, src, length, desthi, destlo);
#endif
}

static void stationary_wavelet_apply16(WaveletType type, int level,
                                       ExtensionType ext,
                                       const float *__restrict src,
                                       size_t length,
                                       float *__restrict desthi,
                                       float *__restrict destlo) {
  int stride = 1 << (level - 1);
#ifdef SIMD
  assert(length > 0);
  assert(src && desthi && destlo);

  if (stride >= 4 || length < 16) {
    stationary_wavelet_apply_na(type, 16 / stride, level, ext, src, length,
                                desthi, destlo);
    return;
  }

  int ilength = (int)length;
  DECLARE_PASSC(16);
  stationary_initialize_highpass_lowpass(type, 16, level, highpassC, lowpassC);
  float src_ext[16];
  initialize_extension(ext, 16, src, length, src_ext);

#ifdef __AVX__
  int simd_end = ilength - 16;
  const __m256 hpvec1 = _mm256_load_ps(highpassC);
  const __m256 lpvec1 = _mm256_load_ps(lowpassC);
  const __m256 hpvec2 = _mm256_load_ps(&highpassC[8]);
  const __m256 lpvec2 = _mm256_load_ps(&lowpassC[8]);
  for (int i = 0; i < simd_end; i += 2) {
    __m256 srcvec1_1 = _mm256_loadu_ps(src + i);
    __m256 srcvec1_2 = _mm256_loadu_ps(src + i + 8);
    __m256 srcvec2_1 = _mm256_loadu_ps(src + i + 1);
    __m256 srcvec2_2 = _mm256_loadu_ps(src + i + 9);

    __m256 vechi1_1 = _mm256_dp_ps(srcvec1_1, hpvec1, 0xFF);
    __m256 veclo1_1 = _mm256_dp_ps(srcvec1_1, lpvec1, 0xFF);
    __m256 vechi1_2 = _mm256_dp_ps(srcvec1_2, hpvec2, 0xFF);
    __m256 veclo1_2 = _mm256_dp_ps(srcvec1_2, lpvec2, 0xFF);
    __m256 vechi2_1 = _mm256_dp_ps(srcvec2_1, hpvec1, 0xFF);
    __m256 veclo2_1 = _mm256_dp_ps(srcvec2_1, lpvec1, 0xFF);
    __m256 vechi2_2 = _mm256_dp_ps(srcvec2_2, hpvec2, 0xFF);
    __m256 veclo2_2 = _mm256_dp_ps(srcvec2_2, lpvec2, 0xFF);

    __m256 vechi1 = _mm256_add_ps(vechi1_1, vechi1_2);
    __m256 veclo1 = _mm256_add_ps(veclo1_1, veclo1_2);
    __m256 vechi2 = _mm256_add_ps(vechi2_1, vechi2_2);
    __m256 veclo2 = _mm256_add_ps(veclo2_1, veclo2_2);
    float reshi1 = vechi1[0] + vechi1[4];
    float reslo1 = veclo1[0] + veclo1[4];
    float reshi2 = vechi2[0] + vechi2[4];
    float reslo2 = veclo2[0] + veclo2[4];
    desthi[i] = reshi1;
    destlo[i] = reslo1;
    desthi[i + 1] = reshi2;
    destlo[i + 1] = reslo2;
  }
#elif defined(__ARM_NEON__)
  int simd_end = ilength - 15;
  const float32x4_t hivec1 = vld1q_f32(highpassC);
  const float32x4_t hivec2 = vld1q_f32(&highpassC[4]);
  const float32x4_t hivec3 = vld1q_f32(&highpassC[8]);
  const float32x4_t hivec4 = vld1q_f32(&highpassC[12]);
  const float32x4_t lovec1 = vld1q_f32(lowpassC);
  const float32x4_t lovec2 = vld1q_f32(&lowpassC[4]);
  const float32x4_t lovec3 = vld1q_f32(&lowpassC[8]);
  const float32x4_t lovec4 = vld1q_f32(&lowpassC[12]);
  for (int i = 0; i < simd_end; i++) {
    float32x4_t srcvec1 = vld1q_f32(src + i);
    float32x4_t srcvec2 = vld1q_f32(src + i + 4);
    float32x4_t srcvec3 = vld1q_f32(src + i + 8);
    float32x4_t srcvec4 = vld1q_f32(src + i + 12);

    float32x4_t vechiadd = vmulq_f32(srcvec1, hivec1);
    float32x4_t vecloadd = vmulq_f32(srcvec1, lovec1);
    vechiadd = vmlaq_f32(vechiadd, srcvec2, hivec2);
    vecloadd = vmlaq_f32(vecloadd, srcvec2, lovec2);
    vechiadd = vmlaq_f32(vechiadd, srcvec3, hivec3);
    vecloadd = vmlaq_f32(vecloadd, srcvec3, lovec3);
    vechiadd = vmlaq_f32(vechiadd, srcvec4, hivec4);
    vecloadd = vmlaq_f32(vecloadd, srcvec4, lovec4);

    float32x2_t vechipair = vadd_f32(vget_high_f32(vechiadd),
                                     vget_low_f32(vechiadd));
    float32x2_t veclopair = vadd_f32(vget_high_f32(vecloadd),
                                     vget_low_f32(vecloadd));

    float32x2_t vecres = vpadd_f32(vechipair, veclopair);

    float reshi = vget_lane_f32(vecres, 0);
    float reslo = vget_lane_f32(vecres, 1);

    desthi[i] = reshi;
    destlo[i] = reslo;
  }
#else
#error This SIMD variant is not supported.
#endif  // #elif defined(__ARM_NEON__)
  // Finish with the extended end
  for (int i = simd_end; i < ilength; i++) {
    float reshi = 0.f, reslo = 0.f;
    for (int j = 0; j < 16; j++) {
      int index = i + j;
      float srcval = index < ilength? src[index] : src_ext[index - ilength];
      reshi += highpassC[j] * srcval;
      reslo += lowpassC[j] * srcval;
    }
    desthi[i] = reshi;
    destlo[i] = reslo;
  }
#else  // #ifdef SIMD
  stationary_wavelet_apply_na(type, 16 / stride, level, ext, src, length,
                              desthi, destlo);
#endif
}

#ifdef SIMD
static void stationary_wavelet_applyN_core(const float* highpassC,
                                           const float* lowpassC, int size,
                                           const float *restrict src,
                                           size_t length,
                                           float *restrict desthi,
                                           float *restrict destlo) {
#ifdef __AVX__
  assert(size % 16 == 0);
  for (int i = 0; i < (int)length - size + 1; i++) {
    __m256 rvechi = _mm256_set1_ps(0), rveclo = _mm256_set1_ps(0);
    for (int j = 0; j < size; j += 16) {
      __m256 srcvec1 = _mm256_loadu_ps(src + i + j);
      __m256 srcvec2 = _mm256_loadu_ps(src + i + j + 8);
      __m256 hivec1 = _mm256_loadu_ps(highpassC + j);
      __m256 lovec1 = _mm256_loadu_ps(lowpassC + j);
      __m256 hivec2 = _mm256_loadu_ps(highpassC + j + 8);
      __m256 lovec2 = _mm256_loadu_ps(lowpassC + j + 8);
      __m256 mulvechi1 = _mm256_mul_ps(srcvec1, hivec1);
      __m256 mulveclo1 = _mm256_mul_ps(srcvec1, lovec1);
      __m256 mulvechi2 = _mm256_mul_ps(srcvec2, hivec2);
      __m256 mulveclo2 = _mm256_mul_ps(srcvec2, lovec2);

      rvechi = _mm256_add_ps(rvechi, mulvechi1);
      rvechi = _mm256_add_ps(rvechi, mulvechi2);
      rveclo = _mm256_add_ps(rveclo, mulveclo1);
      rveclo = _mm256_add_ps(rveclo, mulveclo2);
    }
    rvechi = _mm256_hadd_ps(rvechi, rvechi);
    rvechi = _mm256_hadd_ps(rvechi, rvechi);
    rveclo = _mm256_hadd_ps(rveclo, rveclo);
    rveclo = _mm256_hadd_ps(rveclo, rveclo);

    float reshi = rvechi[0] + rvechi[4];
    float reslo = rveclo[0] + rveclo[4];

    desthi[i] = reshi;
    destlo[i] = reslo;
  }
#elif defined(__ARM_NEON__)
  assert(size % 8 == 0);
  for (int i = 0; i < (int)length - size + 1; i++) {
    float32x4_t rvechi = vdupq_n_f32(0), rveclo = vdupq_n_f32(0);
    for (int j = 0; j < size; j += 8) {
      float32x4_t srcvec1 = vld1q_f32(src + i + j);
      float32x4_t srcvec2 = vld1q_f32(src + i + j + 4);
      float32x4_t hivec1 = vld1q_f32(highpassC + j);
      float32x4_t lovec1 = vld1q_f32(lowpassC + j);
      float32x4_t hivec2 = vld1q_f32(highpassC + j + 4);
      float32x4_t lovec2 = vld1q_f32(lowpassC + j + 4);
      rvechi = vmlaq_f32(rvechi, srcvec1, hivec1);
      rveclo = vmlaq_f32(rveclo, srcvec1, lovec1);
      rvechi = vmlaq_f32(rvechi, srcvec2, hivec2);
      rveclo = vmlaq_f32(rveclo, srcvec2, lovec2);
    }

    float32x2_t vechipair = vadd_f32(vget_high_f32(rvechi),
                                     vget_low_f32(rvechi));
    float32x2_t veclopair = vadd_f32(vget_high_f32(rveclo),
                                     vget_low_f32(rveclo));
    float32x2_t vecres = vpadd_f32(vechipair, veclopair);

    desthi[i] = vget_lane_f32(vecres, 0);
    destlo[i] = vget_lane_f32(vecres, 1);
  }
#else
#error This SIMD variant is not supported.
#endif  // #elif defined(__ARM_NEON__)
}
#endif  // SIMD

static void stationary_wavelet_apply24(WaveletType type, int level,
                                       ExtensionType ext,
                                       const float *__restrict src,
                                       size_t length,
                                       float *__restrict desthi,
                                       float *__restrict destlo) {
  int stride = 1 << (level - 1);
#ifdef SIMD
  assert(length > 0);
  assert(src && desthi && destlo);

  if (stride >= 4 ||
#ifdef __AVX__
      length < 32
#elif defined(__ARM_NEON__)
      length < 24
#else
      0
#endif
  ) {
    stationary_wavelet_apply_na(type, 24 / stride, level, ext, src, length,
                                desthi, destlo);
    return;
  }

  int ilength = (int)length;
  DECLARE_PASSC(24);
  stationary_initialize_highpass_lowpass(type, 24, level, highpassC, lowpassC);
  float src_ext[24];
  initialize_extension(ext, 24, src, length, src_ext);

  int simd_end = ilength - 23;
#ifdef __AVX__
  const __m256 hpvec1 = _mm256_load_ps(highpassC);
  const __m256 lpvec1 = _mm256_load_ps(lowpassC);
  const __m256 hpvec2 = _mm256_load_ps(&highpassC[8]);
  const __m256 lpvec2 = _mm256_load_ps(&lowpassC[8]);
  const __m256 hpvec3 = _mm256_load_ps(&highpassC[16]);
  const __m256 lpvec3 = _mm256_load_ps(&lowpassC[16]);
  for (int i = 0; i < simd_end; i++) {
    __m256 srcvec1 = _mm256_loadu_ps(src + i);
    __m256 srcvec2 = _mm256_loadu_ps(src + i + 8);
    __m256 srcvec3 = _mm256_loadu_ps(src + i + 16);

    __m256 vechi1 = _mm256_mul_ps(srcvec1, hpvec1);
    __m256 veclo1 = _mm256_mul_ps(srcvec1, lpvec1);
    __m256 vechi2 = _mm256_mul_ps(srcvec2, hpvec2);
    __m256 veclo2 = _mm256_mul_ps(srcvec2, lpvec2);
    __m256 vechi3 = _mm256_mul_ps(srcvec3, hpvec3);
    __m256 veclo3 = _mm256_mul_ps(srcvec3, lpvec3);

    __m256 vechi = _mm256_add_ps(vechi1, vechi2);
    vechi = _mm256_add_ps(vechi, vechi3);
    vechi = _mm256_hadd_ps(vechi, vechi);
    vechi = _mm256_hadd_ps(vechi, vechi);
    __m256 veclo = _mm256_add_ps(veclo1, veclo2);
    veclo = _mm256_add_ps(veclo, veclo3);
    veclo = _mm256_hadd_ps(veclo, veclo);
    veclo = _mm256_hadd_ps(veclo, veclo);

    float reshi = vechi[0] + vechi[4];
    float reslo = veclo[0] + veclo[4];

    desthi[i] = reshi;
    destlo[i] = reslo;
  }
#elif defined(__ARM_NEON__)
  stationary_wavelet_applyN_core(highpassC, lowpassC, 24, src, length,
                                 desthi, destlo);
#else
#error This SIMD variant is not supported.
#endif  // #elif defined(__ARM_NEON__)
  // Finish with the extended end
  for (int i = simd_end; i < ilength; i++) {
    float reshi = 0.f, reslo = 0.f;
    for (int j = 0; j < 24; j++) {
      int index = i + j;
      float srcval = index < ilength? src[index] : src_ext[index - ilength];
      reshi += highpassC[j] * srcval;
      reslo += lowpassC[j] * srcval;
    }
    desthi[i] = reshi;
    destlo[i] = reslo;
  }
#else  // #ifdef SIMD
  stationary_wavelet_apply_na(type, 24 / stride, level, ext, src, length,
                              desthi, destlo);
#endif
}

static void stationary_wavelet_applyN(WaveletType type, int size, int level,
                                      ExtensionType ext,
                                      const float *__restrict src,
                                      size_t length,
                                      float *__restrict desthi,
                                      float *__restrict destlo) {
  int stride = 1 << (level - 1);
#ifdef SIMD
  assert(length > 0);
  assert(src && desthi && destlo);

  if (stride >= 4 || length < (size_t)size) {
    stationary_wavelet_apply_na(type, size / stride, level, ext, src, length,
                                desthi, destlo);
    return;
  }

  int ilength = (int)length;
  DECLARE_PASSC(size);
  stationary_initialize_highpass_lowpass(type, size, level,
                                         highpassC, lowpassC);
  float src_ext[size];
  initialize_extension(ext, size, src, length, src_ext);

  stationary_wavelet_applyN_core(highpassC, lowpassC, size, src, length,
                                 desthi, destlo);
  // Finish with the extended end
  for (int i = ilength - size + 1; i < ilength; i++) {
    float reshi = 0.f, reslo = 0.f;
    for (int j = 0; j < size; j++) {
      int index = i + j;
      float srcval = index < ilength? src[index] : src_ext[index - ilength];
      reshi += highpassC[j] * srcval;
      reslo += lowpassC[j] * srcval;
    }
    desthi[i] = reshi;
    destlo[i] = reslo;
  }
#else  // #ifdef SIMD
  stationary_wavelet_apply_na(type, size / stride, level, ext, src, length,
                              desthi, destlo);
#endif
}

void wavelet_apply(WaveletType type, int order, ExtensionType ext,
                   const float *__restrict src, size_t length,
                   float *__restrict desthi, float *__restrict destlo) {
  switch (order) {
    case 2:
      wavelet_apply2(type, ext, src, length, desthi, destlo);
      break;
    case 4:
      wavelet_apply4(type, ext, src, length, desthi, destlo);
      break;
    case 6:
      wavelet_apply6(type, ext, src, length, desthi, destlo);
      break;
    case 8:
      wavelet_apply8(type, ext, src, length, desthi, destlo);
      break;
    case 12:
      wavelet_apply12(type, ext, src, length, desthi, destlo);
      break;
    case 16:
      wavelet_apply16(type, ext, src, length, desthi, destlo);
      break;
    default:
      // TODO(v.markovtsev): implement universal SIMD version
      wavelet_apply_na(type, order, ext, src, length, desthi, destlo);
      break;
  }
}

void stationary_wavelet_apply(WaveletType type, int order, int level,
                              ExtensionType ext,
                              const float *__restrict src, size_t length,
                              float *__restrict desthi,
                              float *__restrict destlo) {
  int size = order * (1 << (level - 1));
  switch (size) {
    case 2:
      stationary_wavelet_apply2(type, level, ext, src, length, desthi, destlo);
      break;
    case 4:
      stationary_wavelet_apply4(type, level, ext, src, length, desthi, destlo);
      break;
    case 6:
      stationary_wavelet_apply6(type, level, ext, src, length, desthi, destlo);
      break;
    case 8:
      stationary_wavelet_apply8(type, level, ext, src, length, desthi, destlo);
      break;
    case 12:
      stationary_wavelet_apply12(type, level, ext, src, length, desthi, destlo);
      break;
    case 16:
      stationary_wavelet_apply16(type, level, ext, src, length, desthi, destlo);
      break;
    case 24:
      stationary_wavelet_apply24(type, level, ext, src, length, desthi, destlo);
      break;
    default:
      stationary_wavelet_applyN(type, size, level, ext, src, length,
                                desthi, destlo);
      break;
  }
}
