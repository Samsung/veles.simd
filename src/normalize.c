/*! @file normalize.c
 *  @brief New file description.
 *  @author Vadim Markovtsev <v.markovtsev@samsung.com>
 *  @version 1.0
 *
 *  @section Notes
 *  This code partially conforms to <a href="http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml">Google C++ Style Guide</a>.
 *
 *  @section Copyright
 *  Copyright 2013 Samsung R&D Institute Russia
 */

#include "inc/simd/normalize.h"
#include <assert.h>
#ifdef __SSE2__
#include <immintrin.h>
#elif defined(__ARM_NEON__)
#include <arm_neon.h>
#endif
#include <simd/memory.h>

#ifdef __ARM_NEON__

static void normalize2D_neon(const uint8_t* src, int src_stride,
                             int width, int height,
                             float* dst, int dst_stride) {
  // Step 1 - get the maximum and minimum
  uint8_t min = src[0], max = src[0];
  uint8x16_t min_vec = vdupq_n_u8(min), max_vec = vdupq_n_u8(max);
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width - 15; x += 16) {
      uint8x16_t vec = vld1q_u8(src + y * src_stride + x);
      min_vec = vminq_u8(vec, min_vec);
      max_vec = vmaxq_u8(vec, max_vec);
    }
    for (int x = width & ~0xF; x < width; x++) {
      float val = src[y * src_stride + x];
      if (val < min) {
        min = val;
      } else if (val > max) {
        max = val;
      }
    }
  }
  // Gather the results
  uint8_t min_arr[16] __attribute__((aligned(64))),
      max_arr[16] __attribute__((aligned(64)));
  vst1q_u8(min_arr, min_vec);
  vst1q_u8(max_arr, max_vec);
  for (int i = 0; i < 16; i++) {
    float val = min_arr[i];
    if (val < min) {
      min = val;
    }
    val = max_arr[i];
    if (val > max) {
      max = val;
    }
  }

  // Step 2 - perform the normalization
  if (max == min) {
    memsetf(dst, width * height, 0);
    return;
  }
  float diff = (max - min) / 2.f;
  min_vec = vdupq_n_u8(min);
  const float32x4_t diff_vec = vdupq_n_f32(1.f / diff);
  const float32x4_t sub_vec = vdupq_n_f32(1.f);
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width - 15; x += 16) {
      uint8x16_t vec = vld1q_u8(src + y * src_stride + x);
      vec = vsubq_u8(vec, min_vec);
      uint8x8_t vec8lo = vget_low_u8(vec);
      uint8x8_t vec8hi = vget_high_u8(vec);
      uint16x8_t vec16lo = vmovl_u8(vec8lo);
      uint16x8_t vec16hi = vmovl_u8(vec8hi);

      uint16x4_t vec16lolo = vget_low_u16(vec16lo);
      uint16x4_t vec16lohi = vget_high_u16(vec16lo);
      uint16x4_t vec16hilo = vget_low_u16(vec16hi);
      uint16x4_t vec16hihi = vget_high_u16(vec16hi);

      uint32x4_t intlolo = vmovl_u16(vec16lolo);
      uint32x4_t intlohi = vmovl_u16(vec16lohi);
      uint32x4_t inthilo = vmovl_u16(vec16hilo);
      uint32x4_t inthihi = vmovl_u16(vec16hihi);

      float32x4_t flolo = vcvtq_f32_u32(intlolo);
      float32x4_t flohi = vcvtq_f32_u32(intlohi);
      float32x4_t fhilo = vcvtq_f32_u32(inthilo);
      float32x4_t fhihi = vcvtq_f32_u32(inthihi);

      flolo = vmulq_f32(flolo, diff_vec);
      flohi = vmulq_f32(flohi, diff_vec);
      fhilo = vmulq_f32(fhilo, diff_vec);
      fhihi = vmulq_f32(fhihi, diff_vec);

      flolo = vsubq_f32(flolo, sub_vec);
      flohi = vsubq_f32(flohi, sub_vec);
      fhilo = vsubq_f32(fhilo, sub_vec);
      fhihi = vsubq_f32(fhihi, sub_vec);

      float* dst_ptr = dst + y * dst_stride + x;
      vst1q_f32(dst_ptr, flolo);
      dst_ptr += 4;
      vst1q_f32(dst_ptr, flohi);
      dst_ptr += 4;
      vst1q_f32(dst_ptr, fhilo);
      dst_ptr += 4;
      vst1q_f32(dst_ptr, fhihi);
      dst_ptr += 4;
    }
    for (int x = width & ~0xF; x < width; x++) {
      dst[y * dst_stride + x] = (src[y * src_stride + x] - min) / diff - 1.0f;
    }
  }
}

#endif


#ifdef __SSE2__

static void normalize2D_sse(const uint8_t* src, int src_stride,
                            int width, int height,
                            float* dst, int dst_stride) {
  // Step 1 - get the maximum and minimum
  uint8_t min = src[0], max = src[0];
  __m128i min_vec = _mm_set1_epi8(min), max_vec = _mm_set1_epi8(max);
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width - 15; x += 16) {
      __m128i vec = _mm_loadu_si128((const __m128i*)(src + y * src_stride + x));
      min_vec = _mm_min_epu8(vec, min_vec);
      max_vec = _mm_max_epu8(vec, max_vec);
    }
    for (int x = width & ~0xF; x < width; x++) {
      float val = src[y * src_stride + x];
      if (val < min) {
        min = val;
      } else if (val > max) {
        max = val;
      }
    }
  }
  // Gather the results
  uint8_t min_arr[16] __attribute__((aligned(64))),
      max_arr[16] __attribute__((aligned(64)));
  _mm_store_si128((__m128i*)min_arr, min_vec);
  _mm_store_si128((__m128i*)max_arr, max_vec);
  for (int i = 0; i < 16; i++) {
    float val = min_arr[i];
    if (val < min) {
      min = val;
    }
    val = max_arr[i];
    if (val > max) {
      max = val;
    }
  }

  // Step 2 - perform the normalization
  if (max == min) {
    memsetf(dst, width * height, 0);
    return;
  }
  float diff = (max - min) / 2.f;
  min_vec = _mm_set1_epi8(min);
  const __m128 diff_vec = _mm_set1_ps(1.f / diff);
  const __m128 sub_vec = _mm_set1_ps(1.f);
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width - 15; x += 16) {
      __m128i vec = _mm_loadu_si128((const __m128i*)(src + y * src_stride + x));
      vec = _mm_subs_epu8(vec, min_vec);
      __m128i intlo = _mm_unpacklo_epi8(vec, _mm_set1_epi8(0));
      __m128i inthi = _mm_unpackhi_epi8(vec, _mm_set1_epi8(0));
      __m128i intlolo = _mm_unpacklo_epi16(intlo, _mm_set1_epi16(0));
      __m128i intlohi = _mm_unpackhi_epi16(intlo, _mm_set1_epi16(0));
      __m128i inthilo = _mm_unpacklo_epi16(inthi, _mm_set1_epi16(0));
      __m128i inthihi = _mm_unpackhi_epi16(inthi, _mm_set1_epi16(0));

      __m128 flolo = _mm_cvtepi32_ps(intlolo);
      __m128 flohi = _mm_cvtepi32_ps(intlohi);
      __m128 fhilo = _mm_cvtepi32_ps(inthilo);
      __m128 fhihi = _mm_cvtepi32_ps(inthihi);

      flolo = _mm_mul_ps(flolo, diff_vec);
      flohi = _mm_mul_ps(flohi, diff_vec);
      fhilo = _mm_mul_ps(fhilo, diff_vec);
      fhihi = _mm_mul_ps(fhihi, diff_vec);

      flolo = _mm_sub_ps(flolo, sub_vec);
      flohi = _mm_sub_ps(flohi, sub_vec);
      fhilo = _mm_sub_ps(fhilo, sub_vec);
      fhihi = _mm_sub_ps(fhihi, sub_vec);

      float* dst_ptr = dst + y * dst_stride + x;
      _mm_storeu_ps(dst_ptr, flolo);
      dst_ptr += 4;
      _mm_storeu_ps(dst_ptr, flohi);
      dst_ptr += 4;
      _mm_storeu_ps(dst_ptr, fhilo);
      dst_ptr += 4;
      _mm_storeu_ps(dst_ptr, fhihi);
      dst_ptr += 4;
    }
    for (int x = width & ~0xF; x < width; x++) {
      dst[y * dst_stride + x] = (src[y * src_stride + x] - min) / diff - 1.0f;
    }
  }
}

#endif

static void normalize2D_novec(const uint8_t* src, int src_stride,
                              int width, int height,
                              float* dst, int dst_stride) {
  // Step 1 - get the maximum and minimum
  uint8_t min = src[0], max = src[0];
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      float val = src[y * src_stride + x];
      if (val < min) {
        min = val;
      } else if (val > max) {
        max = val;
      }
    }
  }

  // Step 2 - perform the normalization
  if (max == min) {
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
        dst[y * dst_stride + x] = 0;
      }
    }
    return;
  }
  float diff = (max - min) / 2.f;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      dst[y * dst_stride + x] = (src[y * src_stride + x] - min) / diff - 1.0f;
    }
  }
}

void normalize2D(int simd, const uint8_t* src, int src_stride,
                 int width, int height, float* dst, int dst_stride) {
  assert(src);
  assert(dst);
  assert(width > 0);
  assert(height > 0);
  assert(src_stride >= width);
  assert(dst_stride >= width);
  if (simd) {
#ifdef __ARM_NEON__
    normalize2D_neon(src, src_stride, width, height, dst, dst_stride);
  } else {
#elif defined(__SSE2__)
    normalize2D_sse(src, src_stride, width, height, dst, dst_stride);
  } else {
#else
  } {
#endif
    normalize2D_novec(src, src_stride, width, height, dst, dst_stride);
  }
}
