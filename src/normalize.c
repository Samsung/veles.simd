/*! @file normalize.c
 *  @brief New file description.
 *  @author Vadim Markovtsev <v.markovtsev@samsung.com>
 *  @version 1.0
 *
 *  @section Notes
 *  This code partially conforms to <a href="http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml">Google C++ Style Guide</a>.
 *
 *  @section Copyright
 *  Copyright © 2013 Samsung R&D Institute Russia
 *
 *  @section License
 *  Licensed to the Apache Software Foundation (ASF) under one
 *  or more contributor license agreements.  See the NOTICE file
 *  distributed with this work for additional information
 *  regarding copyright ownership.  The ASF licenses this file
 *  to you under the Apache License, Version 2.0 (the
 *  "License"); you may not use this file except in compliance
 *  with the License.  You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an
 *  "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 *  KIND, either express or implied.  See the License for the
 *  specific language governing permissions and limitations
 *  under the License.
 */

#define IMPLEMENTATION
#include "inc/simd/normalize.h"
#include <assert.h>
#include <float.h>
#include <simd/instruction_set.h>
#include <simd/memory.h>

#ifdef __ARM_NEON__

static void normalize2D_minmax_neon(uint8_t min, uint8_t max,
                                    const uint8_t* src, int src_stride,
                                    int width, int height,
                                    float* dst, int dst_stride) {
  if (max == min) {
    memsetf(dst, 0, width * height);
    return;
  }
  const uint8x16_t min_vec = vdupq_n_u8(min);
  float diff = (max - min) / 2.f;
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

static void minmax2D_neon(const uint8_t* src, int src_stride,
                          int width, int height,
                          uint8_t* min_ptr, uint8_t* max_ptr) {
  uint8_t min = src[0], max = src[0];
  uint8x16_t min_vec = vdupq_n_u8(min), max_vec = vdupq_n_u8(max);
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width - 15; x += 16) {
      uint8x16_t vec = vld1q_u8(src + y * src_stride + x);
      if (min_ptr) {
        min_vec = vminq_u8(vec, min_vec);
      }
      if (max_ptr) {
        max_vec = vmaxq_u8(vec, max_vec);
      }
    }
    for (int x = width & ~0xF; x < width; x++) {
      float val = src[y * src_stride + x];
      if (min_ptr && val < min) {
        min = val;
      }
      if (max_ptr && val > max) {
        max = val;
      }
    }
  }
  // Gather the results
  uint8_t min_arr[16] __attribute__((aligned(64))),
      max_arr[16] __attribute__((aligned(64)));
  if (min_ptr) {
    vst1q_u8(min_arr, min_vec);
  }
  if (max_ptr) {
    vst1q_u8(max_arr, max_vec);
  }
  for (int i = 0; i < 16; i++) {
    float val = min_arr[i];
    if (min_ptr && val < min) {
      min = val;
    }
    val = max_arr[i];
    if (max_ptr && val > max) {
      max = val;
    }
  }

  if (min_ptr) {
    *min_ptr = min;
  }
  if (max_ptr) {
    *max_ptr = max;
  }
}

static void minmax1D_neon(const float* src, int length,
                          float* min_ptr, float* max_ptr) {
  float min = src[0], max = src[0];
  float32x4_t min_vec = vdupq_n_f32(min), max_vec = vdupq_n_f32(max);
  for (int i = 0; i < length - 3; i += 4) {
    float32x4_t vec = vld1q_f32(src + i);
    if (min_ptr) {
      min_vec = vminq_f32(vec, min_vec);
    }
    if (max_ptr) {
      max_vec = vmaxq_f32(vec, max_vec);
    }
  }
  for (int i = length & ~0x3; i < length; i++) {
    float val = src[i];
    if (min_ptr && val < min) {
      min = val;
    }
    if (max_ptr && val > max) {
      max = val;
    }
  }

  // Gather the results
  float min_arr[4] __attribute__((aligned(64))),
      max_arr[4] __attribute__((aligned(64)));
  if (min_ptr) {
    vst1q_f32(min_arr, min_vec);
  }
  if (max_ptr) {
    vst1q_f32(max_arr, max_vec);
  }
  for (int i = 0; i < 4; i++) {
    float val = min_arr[i];
    if (min_ptr && val < min) {
      min = val;
    }
    val = max_arr[i];
    if (max_ptr && val > max) {
      max = val;
    }
  }

  if (min_ptr) {
    *min_ptr = min;
  }
  if (max_ptr) {
    *max_ptr = max;
  }
}

#endif


#ifdef __SSE2__

static void normalize2D_minmax_sse(uint8_t min, uint8_t max,
                                   const uint8_t* src, int src_stride,
                                   int width, int height,
                                   float* dst, int dst_stride) {
  if (max == min) {
    memsetf(dst, 0, width * height);
    return;
  }
  const __m128i min_vec = _mm_set1_epi8(min);
  float diff = (max - min) / 2.f;
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
    }
    for (int x = width & ~0xF; x < width; x++) {
      dst[y * dst_stride + x] = (src[y * src_stride + x] - min) / diff - 1.0f;
    }
  }
}

static void minmax2D_sse(const uint8_t* src, int src_stride,
                         int width, int height,
                         uint8_t* min_ptr, uint8_t* max_ptr) {
  assert(width > 0 && height > 0);
  uint8_t min = src[0], max = src[0];
  __m128i min_vec = _mm_set1_epi8(min), max_vec = _mm_set1_epi8(max);
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width - 15; x += 16) {
      __m128i vec = _mm_loadu_si128((const __m128i*)(src + y * src_stride + x));
      if (min_ptr) {
        min_vec = _mm_min_epu8(vec, min_vec);
      }
      if (max_ptr) {
        max_vec = _mm_max_epu8(vec, max_vec);
      }
    }
    for (int x = width & ~0xF; x < width; x++) {
      float val = src[y * src_stride + x];
      if (min_ptr && val < min) {
        min = val;
      }
      if (max_ptr && val > max) {
        max = val;
      }
    }
  }
  // Gather the results
  uint8_t min_arr[16] __attribute__((aligned(64))),
      max_arr[16] __attribute__((aligned(64)));
  if (min_ptr) {
    _mm_store_si128((__m128i*)min_arr, min_vec);
  }
  if (max_ptr) {
    _mm_store_si128((__m128i*)max_arr, max_vec);
  }
  for (int i = 0; i < 16; i++) {
    if (min_ptr) {
      float val = min_arr[i];
      if (val < min) {
        min = val;
      }
    }
    if (max_ptr) {
      float val = max_arr[i];
      if (val > max) {
        max = val;
      }
    }
  }

  if (min_ptr) {
    *min_ptr = min;
  }
  if (max_ptr) {
    *max_ptr = max;
  }
}

#ifdef __AVX__
static void minmax1D_avx(const float* src, int length,
                         float* min_ptr, float* max_ptr) {
  assert(length > 0);
  float min = src[0], max = src[0];
  __m256 min_vec = _mm256_set1_ps(min), max_vec = _mm256_set1_ps(max);
  for (int i = 0; i < length - 7; i += 8) {
    __m256 vec = _mm256_loadu_ps(src + i);
    if (min_ptr) {
      min_vec = _mm256_min_ps(vec, min_vec);
    }
    if (max_ptr) {
      max_vec = _mm256_max_ps(vec, max_vec);
    }
  }
  for (int i = length & ~0x7; i < length; i++) {
    float val = src[i];
    if (min_ptr && val < min) {
      min = val;
    }
    if (max_ptr && val > max) {
      max = val;
    }
  }

  // Gather the results
  float min_arr[8] __attribute__((aligned(64))),
      max_arr[8] __attribute__((aligned(64)));
  if (min_ptr) {
    _mm256_store_ps(min_arr, min_vec);
  }
  if (max_ptr) {
    _mm256_store_ps(max_arr, max_vec);
  }
  for (int i = 0; i < 8; i++) {
    if (min_ptr) {
      float val = min_arr[i];
      if (val < min) {
        min = val;
      }
    }
    if (max_ptr) {
      float val = max_arr[i];
      if (val > max) {
        max = val;
      }
    }
  }

  if (min_ptr) {
    *min_ptr = min;
  }
  if (max_ptr) {
    *max_ptr = max;
  }
}
#endif  // __AVX__

#endif  // __SSE2__

static void normalize2D_minmax_novec(uint8_t min, uint8_t max,
                                     const uint8_t* src, int src_stride,
                                     int width, int height,
                                     float* dst, int dst_stride) {
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

static void minmax2D_novec(const uint8_t* src, int src_stride,
                           int width, int height,
                           uint8_t* min_ptr, uint8_t* max_ptr) {
  uint8_t min = src[0], max = src[0];
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      float val = src[y * src_stride + x];
      if (min_ptr && val < min) {
        min = val;
      }
      if (max_ptr && val > max) {
        max = val;
      }
    }
  }
  if (min_ptr) {
    *min_ptr = min;
  }
  if (max_ptr) {
    *max_ptr = max;
  }
}

static void minmax1D_novec(const float* src, int length,
                           float* min_ptr, float* max_ptr) {
  float min = src[0], max = src[0];
  for (int i = 0; i < length; i++) {
    float val = src[i];
    if (min_ptr && val < min) {
      min = val;
    }
    if (max_ptr && val > max) {
      max = val;
    }
  }
  if (min_ptr) {
    *min_ptr = min;
  }
  if (max_ptr) {
    *max_ptr = max;
  }
}

void normalize2D(int simd, const uint8_t* src, int src_stride,
                 int width, int height, float* dst, int dst_stride) {
  uint8_t min, max;
  minmax2D(simd, src, src_stride, width, height, &min, &max);
  normalize2D_minmax(simd, min, max, src, src_stride, width, height,
                     dst, dst_stride);
}

void minmax2D(int simd, const uint8_t* src, int src_stride,
              int width, int height, uint8_t* min, uint8_t* max) {
  assert(src);
  assert(width > 0);
  assert(height > 0);
  assert(src_stride >= width);
  if (!min && !max) {
    return;
  }
  if (simd) {
#ifdef __ARM_NEON__
    minmax2D_neon(src, src_stride, width, height, min, max);
  } else {
#elif defined(__SSE2__)
    minmax2D_sse(src, src_stride, width, height, min, max);
  } else {
#else
  } {
#endif
    minmax2D_novec(src, src_stride, width, height, min, max);
  }
}

void normalize2D_minmax(int simd, uint8_t min, uint8_t max,
                        const uint8_t* src, int src_stride,
                        int width, int height, float* dst, int dst_stride) {
  assert(src);
  assert(dst);
  assert(width > 0);
  assert(height > 0);
  assert(src_stride >= width);
  assert(dst_stride >= width);
  assert(min <= max);
  if (simd) {
#ifdef __ARM_NEON__
    normalize2D_minmax_neon(min, max, src, src_stride, width, height,
                            dst, dst_stride);
  } else {
#elif defined(__SSE2__)
    normalize2D_minmax_sse(min, max, src, src_stride, width, height,
                           dst, dst_stride);
  } else {
#else
  } {
#endif
    normalize2D_minmax_novec(min, max, src, src_stride, width, height,
                             dst, dst_stride);
  }
}

void minmax1D(int simd, const float *src, int length, float *min, float *max) {
  assert(src);
  assert(length > 0);
  if (!min && !max) {
    return;
  }
  if (simd) {
#ifdef __ARM_NEON__
    minmax1D_neon(src, length, min, max);
  } else {
#elif defined(__SSE2__)
    minmax1D_avx(src, length, min, max);
  } else {
#else
  } {
#endif
    minmax1D_novec(src, length, min, max);
  }
}
