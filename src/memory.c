/*! @file memory.cc
 *  @brief Memory routines with SIMD optimization.
 *  @author Markovtsev Vadim <v.markovtsev@samsung.com>
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

#define LIBSIMD_IMPLEMENTATION
#include "inc/simd/memory.h"
#include <assert.h>
#include <simd/instruction_set.h>
#ifndef __USE_XOPEN2K
#define __USE_XOPEN2K
#endif
#include <stdlib.h>
#include <stdint.h>

#ifdef __AVX__
static int align_offset_internal(const void *ptr) {
  uintptr_t addr = (uintptr_t)ptr;
  if ((addr & 31) != 0) {
    return 32 - (addr % 32);
  }
  return 0;
}

int align_complement_f32(const float *ptr) {
  return align_offset_internal(ptr) / 4;
}

int align_complement_i16(const int16_t *ptr) {
  return align_offset_internal(ptr) / 2;
}

int align_complement_u16(const uint16_t *ptr) {
  return align_offset_internal(ptr) / 2;
}

int align_complement_i32(const int32_t *ptr) {
  return align_offset_internal(ptr) / 4;
}

int align_complement_u32(const uint32_t *ptr) {
  return align_offset_internal(ptr) / 4;
}
#endif

void *malloc_aligned_offset(size_t size, int offset) {
  assert(offset >= 0 && offset < 32);
  void *ptr = malloc_aligned(size + offset);
  return (char *)ptr + offset;
}

void *malloc_aligned(size_t size) {
  void *ptr;
#ifndef __ANDROID__
  if (posix_memalign(&ptr, 64, size) != 0) {
    return NULL;
  }
#else
  ptr = memalign(64, size);
#endif
  return ptr;
}

float *mallocf(size_t length) {
  return malloc_aligned(length * sizeof(float));
}

void memsetf(float *ptr, float value, size_t length) {
#ifdef __AVX__
  const __m256 fillvec = _mm256_set1_ps(value);
  size_t startIndex = align_complement_f32(ptr);

  for (size_t i = 0; i < startIndex; i++) {
    ptr[i] = value;
  }

  for (int i = (int)startIndex; i < (int)length - 7; i += 8) {
    _mm256_store_ps(ptr + i, fillvec);
  }

  for (size_t i = startIndex + ((length - startIndex) & ~0x7);
      i < length; i++) {
    ptr[i] = value;
  }
#elif defined(__ARM_NEON__)
  const float32x4_t fillvec = vdupq_n_f32(value);
  for (int i = 0; i < (int)length - 3; i += 4) {
    vst1q_f32(ptr + i, fillvec);
  }
  for (size_t i = (length & ~0x3); i < length; i++) {
    ptr[i] = value;
  }
#else
  for (size_t i = 0; i < length; i++) {
    ptr[i] = value;
  }
#endif
}

float *zeropadding(const float *ptr, size_t length, size_t *newLength) {
  return zeropaddingex(ptr, length, newLength, 0);
}

float *zeropaddingex(const float *ptr, size_t length, size_t *newLength,
                     size_t additionalLength) {
  size_t nl = length;
  int log = 2;
  while (nl >>= 1) {
    log++;
  }
  nl = (1 << log);
  *newLength = nl;
  float *ret = mallocf(nl + additionalLength);
  memcpy(ret, ptr, length * sizeof(float));
  memsetf(ret + length, 0.f, nl - length);
  return ret;
}

float *rmemcpyf(float *__restrict dest,
                const float *__restrict src, size_t length) {
#ifdef __AVX__
  for (int i = 0; i < (int)length - 7; i += 8) {
    __m256 vec = _mm256_loadu_ps(src + i);
    vec = _mm256_permute2f128_ps(vec, vec, 1);
    vec = _mm256_permute_ps(vec, 0x1B);
    _mm256_storeu_ps(dest + length - i - 8, vec);
  }

  for (size_t i = (length & ~0x7); i < length; i++) {
    dest[length - i - 1] = src[i];
  }
#elif defined(__ARM_NEON__)
  for (int i = 0; i < (int)length - 3; i += 4) {
    float32x4_t vec = vld1q_f32(src + i);
    vec = vrev64q_f32(vec);
    vec = vcombine_f32(vget_high_f32(vec), vget_low_f32(vec));
    vst1q_f32(dest + length - i - 4, vec);
  }

  for (size_t i = (length & ~0x3); i < length; i++) {
    dest[length - i - 1] = src[i];
  }
#else
  for (size_t i = 0; i < length; i++) {
    dest[i] = src[length - i - 1];
  }
#endif
  return dest;
}

float *crmemcpyf(float *__restrict dest,
                 const float *__restrict src, size_t length) {
  for (size_t i = 0; i < length; i += 2) {
    dest[i] = src[length - i - 2];
    dest[i + 1] = src[length - i - 1];
  }
  return dest;
}
