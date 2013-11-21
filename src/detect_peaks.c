/*! @file detect_peaks.c
 *  @brief Find extrema in 1D signal.
 *  @author Vadim Markovtsev <v.markovtsev@samsung.com>
 *  @version 1.0
 *
 *  @section Notes
 *  This code partially conforms to <a href="http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml">Google C++ Style Guide</a>.
 *
 *  @section Copyright
 *  Copyright 2013 Samsung R&D Institute Russia
 */

#include "inc/simd/detect_peaks.h"
#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <simd/instruction_set.h>

INLINE void append_peak(int position, float value, int index,
                        ExtremumPoint **results, size_t *allocatedSize) {
  int alloc_size = *allocatedSize;
  if (alloc_size > index) {
    (*results)[index] = (ExtremumPoint) { .position = position,
                                          .value = value };
    return;
  }
  if (alloc_size == 0) {
    alloc_size = 2;
  } else if (alloc_size < (INT_MAX >> 1)) {
    alloc_size <<= 1;
  } else {
    // No more place. This is weird anyway.
    return;
  }
  *allocatedSize = alloc_size;
  *results = realloc(*results, alloc_size * sizeof((*results)[0]));
  (*results)[index] = (ExtremumPoint) { .position = position,
                                          .value = value };
}

INLINE void check_peak(const float *data, int index, ExtremumType type,
                       ExtremumPoint **results, size_t *resultsLength,
                       size_t *allocatedSize) {
  float prev = data[index - 1];
  float curr = data[index];
  float next = data[index + 1];
  float delta1 = curr - prev;
  float delta2 = curr - next;
  if (delta1 * delta2 > 0) {
    if ((delta1 > 0 && (type & kExtremumTypeMaximum) != 0) ||
        (delta1 < 0 && (type & kExtremumTypeMinimum) != 0)) {
      append_peak(index, curr, *resultsLength, results, allocatedSize);
      (*resultsLength)++;
    }
  }
}

void detect_peaks(int simd, const float *data, size_t size, ExtremumType type,
                  ExtremumPoint **results, size_t *resultsLength) {
  assert(data);
  assert(results);
  assert(resultsLength);
  assert(size > 2);
  *resultsLength = 0;
  *results = NULL;
  size_t allocated_size = 0;
  int isize = (int)size;

  if (simd) {
#ifdef __ARM_NEON__
    for (int i = 0; i < isize - 4; i += 4) {
      float32x4_t vec1 = vld1q_f32(data + i);
      float32x4_t vec2 = vld1q_f32(data + i + 1);
      float32x4_t max = vmaxq_f32(vec1, vec2);
      uint32x4_t cmpvec1 = vceqq_f32(max, vec1);
      uint32x4_t cmpvec2 = vceqq_f32(max, vec2);
      uint64x2_t cmpvec1_64 = vpaddlq_u32(cmpvec1);
      uint64x2_t cmpvec2_64 = vpaddlq_u32(cmpvec2);
      if (vgetq_lane_u64(cmpvec1_64, 0) == 0x1FFFFFFFE &&
          vgetq_lane_u64(cmpvec1_64, 1) == 0x1FFFFFFFE) {
        continue;
      }
      if (vgetq_lane_u64(cmpvec2_64, 0) == 0x1FFFFFFFE &&
          vgetq_lane_u64(cmpvec2_64, 1) == 0x1FFFFFFFE) {
        check_peak(data, i + 4, type, results, resultsLength, &allocated_size);
        continue;
      }
      for (int j = i + 1; j < i + 5; j++) {
        check_peak(data, j, type, results, resultsLength, &allocated_size);
      }
    }
    for (int i = (isize - 1) & ~0x3; i < isize - 1; i++) {
      check_peak(data, i, type, results, resultsLength, &allocated_size);
    }
  } else {
#elif defined(__AVX__)
    for (int i = 0; i < isize - 8; i += 8) {
      __m256 vec1 = _mm256_loadu_ps(data + i);
      __m256 vec2 = _mm256_loadu_ps(data + i + 1);
      __m256 max = _mm256_max_ps(vec1, vec2);
      __m256 cmpvec = _mm256_cmp_ps(max, vec1, _CMP_EQ_UQ);
      int cmpres = _mm256_movemask_ps(cmpvec);
      if (cmpres == 0xFF) {
        continue;
      }
      cmpvec = _mm256_cmp_ps(max, vec2, _CMP_EQ_UQ);
      cmpres = _mm256_movemask_ps(cmpvec);
      if (cmpres == 0xFF) {
        check_peak(data, i + 8, type, results, resultsLength, &allocated_size);
        continue;
      }
      for (int j = i + 1; j < i + 9; j++) {
        check_peak(data, j, type, results, resultsLength, &allocated_size);
      }
    }
    for (int i = (isize - 1) & ~0x7; i < isize - 1; i++) {
      check_peak(data, i, type, results, resultsLength, &allocated_size);
    }
  } else {
#else
  } {
#endif
    for (int i = 1; i < isize - 1; i++) {
      check_peak(data, i, type, results, resultsLength, &allocated_size);
    }
  }
}
