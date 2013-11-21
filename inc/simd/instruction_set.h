/*! @file instruction_set.h
 *  @brief Selects the right include file for the available SIMD instruction set.
 *  @author Vadim Markovtsev <v.markovtsev@samsung.com>
 *  @version 1.0
 *
 *  @section Notes
 *  This code partially conforms to <a href="http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml">Google C++ Style Guide</a>.
 *
 *  @section Copyright
 *  Copyright 2013 Samsung R&D Institute Russia
 */

#ifndef INC_SIMD_INSTRUCTION_SET_H_
#define INC_SIMD_INSTRUCTION_SET_H_

#include <simd/common.h>

#ifdef __SSE3__
#ifdef __AVX__
#include <immintrin.h>
#define __m256_get_ps(vec, index) vec[index]
#else
#include <simd/avxintrin-emu.h>
#define __AVX__
#endif
#ifdef __cplusplus

SIMD_API_BEGIN

#ifdef __AVX__

float _mm256_get_ps(__m256 vector, int index) __attribute__((optimize(2)));

static __attribute__((always_inline)) inline unsigned long long __xgetbv() {
#if defined(__GNUC__) && __GNUC_PREREQ(4, 4)
  unsigned int index = 0;
  unsigned int eax, edx;
  __asm__ __volatile__("xgetbv" : "=a"(eax), "=d"(edx) : "c"(index));
  return ((unsigned long long)edx << 32) | eax;
#else
  return 0;
#endif
}

#endif  // __AVX__

SIMD_API_END
#else
#ifndef _mm256_get_ps_IMPLEMENTATION
#define _mm256_get_ps(vec, index) __m256_get_ps(vec, index)
#endif  // _mm256_get_ps_IMPLEMENTATION
#endif  // __cplusplus
#elif defined(__ARM_NEON__)
#include <arm_neon.h>
#endif

#endif  // INC_SIMD_INSTRUCTION_SET_H_
