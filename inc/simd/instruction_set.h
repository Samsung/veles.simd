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

#ifdef __SSE3__
#ifdef __AVX__
#include <immintrin.h>
#define __m256_get_ps(vec, index) vec[index]
#else
#include <simd/avxintrin-emu.h>
#define __AVX__
#endif
#elif defined(__ARM_NEON__)
#include <arm_neon.h>
#endif

#ifdef __AVX__

#ifndef __xgetbv
static __attribute__((always_inline)) inline unsigned long long __xgetbv() {
#if defined(__GNUC__) && __GNUC_PREREQ(4, 4)
  unsigned int index = 0;
  unsigned int eax, edx;
  __asm__ __volatile__("xgetbv" : "=a"(eax), "=d"(edx) : "c"(index));
#ifdef __cplusplus
  return (static_cast<unsigned long long>(edx) << 32) | eax;
#else
  return ((unsigned long long)edx << 32) | eax;
#endif
#else
  return 0;
#endif
}
#endif  // __xgetbv

#if defined(__cplusplus) && __GNUC__ == 4 && __GNUC_MINOR__ < 8

__attribute__((optimize(2), always_inline)) inline float
_mm256_get_ps(__m256 vector, int index) {
  union __m128_buggy_gxx_up_to_4_7 {
      __m256 v;
      float e[8];
  };
  __m128_buggy_gxx_up_to_4_7 vfix;
  vfix.v = vector;
  return vfix.e[index];
}

#else

#define _mm256_get_ps(vec, index) __m256_get_ps(vec, index)

#endif  // defined(__cplusplus) && __GNUC__ == 4 && __GNUC_MINOR__ < 8

#endif

#endif  // INC_SIMD_INSTRUCTION_SET_H_
