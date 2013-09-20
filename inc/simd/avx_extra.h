/*! @file avx_extra.h
 *  @brief AVX ElementAt.
 *  @author Markovtsev Vadim <v.markovtsev@samsung.com>
 *  @version 1.0
 *
 *  @section Notes
 *  This code partially conforms to <a href="http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml">Google C++ Style Guide</a>.
 *
 *  @section Copyright
 *  Copyright 2013 Samsung R&D Institute Russia
 */

#ifndef INC_SIMD_AVX_EXTRA_H_
#define INC_SIMD_AVX_EXTRA_H_

#ifdef __cplusplus
extern "C" {
#endif

#if __GNUC__ >= 4
#pragma GCC visibility push(default)
#endif

#ifdef __AVX__

#include <immintrin.h>
#include "src/config.h"

float ElementAt(__m256 vector, int index) __attribute__((optimize(2)));

INLINE unsigned long long __xgetbv() {
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

#if __GNUC__ >= 4
#pragma GCC visibility pop
#endif

#ifdef __cplusplus
}
#endif

#endif  // INC_SIMD_AVX_EXTRA_H_
