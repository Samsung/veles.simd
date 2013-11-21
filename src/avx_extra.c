/*! @file avx_extra.c
 *  @brief AVX _mm256_get_ps.
 *  @author Markovtsev Vadim <v.markovtsev@samsung.com>
 *  @version 1.0
 *
 *  @section Notes
 *  This code partially conforms to <a href="http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml">Google C++ Style Guide</a>.
 *
 *  @section Copyright
 *  Copyright 2013 Samsung R&D Institute Russia
 */

#define _mm256_get_ps_IMPLEMENTATION
#include "inc/simd/instruction_set.h"

#ifdef __AVX__
float _mm256_get_ps(__m256 vector, int index) {
  return __m256_get_ps(vector, index);
}
#endif

