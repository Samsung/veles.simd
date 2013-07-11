/*! @file matrix.h
 *  @brief New file description.
 *  @author Markovtsev Vadim <v.markovtsev@samsung.com>
 *  @version 1.0
 *
 *  @section Notes
 *  This code partially conforms to <a href="http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml">Google C++ Style Guide</a>.
 *
 *  @section Copyright
 *  Copyright 2013 Samsung R&D Institute Russia
 */

#include <stddef.h>
#include <simd/attributes.h>

#ifndef INC_SIMD_MATRIX_H_
#define INC_SIMD_MATRIX_H_

#ifdef __cplusplus
extern "C" {
#endif

void matrix_multiply(int simd, const float *m1, const float *m2,
                     size_t w1, size_t h1, size_t w2, size_t h2,
                     float *res) NOTNULL(2,3,8);

void matrix_multiply_transposed(int simd, const float *m1, const float *m2,
                                size_t w1, size_t h1, size_t w2, size_t h2,
                                float *res) NOTNULL(2,3,8);

#ifdef __cplusplus
}
#endif

#endif  // INC_SIMD_MATRIX_H_
