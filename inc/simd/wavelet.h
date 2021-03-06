/*! @file wavelet.h
 *  @brief Wavelet low level functions.
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

#ifndef INC_SIMD_WAVELET_H_
#define INC_SIMD_WAVELET_H_

#include <stddef.h>
#include <simd/common.h>
#include <simd/attributes.h>
#include <simd/wavelet_types.h>

SIMD_API_BEGIN

/// @brief Checks if the specified wavelet order is supported by the engine.
/// @param type The wavelet type.
/// @param order The order of the wavelet to check.
/// @return 1 if order is valid, otherwise 0.
int wavelet_validate_order(WaveletType type, int order);

/// @brief If AVX SIMD instruction set is used, copies the source signal to
/// a new array with specific excessive aligned format, else just returns src.
/// @param order The number of applied wavelet's coefficients.
/// @param src The source signal.
/// @param length The length of source in float-s.
/// @return A newly allocated memory block which should be disposed with
/// free() if AVX SIMD is activated; otherwise, src.
/// @note It should be used with wavelet_apply() or wavelet_apply_na().
float *wavelet_prepare_array(int order, const float *src, size_t length)
    NOTNULL(2)
#ifdef __AVX__
    MALLOC WARN_UNUSED_RESULT
#endif
;  // NOLINT(*)

/// @brief Allocates the array which is capable of storing one of the two
/// wavelet splitted parts of a signal (that is, desthi or destlo in
/// wavelet_apply()).
/// @param order The number of applied wavelet's coefficients.
/// @param sourceLength The logical length of the source array
/// (in float-s, not in bytes).
/// @return A newly allocated memory block which should be disposed with
/// free().
/// @note It should be used with wavelet_apply() or wavelet_apply_na().
float *wavelet_allocate_destination(int order, size_t sourceLength)
    MALLOC WARN_UNUSED_RESULT;

/// @brief Splits src into four regions to reuse the allocated memory.
/// @note It should be used with wavelet_apply() or wavelet_apply_na().
void wavelet_recycle_source(int order, float *src, size_t length,
                            float **desthihi, float **desthilo,
                            float **destlohi, float **destlolo)
    NOTNULL(2, 4, 5, 6, 7);

/// @brief Performs a single wavelet transform on series of real numbers.
/// @param type The wavelet type.
/// @param order The order of the wavelet to apply.
/// For example, order = 6 means 6 coefficients.
/// @param ext The way to extend the signal.
/// @param src An array of floating point numbers to transform.
/// @param length The logical length of src (in float-s, not in bytes).
/// @param desthi The high frequency part of result (highpass). It must be at
/// least of size length/2.
/// @param destlo The low frequency part of result (lowpass). It must be at
/// least of size length/2.
/// @details Daubechies wavelet of order 8 is used by default.
/// @pre length must be greater than or equal to 8.
/// @pre length must be even.
void wavelet_apply(WaveletType type, int order, ExtensionType ext,
                   const float *__restrict src, size_t length,
                   float *__restrict desthi, float *__restrict destlo)
    NOTNULL(4, 6, 7);

/// @brief Performs a single wavelet transform on series of real numbers.
/// (no SIMD acceleration is used).
/// @param type The wavelet type.
/// @param order The order of the wavelet to apply.
/// For example, order = 6 means 6 coefficients.
/// @param ext The way to extend the signal.
/// @param src An array of floating point numbers to transform.
/// @param length The logical length of src (in float-s, not in bytes).
/// @param desthi The high frequency part of result (highpass). It must be at
/// least of size length/2.
/// @param destlo The low frequency part of result (lowpass). It must be at
/// least of size length/2.
/// @pre length must be greater than or equal to (order * 2).
/// @pre length must be even.
void wavelet_apply_na(WaveletType type, int order, ExtensionType ext,
                      const float *__restrict src, size_t length,
                      float *__restrict desthi, float *__restrict destlo)
    NOTNULL(4, 6, 7);

/// @brief Performs a single stationary (undecimated) wavelet transform
/// on series of real numbers.
/// @param type The wavelet type.
/// @param order The order of the wavelet to apply.
/// For example, order = 6 means 6 coefficients.
/// @param level The current decomposition level.
/// @param ext The way to extend the signal.
/// @param src An array of floating point numbers to transform.
/// @param length The logical length of src (in float-s, not in bytes).
/// @param desthi The high frequency part of result (highpass). It must be at
/// least of size length.
/// @param destlo The low frequency part of result (lowpass). It must be at
/// least of size length.
/// @details Daubechies wavelet of order 8 is used by default.
/// @pre length must be greater than or equal to 8.
/// @pre length must be even.
void stationary_wavelet_apply(WaveletType type, int order, int level,
                              ExtensionType ext,
                              const float *__restrict src, size_t length,
                              float *__restrict desthi,
                              float *__restrict destlo)
    NOTNULL(5, 7, 8);

/// @brief Performs a single stationary (undecimated) wavelet transform
/// on series of real numbers (no SIMD acceleration is used).
/// @param type The wavelet type.
/// @param order The order of the wavelet to apply.
/// For example, order = 6 means 6 coefficients.
/// @param level The current decomposition level.
/// @param ext The way to extend the signal.
/// @param src An array of floating point numbers to transform.
/// @param length The logical length of src (in float-s, not in bytes).
/// @param desthi The high frequency part of result (highpass). It must be at
/// least of size length.
/// @param destlo The low frequency part of result (lowpass). It must be at
/// least of size length.
/// @pre length must be greater than or equal to (order * 2).
/// @pre length must be even.
void stationary_wavelet_apply_na(WaveletType type, int order, int level,
                                 ExtensionType ext,
                                 const float *__restrict src, size_t length,
                                 float *__restrict desthi,
                                 float *__restrict destlo)
    NOTNULL(5, 7, 8);

SIMD_API_END

#endif  // INC_SIMD_WAVELET_H_
