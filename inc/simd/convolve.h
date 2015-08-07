/*! @file convolve.h
 *  @brief Calculates the linear convolution of two signals.
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

#ifndef INC_SIMD_CONVOLVE_H_
#define INC_SIMD_CONVOLVE_H_

#include <stddef.h>
#include <simd/common.h>
#include <simd/attributes.h>
#include <simd/convolve_structs.h>

SIMD_API_BEGIN

typedef struct ConvolutionFFTHandle ConvolutionFFTHandle;

/// @brief Prepares for the calculation of linear convolution of two signals
/// using the FFT method.
/// @param xLength The length of the first array in float-s.
/// @param hLength The length of the second array in float-s.
/// @return The handle for convolve_fft().
ConvolutionFFTHandle convolve_fft_initialize(size_t xLength, size_t hLength);

/// @brief Calculates the linear convolution of two signals using
/// the FFT method.
/// @param handle The structure obtained from convolve_fft_initialize().
/// @param x The first signal (long one).
/// @param h The second signal (short one).
/// @param result The resulting signal of length xLength.
/// @note result and x or h may be the same arrays.
void convolve_fft(ConvolutionFFTHandle handle,
                  const float *x, const float *h,
                  float *result) NOTNULL(2, 3, 4);

/// @brief Frees any resources allocated by convolve_fft_initialize().
/// @param handle The structure obtained from convolve_fft_initialize().
void convolve_fft_finalize(ConvolutionFFTHandle handle);

typedef struct ConvolutionOverlapSaveHandle ConvolutionOverlapSaveHandle;

/// @brief Prepares for the calculation of linear convolution of two signals
/// using the overlap-save method.
/// @param xLength The length of the first array in float-s.
/// @param hLength The length of the second array in float-s.
/// @return The handle for convolve_overlap_save().
ConvolutionOverlapSaveHandle convolve_overlap_save_initialize(
    size_t xLength, size_t hLength);

/// @brief Calculates the linear convolution of two signals using
/// the overlap-save method.
/// @param handle The structure obtained from convolve_overlap_save_initialize().
/// @param x The first signal (long one).
/// @param h The second signal (short one).
/// @param result The resulting signal of length xLength.
/// @note result and x or h may be the same arrays.
void convolve_overlap_save(ConvolutionOverlapSaveHandle handle,
                           const float *x,
                           const float *h,
                           float *result) NOTNULL(2, 3, 4);

/// @brief Frees any resources allocated by convolve_overlap_save_initialize().
/// @param handle The structure obtained from convolve_overlap_save_initialize().
void convolve_overlap_save_finalize(ConvolutionOverlapSaveHandle handle);

/// @brief Calculates the linear convolution of two signals using
/// the "brute force" method.
/// @param simd Value indicating whether to use SIMD acceleration or not.
/// @param x The first signal (long one).
/// @param xLength The length of the first array in float-s.
/// @param h The second signal (short one).
/// @param hLength The length of the second array in float-s.
/// @param result The resulting signal of length xLength.
void convolve_simd(int simd,
                   const float *__restrict x, size_t xLength,
                   const float *__restrict h, size_t hLength,
                   float *__restrict result) NOTNULL(2, 4, 6);

typedef struct ConvolutionHandle ConvolutionHandle;

/// @brief Prepares for the calculation of linear convolution of two signals
/// using the best method.
/// @param xLength The length of the first array in float-s.
/// @param hLength The length of the second array in float-s.
/// @return The handle for convolve().
ConvolutionHandle convolve_initialize(size_t xLength, size_t hLength);

/// @brief Calculates the linear convolution of two signals using
/// the best method.
/// @param handle The structure obtained from convolve_initialize().
/// @param x The first signal (longer).
/// @param h The second signal (shorter).
/// @param result The resulting signal of length xLength + hLength - 1.
void convolve(ConvolutionHandle handle,
              const float *__restrict x, const float *__restrict h,
              float *__restrict result) NOTNULL(2, 3, 4);

/// @brief Frees any resources allocated by convolve_overlap_initialize().
/// @param handle The structure obtained from convolve_overlap_initialize().
void convolve_finalize(ConvolutionHandle handle);

SIMD_API_END

#endif  // INC_SIMD_CONVOLVE_H_
