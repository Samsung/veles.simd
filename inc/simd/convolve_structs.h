/*! @file convolute_structs.h
 *  @brief New file description.
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

#ifndef INC_SIMD_CONVOLVE_STRUCTS_H_
#define INC_SIMD_CONVOLVE_STRUCTS_H_

#include <stddef.h>
#include <simd/common.h>

SIMD_API_BEGIN

struct ConvolutionOverlapSaveHandle {
  void *fft_plan;
  void *fft_inverse_plan;
  float *fft_boiler_plate;
  float *H;
  size_t x_length;
  size_t h_length;
  int *L;
  int reverse;
};

struct ConvolutionFFTHandle {
  void *fft_plan;
  void *fft_inverse_plan;
  int *M;
  int x_length;
  int h_length;
  float **inputs;
  int reverse;
};

typedef enum {
  kConvolutionAlgorithmBruteForce,
  kConvolutionAlgorithmFFT,
  kConvolutionAlgorithmOverlapSave
} ConvolutionAlgorithm;

struct ConvolutionHandle {
  ConvolutionAlgorithm algorithm;
  int x_length;
  int h_length;
  union {
    struct ConvolutionFFTHandle fft;
    struct ConvolutionOverlapSaveHandle os;
  } handle;
};

SIMD_API_END

#endif  // INC_SIMD_CONVOLVE_STRUCTS_H_
