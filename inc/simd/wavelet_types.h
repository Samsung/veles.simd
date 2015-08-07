/*! @file wavelet_types.h
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

#ifndef INC_SIMD_WAVELET_TYPES_H_
#define INC_SIMD_WAVELET_TYPES_H_

#include <simd/common.h>

SIMD_API_BEGIN

typedef enum {
  WAVELET_TYPE_DAUBECHIES,
  WAVELET_TYPE_COIFLET,
  WAVELET_TYPE_SYMLET
} WaveletType;

typedef enum {
  /// @brief 1 2 3 | 1 2 3
  EXTENSION_TYPE_PERIODIC,
  /// @brief 1 2 3 | 3 2 1
  EXTENSION_TYPE_MIRROR,
  /// @brief 1 2 3 | 3 3 3
  EXTENSION_TYPE_CONSTANT,
  /// @brief 1 2 3 | 0 0 0
  EXTENSION_TYPE_ZERO,
} ExtensionType;

SIMD_API_END


#endif  // INC_SIMD_WAVELET_TYPES_H_
