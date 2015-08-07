/*! @file normalize.h
 *  @brief Image plane normalization routines.
 *  @author Vadim Markovtsev <v.markovtsev@samsung.com>
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

#ifndef INC_SIMD_NORMALIZE_H_
#define INC_SIMD_NORMALIZE_H_

#include <stdint.h>
#include <simd/common.h>

SIMD_API_BEGIN

/// @brief Performs the plane normalization [min, max] -> [-1, 1]. Minimum
/// and maximum is determined from the array.
/// @param simd Value indicating whether to use available SIMD acceleration.
/// @param src The source byte array, stored in row-major format.
/// @param src_stride The stride (the actual width) of the plane.
/// @param width The width of the plane.
/// @param height The height of the plane.
/// @param dst The resulting floating point array.
/// @param dst_stride The stride of dst.
void normalize2D(int simd, const uint8_t *src, int src_stride,
                 int width, int height, float *dst, int dst_stride)
    NOTNULL(2, 6);

/// @brief Finds the minimum and the maximum value in the specified array.
/// @param simd Value indicating whether to use available SIMD acceleration.
/// @param src The source byte array, stored in row-major format.
/// @param src_stride The stride (the actual width) of the plane.
/// @param width The width of the plane.
/// @param height The height of the plane.
/// @param min The pointer to the resulting minimum. If NULL, minimum is not
/// calculated.
/// @param max The pointer to the resulting maximum. If NULL, maximum is not
/// calculated.
void minmax2D(int simd, const uint8_t *src, int src_stride,
              int width, int height, uint8_t *min, uint8_t *max)
    NOTNULL(2);

/// @brief Performs the plane normalization [min, max] -> [-1, 1].
/// @param simd Value indicating whether to use available SIMD acceleration.
/// @param min The precalculated minimum value.
/// @param max The precalculated maximum value.
/// @param src The source byte array, stored in row-major format.
/// @param src_stride The stride (the actual width) of the plane.
/// @param width The width of the plane.
/// @param height The height of the plane.
/// @param dst The resulting floating point array.
/// @param dst_stride The stride of dst.
void normalize2D_minmax(int simd, uint8_t min, uint8_t max,
                        const uint8_t *src, int src_stride,
                        int width, int height, float *dst, int dst_stride)
    NOTNULL(4, 8);

/// @brief Finds the minimum and the maximum value in the specified array.
/// @param simd Value indicating whether to use available SIMD acceleration.
/// @param src The source floating point array.
/// @param length The size of the array (in float-s, not in bytes).
/// @param min The pointer to the resulting minimum. If NULL, minimum is not
/// calculated.
/// @param max The pointer to the resulting maximum. If NULL, maximum is not
/// calculated.
void minmax1D(int simd, const float *src, int length, float *min,
              float *max) NOTNULL(2);

SIMD_API_END

#endif  // INC_SIMD_NORMALIZE_H_
