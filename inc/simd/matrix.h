/*! @file matrix.h
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

#ifndef INC_SIMD_MATRIX_H_
#define INC_SIMD_MATRIX_H_

#include <stddef.h>
#include <simd/common.h>
#include <simd/attributes.h>

SIMD_API_BEGIN

/// @brief Sums two matrices.
/// @param simd Value which indicates whether to use SIMD acceleration or not.
/// @param m1 The first matrix.
/// @param m2 The second matrix.
/// @param w The width of the matrices (the number of columns).
/// @param h The height of the matrices (the number of rows).
/// @param res The resulting matrix of the same size.
void matrix_add(int simd, const float *m1, const float *m2,
                size_t w, size_t h, float *res) NOTNULL(2,3,6);
/// @brief Subtracts the second matrix from the first one.
/// @param simd Value which indicates whether to use SIMD acceleration or not.
/// @param m1 The first matrix.
/// @param m2 The second matrix.
/// @param w The width of the matrices (the number of columns).
/// @param h The height of the matrices (the number of rows).
/// @param res The resulting matrix of the same size.
void matrix_sub(int simd, const float *m1, const float *m2,
                size_t w, size_t h, float *res) NOTNULL(2,3,6);
/// @brief Multiplies two matrices.
/// @param simd Value which indicates whether to use SIMD acceleration or not.
/// @param m1 The first matrix in row-major format.
/// @param m2 The seconds matrix in row-major format.
/// @param w1 The width of the first matrix (the number of columns).
/// @param h1 The height of the first matrix (the number of rows).
/// @param w2 The width of the second matrix (the number of columns).
/// @param h2 The height of the second matrix (the number of rows).
/// @param res The resulting matrix, of size w2 x h1.
/// @pre w1 must be equal to h2.
void matrix_multiply(int simd, const float *m1, const float *m2,
                     size_t w1, size_t h1, size_t w2, size_t h2,
                     float *res) NOTNULL(2,3,8);

/// @brief Multiplies two matrices, the second one being stored transposed.
/// @param simd Value which indicates whether to use SIMD acceleration or not.
/// @param m1 The first matrix in row-major format.
/// @param m2 The seconds matrix in row-major format.
/// @param w1 The width of the first matrix (the number of columns).
/// @param h1 The height of the first matrix (the number of rows).
/// @param w2 The width of the second (transposed) matrix
/// (the number of columns). This value equals to the number of rows
/// in original matrix.
/// @param h2 The height of the second (transposed) matrix
/// (the number of rows). This value equals to the number of columns
/// in original matrix.
/// @param res The resulting matrix, of size h2 x h1.
/// @pre w1 must be equal to w2.
/// @note This function typically performs 10% faster than matrix_multiply().
void matrix_multiply_transposed(int simd, const float *m1, const float *m2,
                                size_t w1, size_t h1, size_t w2, size_t h2,
                                float *res) NOTNULL(2,3,8);

SIMD_API_END

#endif  // INC_SIMD_MATRIX_H_
