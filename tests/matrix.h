/*! @file matrix.h
 *  @brief Matrix unit tests
 *  @author Ernesto Sanches <ernestosanches@gmail.com>
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

#ifndef TESTS_SIMD_MATRIX_H_
#define TESTS_SIMD_MATRIX_H_

#define GTEST_HAS_TR1_TUPLE 1
#include <tuple>
#include <memory>
#include <gtest/gtest.h>

/** @brief Parameterized test for validation of SIMD versions
 *  of matrix functions by comparison of their results with
 *  classic versions.
 *  @details Test parameters in the tuple:
 *  - Tuple containing matrix dimensions:
 *  -- width of the first matrix
 *  -- height of the first matrix
 *  -- width of the second matrix
 *  -- height of the second matrix
 *  - Tuple containing function information
 *  -- function object to be tested
 *  -- Flag that indicates that function accepts transposed second matrix.
 */
class MatrixTest : public ::testing::TestWithParam<
      std::tuple<
        std::tuple<int, int, int, int>,
        std::tuple<
          std::function<void(int, float*, float*, int, int, int, int, float*)>,
          bool>>> {
 protected:
  virtual void SetUp() override;
  void CallFunction(int simd, float* result);
  void CompareResults();

  const std::tuple<int, int, int, int>& dimensions() {
    return std::get<0>(GetParam());
  }
  bool transposed() {
    return std::get<1>(std::get<1>(GetParam()));
  }

  typedef std::shared_ptr<float> FloatPtr;
  FloatPtr m1_;
  FloatPtr m2_;
  FloatPtr res_base_;
  FloatPtr res_simd_;
};

#endif  // TESTS_SIMD_MATRIX_H_
