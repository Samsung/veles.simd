/*! @file detect_peaks.cc
 *  @brief Tests for detect_peaks().
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


#include <simd/detect_peaks.h>
#include <cmath>
#include <gtest/gtest.h>

class DetectPeaksTest : public ::testing::TestWithParam<bool> {
 protected:
  bool is_simd() {
    return GetParam();
  }
};

TEST_P(DetectPeaksTest, sin) {
  size_t length = 4000;
  float array[length];
  for (size_t i = 0; i < length; i++) {
    array[i] = sinf(i * M_PI / 100);
  }
  ExtremumPoint *points;
  size_t points_count;
  detect_peaks(is_simd(), array, length, kExtremumTypeMaximum, &points,
               &points_count);
  ASSERT_EQ(20U, points_count);
  for (int i = 0; i < 20; i++) {
    ASSERT_EQ(i * 200 + 50, points[i].position) << i;
    ASSERT_FLOAT_EQ(1.f, points[i].value) << i;
  }
  free(points);
  detect_peaks(is_simd(), array, length, kExtremumTypeMinimum, &points,
               &points_count);
  ASSERT_EQ(20U, points_count);
  for (int i = 0; i < 20; i++) {
    ASSERT_EQ(i * 200 + 150, points[i].position) << i;
    ASSERT_FLOAT_EQ(-1.f, points[i].value) << i;
  }
  free(points);
  detect_peaks(is_simd(), array, length, kExtremumTypeBoth, &points,
               &points_count);
  ASSERT_EQ(40U, points_count);
  for (int i = 0; i < 40; i++) {
    ASSERT_EQ((i / 2) * 200 + 50 + 100 * (i % 2), points[i].position) << i;
    ASSERT_FLOAT_EQ(1.f - (i % 2) * 2, points[i].value) << i;
  }
  free(points);
}

TEST_P(DetectPeaksTest, nasty_peaks) {
  size_t length = 101;
  float array[length];
  for (size_t i = 0; i < length; i++) {
    array[i] = 0;
  }
  array[7] = 1;
  array[16] = 1;
  array[97] = 1;
  array[99] = 1;
  ExtremumPoint *points;
  size_t points_count;
  detect_peaks(is_simd(), array, length, kExtremumTypeMaximum, &points,
               &points_count);
  ASSERT_EQ(4U, points_count);
  ASSERT_EQ(7, points[0].position) << "0";
  ASSERT_EQ(16, points[1].position) << "1";
  ASSERT_EQ(97, points[2].position) << "2";
  ASSERT_EQ(99, points[3].position) << "3";
  for (int i = 0; i < 4; i++) {
    ASSERT_FLOAT_EQ(1.f, points[i].value) << i;
  }
  free(points);
}

INSTANTIATE_TEST_CASE_P(DetectPeaksTests, DetectPeaksTest, ::testing::Bool());

#include "tests/google/src/gtest_main.cc"
