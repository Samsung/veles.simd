/*! @file mathfun.cc
 *  @brief New file description.
 *  @author Ernesto Sanches <ernestosanches@gmail.com>
 *  @version 1.0
 *
 *  @section Notes
 *  This code partially conforms to <a href="http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml">Google C++ Style Guide</a>.
 *
 *  @section Copyright
 *  Copyright 2013 Samsung R&D Institute Russia
 */

#define GTEST_HAS_TR1_TUPLE 1
#include <tuple>
#include <cmath>
#include <memory>
#include <functional>
#include <cstdio>
#include <gtest/gtest.h>
#include <simd/mathfun.h>
#include <simd/memory.h>

typedef std::function<void(int, const float *, size_t, float *)> TestedFunc;
typedef std::function<float(float)> ReferenceFunc;

class MathTest : public ::testing::TestWithParam<
  std::tuple<bool, size_t, std::tuple<TestedFunc, ReferenceFunc>>> {
 protected:
  virtual void SetUp() override {
    std::tuple<TestedFunc, ReferenceFunc> funcs;
    std::tie(is_simd_, length_, funcs) = GetParam();
    std::tie(tested_func_, reference_func_) = funcs;
  }

  bool is_simd_;
  size_t length_;
  TestedFunc tested_func_;
  ReferenceFunc reference_func_;
};

TEST_P(MathTest, MathTestCase) {
  std::unique_ptr<float[], void(*)(void*)> data(mallocf(length_), std::free);
  for (size_t i = 0; i < length_; ++i) {
    data[i] = (float(i % 19) - float(i % 6)) * 1.1f;
  }
  std::unique_ptr<float[], void(*)(void*)> test_result(mallocf(length_),
                                                       std::free);
  tested_func_(is_simd_, data.get(), length_, test_result.get());
  for (size_t i = 0; i < length_; ++i) {
    float reference_result = reference_func_(data[i]);
    if(std::isfinite(reference_result) && std::isfinite(test_result[i])) {
      ASSERT_FLOAT_EQ(reference_result, test_result[i]) << "i = " << i;
    }
  }
}

INSTANTIATE_TEST_CASE_P(
    Math, MathTest,
    ::testing::Combine(
        ::testing::Bool(),
        ::testing::Values(1, 3, 64, 199),
        ::testing::Values(
            std::make_tuple(sin_psv, sin),
            std::make_tuple(cos_psv, cos),
            std::make_tuple(log_psv, log),
            std::make_tuple(exp_psv, exp))));

#include "tests/google/src/gtest_main.cc"

