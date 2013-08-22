/*! @file matrix.cc
 *  @brief Matrix unit tests
 *  @author Markovtsev Vadim <v.markovtsev@samsung.com>
 *  @version 1.0
 *
 *  @section Notes
 *  This code partially conforms to <a href="http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml">Google C++ Style Guide</a>.
 *
 *  @section Copyright
 *  Copyright 2013 Samsung R&D Institute Russia
 */

#include <simd/memory.h>
#include <simd/matrix.h>
#include "tests/matrix.h"

void MatrixTest::SetUp() {
  int w1, h1, w2, h2;
  std::tie(w1, h1, w2, h2) = dimensions();
  res_base_ = FloatPtr(mallocf(w2 * h1), std::free);
  res_simd_ = FloatPtr(mallocf(w2 * h1), std::free);
  if(transposed()) {
    std::swap(w2, h2);
  }
  m1_ = FloatPtr(mallocf(w1 * h1), std::free);
  m2_ = FloatPtr(mallocf(w2 * h2), std::free);
  for (int i = 0; i < h1; i++) {
    for (int j = 0; j < w1; j++) {
      m1_.get()[i * w1 + j] = - j % 17 + i % 6;
    }
  }
  for (int i = 0; i < h2; i++) {
    for (int j = 0; j < w2; j++) {
      m2_.get()[i * w2 + j] = j % 15 - i % 5;
    }
  }
}

void MatrixTest::CallFunction(int simd, float* result) {
  int w1, h1, w2, h2;
  std::tie(w1, h1, w2, h2) = dimensions();
  if (transposed()) {
    std::swap(w2, h2);
  }
  std::get<0>(std::get<1>(GetParam()))(
      simd, m1_.get(), m2_.get(), w1, h1, w2, h2, result);
}

void MatrixTest::CompareResults() {
  int h1, w2;
  std::tie(std::ignore, h1, w2, std::ignore) = dimensions();
  for (int i = 0; i < h1; i++) {
    for (int j = 0; j < w2; j++) {
      ASSERT_NEAR(res_base_.get()[i * w2 + j],
                  res_simd_.get()[i * w2 + j], 0.1) << i << " " << j;
    }
  }
}

/** @brief Wrapper to instantiate parameterized test with matrix_add function
 *  using common signature
 */
void matrix_add_wrapper(int simd, const float *m1, const float *m2,
                        size_t w, size_t h, size_t, size_t,
                        float *res) {
  matrix_add(simd, m1, m2, w, h, res);
}


TEST_P(MatrixTest, SIMD) {
  CallFunction(0, res_base_.get());
  CallFunction(1, res_simd_.get());
  CompareResults();
}

INSTANTIATE_TEST_CASE_P(
    Common, MatrixTest,
    ::testing::Combine(
        ::testing::Values(
            std::make_tuple(1, 1, 1, 1),
            std::make_tuple(3, 3, 3, 3),
            std::make_tuple(99, 99, 99, 99)
        ),
        ::testing::Values(
            std::make_tuple(matrix_add_wrapper, false),
            std::make_tuple(matrix_multiply, false),
            std::make_tuple(matrix_multiply_transposed, true)
        )
    ));

INSTANTIATE_TEST_CASE_P(
    Add, MatrixTest,
    ::testing::Values<MatrixTest::ParamType>(
      std::make_tuple(
          std::make_tuple(125, 299, 125, 299),
          std::make_tuple(matrix_add_wrapper, false)
      )
    ));

INSTANTIATE_TEST_CASE_P(
    Multiply, MatrixTest,
    ::testing::Combine(
        ::testing::Values(
            std::make_tuple(128, 300, 1000, 128),
            std::make_tuple(125, 299, 999, 125)
        ),
        ::testing::Values(
            std::make_tuple(matrix_multiply, false),
            std::make_tuple(matrix_multiply_transposed, true)
        )
    ));


TEST(Add, Validate) {
  float m1[6] = { 1, 2, 3,
                 -2, 0, 4 };
  float m2[6] = { 0, 1, 3,
                  5, -1, 2 };
  float res[6];
  float res_valid[6] = { 1, 3, 6,
                         3, -1, 6 };
  matrix_add(0, m1, m2, 3, 2, res);
  for (int i = 0; i < 6; i++) {
    ASSERT_NEAR(res[i], res_valid[i], 0.01);
  }
}

TEST(Multiply, Validate) {
  float m1[6] = { 1, 2, 3,
                 -2, 0, 4 };
  float m2[12] = { 0, 1, 3, -2,
                   5, -1, 2, 4,
                  -3, 0, -4, 2 };
  float res[8];
  float res_valid[8] = { 1, -1, -5, 12,
                       -12, -2,-22, 12 };
  matrix_multiply(0, m1, m2, 3, 2, 4, 3, res);
  for (int i = 0; i < 8; i++) {
    ASSERT_NEAR(res[i], res_valid[i], 0.01);
  }
}

TEST(MultiplyTransposed, Validate) {
  float m1[6] = { 1, 2, 3,
                 -2, 0, 4 };
  float m2[12] = { 0, 5, -3,
                   1,-1,  0,
                   3, 2, -4,
                  -2, 4,  2 };
  float res[8];
  float res_valid[8] = { 1, -1, -5, 12,
                       -12, -2,-22, 12 };
  matrix_multiply_transposed(0, m1, m2, 3, 2, 3, 4, res);
  for (int i = 0; i < 8; i++) {
    ASSERT_NEAR(res[i], res_valid[i], 0.01);
  }
}


#define TEST_NAME matrix_multiply
#define ITER_COUNT 10
#define CUSTOM_CODE_PRE  \
  float *m1 = mallocf(256 * 300);\
  for (int i = 0; i < 300; i++) {\
    for (int j = 0; j < 256; j++) {\
      m1[i * 256 + j] = -j % 17 + i % 6;\
    }\
  }\
  float *m2 = mallocf(1000 * 256);\
  for (int i = 0; i < 256; i++) {\
    for (int j = 0; j < 1000; j++) {\
      m2[i * 1000 + j] = j % 15 - i % 5;\
    }\
  }\
\
  float *res = mallocf(1000 * 300);\

#define CUSTOM_CODE_POST \
  free(m1); free(m2); free(res);

#define CUSTOM_FUNC_PEAK(u1,u2) matrix_multiply(1, m1, m2, 256, 300, 1000, 256,\
                                                res)
#define CUSTOM_FUNC_BASELINE(u1,u2) matrix_multiply(0, m1, m2, 256, 300, 1000,\
                                                    256, res)
#include "tests/benchmark.inc"

#undef TEST_NAME
#define TEST_NAME matrix_multiply_transposed
#undef CUSTOM_CODE_PRE
#define CUSTOM_CODE_PRE  \
  float *m1 = mallocf(256 * 300);\
  for (int i = 0; i < 300; i++) {\
    for (int j = 0; j < 256; j++) {\
      m1[i * 256 + j] = -j % 17 + i % 6;\
    }\
  }\
  float *m2 = mallocf(256 * 1000);\
  for (int i = 0; i < 1000; i++) {\
    for (int j = 0; j < 256; j++) {\
      m2[i * 256 + j] = j % 15 - i % 5;\
    }\
  }\
\
  float *res = mallocf(1000 * 300);\

#undef CUSTOM_FUNC_PEAK
#define CUSTOM_FUNC_PEAK(u1,u2) matrix_multiply_transposed(1, m1, m2, 256, 300,\
                                                           256, 1000, res)
#undef CUSTOM_FUNC_BASELINE
#define CUSTOM_FUNC_BASELINE(u1,u2) matrix_multiply_transposed(0, m1, m2, 256, \
                                                               300, 256, 1000, \
                                                               res)
#include "tests/benchmark.inc"

#undef TEST_NAME
#define TEST_NAME matrix_multiply_straight_vs_transposed
#undef ITER_COUNT
#define ITER_COUNT 40
#undef CUSTOM_CODE_PRE
#define CUSTOM_CODE_PRE  \
  float *m1 = mallocf(256 * 300);\
  for (int i = 0; i < 300; i++) {\
    for (int j = 0; j < 256; j++) {\
      m1[i * 256 + j] = -j % 17 + i % 6;\
    }\
  }\
  float *m2 = mallocf(256 * 1000);\
  for (int i = 0; i < 1000; i++) {\
    for (int j = 0; j < 256; j++) {\
      m2[i * 256 + j] = j % 15 - i % 5;\
    }\
  }\
\
  float *res = mallocf(1000 * 300);\

#undef CUSTOM_FUNC_PEAK
#define CUSTOM_FUNC_PEAK(u1,u2) matrix_multiply_transposed(1, m1, m2, 256, 300,\
                                                           256, 1000, res)
#undef CUSTOM_FUNC_BASELINE
#define CUSTOM_FUNC_BASELINE(u1,u2) matrix_multiply(1, m1, m2, 256, \
                                                    300, 1000, 256, res)
#include "tests/benchmark.inc"

#include "tests/google/src/gtest_main.cc"
