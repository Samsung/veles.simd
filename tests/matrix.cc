/*! @file matrix.cc
 *  @brief New file description.
 *  @author Markovtsev Vadim <v.markovtsev@samsung.com>
 *  @version 1.0
 *
 *  @section Notes
 *  This code partially conforms to <a href="http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml">Google C++ Style Guide</a>.
 *
 *  @section Copyright
 *  Copyright 2013 Samsung R&D Institute Russia
 */


#include <gtest/gtest.h>
#include <simd/matrix.h>
#include <simd/memory.h>

TEST(Multiply, SIMD) {
  float *m1 = mallocf(128 * 300);
  for (int i = 0; i < 300; i++) {
    for (int j = 0; j < 128; j++) {
      m1[i * 128 + j] = -j % 17 + i % 6;
    }
  }
  float *m2 = mallocf(1000 * 128);
  for (int i = 0; i < 128; i++) {
    for (int j = 0; j < 1000; j++) {
      m2[i * 1000 + j] = j % 15 - i % 5;
    }
  }

  float *res_base = mallocf(1000 * 300);
  float *res_simd = mallocf(1000 * 300);
  matrix_multiply(0, m1, m2, 128, 300, 1000, 128, res_base);
  matrix_multiply(1, m1, m2, 128, 300, 1000, 128, res_simd);

  for (int i = 0; i < 300; i++) {
    for (int j = 0; j < 1000; j++) {
      ASSERT_NEAR(res_base[i * 1000 + j], res_simd[i * 1000 + j], 0.1);
    }
  }

  free(m1);
  free(m2);
  free(res_base);
  free(res_simd);
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

#define CUSTOM_FUNC_PEAK(u1,u2) matrix_multiply(1, m1, m2, 256, 300, 1000, 256,\
                                                res)
#define CUSTOM_FUNC_BASELINE(u1,u2) matrix_multiply(0, m1, m2, 256, 300, 1000,\
                                                    256, res)
#include "tests/benchmark.inc"

#include "tests/google/src/gtest_main.cc"
