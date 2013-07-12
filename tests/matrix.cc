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

TEST(Multiply, SIMDUglyLength) {
  float *m1 = mallocf(125 * 299);
  for (int i = 0; i < 299; i++) {
    for (int j = 0; j < 125; j++) {
      m1[i * 125 + j] = -j % 17 + i % 6;
    }
  }
  float *m2 = mallocf(999 * 125);
  for (int i = 0; i < 125; i++) {
    for (int j = 0; j < 999; j++) {
      m2[i * 999 + j] = j % 15 - i % 5;
    }
  }

  float *res_base = mallocf(999 * 299);
  float *res_simd = mallocf(999 * 299);
  matrix_multiply(0, m1, m2, 125, 299, 999, 125, res_base);
  matrix_multiply(1, m1, m2, 125, 299, 999, 125, res_simd);

  for (int i = 0; i < 299; i++) {
    for (int j = 0; j < 999; j++) {
      ASSERT_NEAR(res_base[i * 999 + j], res_simd[i * 999 + j], 0.1);
    }
  }

  free(m1);
  free(m2);
  free(res_base);
  free(res_simd);
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

TEST(MultiplyTransposed, SIMD) {
  float *m1 = mallocf(128 * 300);
  for (int i = 0; i < 300; i++) {
    for (int j = 0; j < 128; j++) {
      m1[i * 128 + j] = -j % 17 + i % 6;
    }
  }
  float *m2 = mallocf(128 * 1000);
  for (int i = 0; i < 1000; i++) {
    for (int j = 0; j < 128; j++) {
      m2[i * 128 + j] = j % 15 - i % 5;
    }
  }

  float *res_base = mallocf(1000 * 300);
  float *res_simd = mallocf(1000 * 300);
  matrix_multiply_transposed(0, m1, m2, 128, 300, 128, 1000, res_base);
  matrix_multiply_transposed(1, m1, m2, 128, 300, 128, 1000, res_simd);

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

TEST(MultiplyTransposed, SIMDUglyLength) {
  float *m1 = mallocf(125 * 299);
  for (int i = 0; i < 299; i++) {
    for (int j = 0; j < 125; j++) {
      m1[i * 125 + j] = -j % 17 + i % 6;
    }
  }
  float *m2 = mallocf(125 * 999);
  for (int i = 0; i < 999; i++) {
    for (int j = 0; j < 125; j++) {
      m2[i * 125 + j] = j % 15 - i % 5;
    }
  }

  float *res_base = mallocf(999 * 299);
  float *res_simd = mallocf(999 * 299);
  matrix_multiply_transposed(0, m1, m2, 125, 299, 125, 999, res_base);
  matrix_multiply_transposed(1, m1, m2, 125, 299, 125, 999, res_simd);

  for (int i = 0; i < 299; i++) {
    for (int j = 0; j < 999; j++) {
      ASSERT_NEAR(res_base[i * 999 + j], res_simd[i * 999 + j], 0.1);
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
