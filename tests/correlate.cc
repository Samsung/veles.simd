/*! @file correlate.cc
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
#ifndef NO_FFTF
#include <math.h>
#include <simd/correlate.h>
#include <simd/memory.h>
#include <simd/arithmetic-inl.h>
#include <fftf/api.h>

void cross_correlate_reference(const float *__restrict x, size_t xLength,
                         const float *__restrict h, size_t hLength,
                         float *__restrict result) {
  cross_correlate_simd(0, x, xLength, h, hLength, result);
}

void DebugPrintConvolution(const char* name, const float* vec) {
  printf("%s\t", name);
  for (int i = 0; i < 40; i++) {
    printf("%f  ", vec[i]);
  }
  printf("\n");
}

TEST(correlate, ross_correlate_reference) {
  float x[] = { 1, 2, 3, 4, 5, 6, 7, 8 };
  float y[] = { 10, 9, 8, 7 };
  float z[11];
  cross_correlate_reference(x, sizeof(x)/sizeof(float),
                            y, sizeof(y)/sizeof(float),
                            z);
  ASSERT_NEAR(z[0], 7, 0.0001f);
  ASSERT_NEAR(z[1], 22, 0.0001f);
  ASSERT_NEAR(z[2], 46, 0.0001f);
  ASSERT_NEAR(z[3], 80, 0.0001f);
  ASSERT_NEAR(z[4], 114, 0.0001f);
  ASSERT_NEAR(z[5], 148, 0.0001f);
  ASSERT_NEAR(z[6], 182, 0.0001f);
  ASSERT_NEAR(z[7], 216, 0.0001f);
  ASSERT_NEAR(z[8], 187, 0.0001f);
  ASSERT_NEAR(z[9], 142, 0.0001f);
  ASSERT_NEAR(z[10], 80, 0.0001f);
}

TEST(correlate, cross_correlate_fft) {
  const int xlen = 1020;
  const int hlen = 50;

  float x[xlen];
  for (int i = 0; i < xlen; i++) {
    x[i] = sinf(i) * 100;
  }
  float h[hlen];
  for (int i = 0; i < hlen; i++) {
    h[i] = i / (hlen - 1.0f);
  }

  float verif[xlen + hlen - 1];
  cross_correlate_reference(x, xlen, h, hlen, verif);
  DebugPrintConvolution("REFERENCE", verif);

  float res[xlen + hlen - 1];
  auto handle = cross_correlate_fft_initialize(xlen, hlen);
  cross_correlate_fft(handle, x, h, res);
  cross_correlate_fft_finalize(handle);
  DebugPrintConvolution("FFT\t", res);

  int firstDifferenceIndex = -1;
  for (int i = 0; i < xlen + hlen - 1; i++) {
    float delta = res[i] - verif[i];
    if (delta * delta > 1E-6 && firstDifferenceIndex == -1) {
      firstDifferenceIndex = i;
    }
  }
  ASSERT_EQ(-1, firstDifferenceIndex);
}

TEST(correlate, cross_correlate_overlap_save) {
  const int xlen = 1021;
  const int hlen = 50;

  float x[xlen];
  for (int i = 0; i < xlen; i++) {
    x[i] = sinf(i) * 100;
  }
  float h[hlen];
  for (int i = 0; i < hlen; i++) {
    h[i] = i / (hlen- 1.0f);
  }

  float verif[xlen + hlen - 1];
  cross_correlate_reference(x, xlen, h, hlen, verif);
  DebugPrintConvolution("REFERENCE", verif);

  float res[xlen + hlen - 1];
  auto handle = cross_correlate_overlap_save_initialize(xlen, hlen);
  cross_correlate_overlap_save(handle, x, h, res);
  cross_correlate_overlap_save_finalize(handle);
  DebugPrintConvolution("OVERLAP-SAVE", res);

  int firstDifferenceIndex = -1;
  for (int i = 0; i < xlen + hlen - 1; i++) {
    float delta = res[i] - verif[i];
    if (delta * delta > 1E-6 && firstDifferenceIndex == -1) {
      firstDifferenceIndex = i;
    }
  }
  ASSERT_EQ(-1, firstDifferenceIndex);
}

TEST(correlate, cross_correlate_simd) {
  const int xlen = 1024;
  const int hlen = 50;

  float x[xlen];
  for (int i = 0; i < xlen; i++) {
    x[i] = sinf(i) * 100;
  }
  float h[hlen];
  for (int i = 0; i < hlen; i++) {
    h[i] = i / (hlen - 1.0f);
  }

  float verif[xlen + hlen - 1];
  cross_correlate_reference(x, xlen, h, hlen, verif);

  float res[xlen + hlen - 1];
  cross_correlate_simd(true, x, xlen, h, hlen, res);

  int firstDifferenceIndex = -1;
  for (int i = 0; i < xlen + hlen - 1; i++) {
    float delta = res[i] - verif[i];
    if (delta * delta > 1E-6 && firstDifferenceIndex == -1) {
      firstDifferenceIndex = i;
    }
  }
  ASSERT_EQ(-1, firstDifferenceIndex);
}

#endif

#include "tests/google/src/gtest_main.cc"
