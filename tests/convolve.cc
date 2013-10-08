/*! @file convolve.cc
 *  @brief Tests for src/primitives/convolve.cc.
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
#include <simd/convolve.h>
#include <simd/memory.h>
#include <simd/arithmetic-inl.h>
#include <fftf/api.h>

void convolve_reference(const float *__restrict x, size_t xLength,
                        const float *__restrict h, size_t hLength,
                        float *__restrict result) {
  convolve_simd(false, x, xLength, h, hLength, result);
}

void DebugPrintConvolution(const char* name, const float* vec) {
  printf("%s\t", name);
  for (int i = 0; i < 40; i++) {
    printf("%f  ", vec[i]);
  }
  printf("\n");
}

TEST(convolve, convolve_reference) {
  float x[] = { 1, 2, 3, 4, 5, 6, 7, 8 };
  float y[] = { 10, 9, 8, 7 };
  float z[11];
  convolve_reference(x, sizeof(x) / sizeof(float),
                     y, sizeof(y) / sizeof(float),
                     z);
  ASSERT_NEAR(z[0], 10, 0.0001f);
  ASSERT_NEAR(z[1], 29, 0.0001f);
  ASSERT_NEAR(z[2], 56, 0.0001f);
  ASSERT_NEAR(z[3], 90, 0.0001f);
  ASSERT_NEAR(z[4], 124, 0.0001f);
  ASSERT_NEAR(z[5], 158, 0.0001f);
  ASSERT_NEAR(z[6], 192, 0.0001f);
  ASSERT_NEAR(z[7], 226, 0.0001f);
  ASSERT_NEAR(z[8], 170, 0.0001f);
  ASSERT_NEAR(z[9], 113, 0.0001f);
  ASSERT_NEAR(z[10], 56, 0.0001f);
}

TEST(convolve, convolve_fft) {
  const int xlen = 1020;
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
  convolve_reference(x, xlen, h, hlen, verif);
  DebugPrintConvolution("REFERENCE", verif);

  float res[xlen + hlen - 1];
  auto handle = convolve_fft_initialize(xlen, hlen);
  convolve_fft(handle, x, h, res);
  convolve_fft_finalize(handle);
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

TEST(convolve, convolve_overlap_save) {
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
  convolve_reference(x, xlen, h, hlen, verif);
  DebugPrintConvolution("REFERENCE", verif);

  float res[xlen + hlen - 1];
  auto handle = convolve_overlap_save_initialize(xlen, hlen);
  convolve_overlap_save(handle, x, h, res);
  convolve_overlap_save_finalize(handle);
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

TEST(convolve, convolve_simd) {
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
  convolve_reference(x, xlen, h, hlen, verif);

  float res[xlen + hlen - 1];
  convolve_simd(true, x, xlen, h, hlen, res);

  int firstDifferenceIndex = -1;
  for (int i = 0; i < xlen + hlen - 1; i++) {
    float delta = res[i] - verif[i];
    if (delta * delta > 1E-6 && firstDifferenceIndex == -1) {
      firstDifferenceIndex = i;
    }
  }
  ASSERT_EQ(-1, firstDifferenceIndex);
}

float BenchmarkH[512] = { 1.f };
float BenchmarkResult[10000];

#define TEST_NAME convolve_simd_50
#ifndef __arm__
#define ITER_COUNT 50000
#else
#define ITER_COUNT 1000
#endif
#define BENCH_FUNC convolve_simd
#define NO_OUTPUT
#define EXTRA_PARAM BenchmarkH, sizeof(BenchmarkH) / sizeof(BenchmarkH[0]), \
  BenchmarkResult
#include "tests/benchmark.inc"

#undef EXTRA_PARAM
#undef ITER_COUNT
#undef TEST_NAME
#define TEST_NAME convolve_simd_500
#ifndef __arm__
#define ITER_COUNT 10000
#else
#define ITER_COUNT 1000
#endif
#define EXTRA_PARAM BenchmarkH, sizeof(BenchmarkH) / sizeof(BenchmarkH[0]), \
  BenchmarkResult
#include "tests/benchmark.inc"

ConvolutionFFTHandle fftHandle;

#undef OPT_STRING
#define OPT_STRING "FFT"
#undef LENGTH
#define LENGTH 512
#define CUSTOM_FUNC_BASELINE(x, xLength) convolve_simd(\
    true, x, xLength, BenchmarkH, 512, BenchmarkResult)
#define CUSTOM_FUNC_PEAK(x, xLength) convolve_fft(\
    fftHandle, x, BenchmarkH, BenchmarkResult)
#undef ITER_COUNT
#undef TEST_NAME
#define TEST_NAME convolve_simd_vs_convolve_fft_512_512
#ifndef __arm__
#define ITER_COUNT 25000
#else
#define ITER_COUNT 2000
#endif
#define CUSTOM_CODE_PRE { fftHandle = convolve_fft_initialize(512, 512); }
#define CUSTOM_CODE_POST { convolve_fft_finalize(fftHandle); }
#include "tests/benchmark.inc"

#ifndef __arm__
#undef LENGTH
#define LENGTH 256
#undef CUSTOM_FUNC_BASELINE
#define CUSTOM_FUNC_BASELINE(x, xLength) convolve_simd(\
    true, x, xLength, BenchmarkH, 256, BenchmarkResult)
#undef CUSTOM_FUNC_PEAK
#define CUSTOM_FUNC_PEAK(x, xLength) convolve_fft(\
    fftHandle, x, BenchmarkH, BenchmarkResult)
#undef ITER_COUNT
#undef TEST_NAME
#define TEST_NAME convolve_simd_vs_convolve_fft_256_256
#define ITER_COUNT 60000
#undef CUSTOM_CODE_PRE
#define CUSTOM_CODE_PRE { fftHandle = convolve_fft_initialize(256, 256); }
#include "tests/benchmark.inc"

#undef LENGTH
#define LENGTH 128
#undef CUSTOM_FUNC_BASELINE
#define CUSTOM_FUNC_BASELINE(x, xLength) convolve_simd(\
    true, x, xLength, BenchmarkH, 128, BenchmarkResult)
#undef CUSTOM_FUNC_PEAK
#define CUSTOM_FUNC_PEAK(x, xLength) convolve_fft(\
    fftHandle, x, BenchmarkH, BenchmarkResult)
#undef ITER_COUNT
#undef TEST_NAME
#define TEST_NAME convolve_simd_vs_convolve_fft_128_128
#define ITER_COUNT 80000
#undef CUSTOM_CODE_PRE
#define CUSTOM_CODE_PRE { fftHandle = convolve_fft_initialize(128, 128); }
#include "tests/benchmark.inc"

#else

#undef LENGTH
#define LENGTH 64
#undef CUSTOM_FUNC_BASELINE
#define CUSTOM_FUNC_BASELINE(x, xLength) convolve_simd(\
    true, x, xLength, BenchmarkH, 64, BenchmarkResult)
#undef CUSTOM_FUNC_PEAK
#define CUSTOM_FUNC_PEAK(x, xLength) convolve_fft(\
    fftHandle, x, BenchmarkH, BenchmarkResult)
#undef ITER_COUNT
#undef TEST_NAME
#define TEST_NAME convolve_simd_vs_convolve_fft_64_64
#define ITER_COUNT 80000
#undef CUSTOM_CODE_PRE
#define CUSTOM_CODE_PRE { fftHandle = convolve_fft_initialize(64, 64); }
#include "tests/benchmark.inc"

#undef LENGTH
#define LENGTH 50
#undef CUSTOM_FUNC_BASELINE
#define CUSTOM_FUNC_BASELINE(x, xLength) convolve_simd(\
    true, x, xLength, BenchmarkH, 50, BenchmarkResult)
#undef CUSTOM_FUNC_PEAK
#define CUSTOM_FUNC_PEAK(x, xLength) convolve_fft(\
    fftHandle, x, BenchmarkH, BenchmarkResult)
#undef ITER_COUNT
#undef TEST_NAME
#define TEST_NAME convolve_simd_vs_convolve_fft_50_50
#define ITER_COUNT 80000
#undef CUSTOM_CODE_PRE
#define CUSTOM_CODE_PRE { fftHandle = convolve_fft_initialize(50, 50); }
#include "tests/benchmark.inc"

#undef LENGTH
#define LENGTH 32
#undef CUSTOM_FUNC_BASELINE
#define CUSTOM_FUNC_BASELINE(x, xLength) convolve_simd(\
    true, x, xLength, BenchmarkH, 32, BenchmarkResult)
#undef CUSTOM_FUNC_PEAK
#define CUSTOM_FUNC_PEAK(x, xLength) convolve_fft(\
    fftHandle, x, BenchmarkH, BenchmarkResult)
#undef ITER_COUNT
#undef TEST_NAME
#define TEST_NAME convolve_simd_vs_convolve_fft_32_32
#define ITER_COUNT 80000
#undef CUSTOM_CODE_PRE
#define CUSTOM_CODE_PRE { fftHandle = convolve_fft_initialize(32, 32); }
#include "tests/benchmark.inc"

#endif


#ifdef __AVX__
#undef LENGTH
#define LENGTH 350
#undef CUSTOM_FUNC_BASELINE
#define CUSTOM_FUNC_BASELINE(x, xLength) convolve_simd(\
    true, x, xLength, BenchmarkH, 350, BenchmarkResult)
#undef CUSTOM_FUNC_PEAK
#define CUSTOM_FUNC_PEAK(x, xLength) convolve_fft(\
    fftHandle, x, BenchmarkH, BenchmarkResult)
#undef ITER_COUNT
#undef TEST_NAME
#define TEST_NAME convolve_simd_vs_convolve_fft_350_350
#define ITER_COUNT 40000
#undef CUSTOM_CODE_PRE
#define CUSTOM_CODE_PRE { fftHandle = convolve_fft_initialize(350, 350); }
#include "tests/benchmark.inc"
#endif

ConvolutionOverlapSaveHandle osHandle;

#undef OPT_STRING
#define OPT_STRING "Overlap-Save"
#undef LENGTH
#define LENGTH 1000
#undef CUSTOM_FUNC_BASELINE
#define CUSTOM_FUNC_BASELINE(x, xLength) convolve_fft(\
    fftHandle, x, BenchmarkH, BenchmarkResult)
#undef CUSTOM_FUNC_PEAK
#define CUSTOM_FUNC_PEAK(x, xLength) convolve_overlap_save(\
    osHandle, x, BenchmarkH, BenchmarkResult)
#undef ITER_COUNT
#undef TEST_NAME
#define TEST_NAME convolve_fft_vs_convolve_overlap_save_1000_50
#define ITER_COUNT 30000
#undef CUSTOM_CODE_PRE
#define CUSTOM_CODE_PRE { osHandle = convolve_overlap_save_initialize(1000, 50); \
                          fftHandle = convolve_fft_initialize(1000, 50); }
#undef CUSTOM_CODE_POST
#define CUSTOM_CODE_POST { convolve_overlap_save_finalize(osHandle); \
                           convolve_fft_finalize(fftHandle); }
#include "tests/benchmark.inc"

#undef LENGTH
#define LENGTH 2000
#undef CUSTOM_FUNC_BASELINE
#define CUSTOM_FUNC_BASELINE(x, xLength) convolve_fft(\
    fftHandle, x, BenchmarkH, BenchmarkResult)
#undef CUSTOM_FUNC_PEAK
#define CUSTOM_FUNC_PEAK(x, xLength) convolve_overlap_save(\
    osHandle, x, BenchmarkH, BenchmarkResult)
#undef ITER_COUNT
#undef TEST_NAME
#define TEST_NAME convolve_fft_vs_convolve_overlap_save_2000_950
#define ITER_COUNT 10000
#undef CUSTOM_CODE_PRE
#define CUSTOM_CODE_PRE { osHandle = convolve_overlap_save_initialize(2000, 950); \
                          fftHandle = convolve_fft_initialize(2000, 950); }
#undef CUSTOM_CODE_POST
#define CUSTOM_CODE_POST { convolve_overlap_save_finalize(osHandle); \
                           convolve_fft_finalize(fftHandle); }
#include "tests/benchmark.inc"

#undef LENGTH
#define LENGTH 1000
#undef CUSTOM_FUNC_BASELINE
#define CUSTOM_FUNC_BASELINE(x, xLength) convolve_simd(\
    true, x, xLength, BenchmarkH, 50, BenchmarkResult)
#undef CUSTOM_FUNC_PEAK
#define CUSTOM_FUNC_PEAK(x, xLength) convolve_overlap_save(\
    osHandle, x, BenchmarkH, BenchmarkResult)
#undef ITER_COUNT
#undef TEST_NAME
#define TEST_NAME convolve_simd_vs_convolve_overlap_save_1000_50
#define ITER_COUNT 20000
#undef CUSTOM_CODE_PRE
#define CUSTOM_CODE_PRE { osHandle = convolve_overlap_save_initialize(1000, 50); }
#undef CUSTOM_CODE_POST
#define CUSTOM_CODE_POST { convolve_overlap_save_finalize(osHandle); }
#include "tests/benchmark.inc"

#undef LENGTH
#define LENGTH 200
#undef CUSTOM_FUNC_BASELINE
#define CUSTOM_FUNC_BASELINE(x, xLength) convolve_simd(\
    true, x, xLength, BenchmarkH, 50, BenchmarkResult)
#undef CUSTOM_FUNC_PEAK
#define CUSTOM_FUNC_PEAK(x, xLength) convolve_overlap_save(\
    osHandle, x, BenchmarkH, BenchmarkResult)
#undef ITER_COUNT
#undef TEST_NAME
#define TEST_NAME convolve_simd_vs_convolve_overlap_save_200_50
#define ITER_COUNT 80000
#undef CUSTOM_CODE_PRE
#define CUSTOM_CODE_PRE { osHandle = convolve_overlap_save_initialize(200, 50); }
#undef CUSTOM_CODE_POST
#define CUSTOM_CODE_POST { convolve_overlap_save_finalize(osHandle); }
#include "tests/benchmark.inc"

#endif  // #ifndef NO_FFTF

#include "tests/google/src/gtest_main.cc"
