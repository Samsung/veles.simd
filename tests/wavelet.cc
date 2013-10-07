/*! @file wavelet.cc
 *  @brief Tests for src/primitives/wavelet.h.
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
#include <chrono>
#include <simd/arithmetic-inl.h>
#include <simd/memory.h>
#include <simd/wavelet.h>
#define WAVELET_INTERNAL_USE
#include "src/daubechies.h"
#include "src/coiflets.h"
#include "src/symlets.h"
#include "tests/make_unique.h"

const int BENCHMARK_LENGTH = 100000;

TEST(Wavelet, wavelet_prepare_array) {
  float array[512];
  int length = sizeof(array) / sizeof(float);  // NOLINT(*)
  for (int i = 0; i < length; i++) {
    array[i] = i;
  }
  auto res = wavelet_prepare_array(8, array, length);

#ifdef __AVX__
  ASSERT_EQ(0, align_complement_f32(res));
  ASSERT_EQ(0, memcmp(array, res, length * sizeof(float)));  // NOLINT(*)
  int checkSize = (length - 8) * sizeof(float);  // NOLINT(*)
  ASSERT_EQ(0, memcmp(array + 2, res + length, checkSize));
  ASSERT_EQ(0, memcmp(array + 4, res + length * 2 - 8, checkSize));
  ASSERT_EQ(0, memcmp(array + 6, res + length * 3 - 16, checkSize));
#else
  ASSERT_EQ(0, memcmp(res, array, sizeof(array)));
#endif
  free(res);

  res = wavelet_prepare_array(4, array, length);

#ifdef __AVX__
  ASSERT_EQ(0, align_complement_f32(res));
  ASSERT_EQ(0, memcmp(array, res, length * sizeof(float)));  // NOLINT(*)
  ASSERT_EQ(0, memcmp(array + 2, res + length, checkSize));
#else
  ASSERT_EQ(0, memcmp(res, array, sizeof(array)));
#endif
  free(res);
}

TEST(Wavelet, wavelet_allocate_destination) {
  auto dest = wavelet_allocate_destination(8, 512);
#ifdef __AVX__
  ASSERT_EQ(0, align_complement_f32(dest));
#endif
  free(dest);
}

#define EPSILON 0.0005f

#define ASSERT_EQF(a, b) ASSERT_NEAR(a, b, EPSILON)

TEST(Wavelet, wavelet_apply_na) {
  float array[32], desthi[16], destlo[16];
  int length = sizeof(array) / sizeof(float);  // NOLINT(*)
  for (int i = 0; i < length; i++) {
    array[i] = i;
  }
  wavelet_apply_na(WAVELET_TYPE_DAUBECHIES, 8, array, length,
                   desthi, destlo);
  int index = 5;
  float vhi = 0.f, vlo = 0.f;
  for (int i = 0; i < 8; i++) {
    vlo += array[index * 2 + i] * kDaubechiesF[3][i];
    vhi += array[index * 2 + i] * kDaubechiesF[3][8 - i - 1] * (i & 1 ? -1 : 1);
  }
  ASSERT_EQF(vlo, destlo[index]);
  ASSERT_EQF(vhi, desthi[index]);
  vhi = vlo = 0.f;
  for (int i = 0; i < 8; i++) {
    float value = i < 2? array[30 + i] :  array[i - 2];
    vlo += value * kDaubechiesF[3][i];
    vhi += value * kDaubechiesF[3][8 - i - 1] * (i & 1 ? -1 : 1);
  }
  ASSERT_EQF(vlo, destlo[15]);
  ASSERT_EQF(vhi, desthi[15]);

  wavelet_apply_na(WAVELET_TYPE_DAUBECHIES, 8, array, 8, desthi, destlo);
}

TEST(Wavelet, stationary_wavelet_apply_na) {
  int length = 32;
  float array[length], desthi[length], destlo[length];
  for (int i = 0; i < length; i++) {
    array[i] = i;
  }
  stationary_wavelet_apply_na(WAVELET_TYPE_DAUBECHIES, 8, array, length,
                              desthi, destlo);
  int index = 5;
  float vhi = 0.f, vlo = 0.f;
  for (int i = 0; i < 8; i++) {
    vlo += array[index + i] * kDaubechiesF[3][i];
    vhi += array[index + i] * kDaubechiesF[3][8 - i - 1] * (i & 1 ? -1 : 1);
  }
  ASSERT_EQF(vlo, destlo[index]);
  ASSERT_EQF(vhi, desthi[index]);
  vhi = vlo = 0.f;
  for (int i = 0; i < 8; i++) {
    float value = i < 2? array[30 + i] :  array[i - 2];
    vlo += value * kDaubechiesF[3][i];
    vhi += value * kDaubechiesF[3][8 - i - 1] * (i & 1 ? -1 : 1);
  }
  ASSERT_EQF(vlo, destlo[30]);
  ASSERT_EQF(vhi, desthi[30]);

  stationary_wavelet_apply_na(WAVELET_TYPE_DAUBECHIES, 8, array, 8,
                              desthi, destlo);
}

class WaveletTest : public ::testing::TestWithParam<
      std::tuple<WaveletType, int>> {
 public:
  typedef std::unique_ptr<float, decltype(&std::free)> FloatPtr;

 protected:
  WaveletTest() : length_(512), array_(nullptr, std::free),
      prep_(nullptr, std::free), desthi_(nullptr, std::free),
      destlo_(nullptr, std::free) {
  }

  virtual void SetUp() override {
    array_ = std::uniquify(mallocf(length_), std::free);
    for (size_t i = 0; i < length_; i++) {
      array_.get()[i] = i;
    }
    prep_ = std::uniquify(wavelet_prepare_array(8, array_.get(), length_),
                          std::free);
    desthi_ = std::uniquify(wavelet_allocate_destination(8, length_),
                            std::free);
    destlo_ = std::uniquify(wavelet_allocate_destination(8, length_),
                            std::free);
  }

  size_t length_;
  FloatPtr array_;
  FloatPtr prep_;
  FloatPtr desthi_;
  FloatPtr destlo_;
};

class StationaryWaveletTest : public WaveletTest {
 protected:
  virtual void SetUp() override {
    array_ = std::uniquify(mallocf(length_), std::free);
    for (size_t i = 0; i < length_; i++) {
      array_.get()[i] = i;
    }
    desthi_ = std::uniquify(mallocf(length_), std::free);
    destlo_ = std::uniquify(mallocf(length_), std::free);
  }
};

TEST_P(WaveletTest, wavelet_apply) {
  wavelet_apply(std::get<0>(GetParam()), std::get<1>(GetParam()),
                prep_.get(), length_, desthi_.get(), destlo_.get());
  float validdesthi[length_ / 2], validdestlo[length_ / 2];
  wavelet_apply_na(std::get<0>(GetParam()), std::get<1>(GetParam()),
                   array_.get(), length_, validdesthi, validdestlo);
  for (size_t i = 0; i < length_ / 2; i++) {
    ASSERT_EQF(validdesthi[i], desthi_.get()[i]) << "i = " << i;
    ASSERT_EQF(validdestlo[i], destlo_.get()[i]) << "i = " << i;
  }
}

TEST_P(StationaryWaveletTest, stationary_wavelet_apply) {
  stationary_wavelet_apply(std::get<0>(GetParam()), std::get<1>(GetParam()),
                           array_.get(), length_, desthi_.get(), destlo_.get());
  float validdesthi[length_], validdestlo[length_];
  stationary_wavelet_apply_na(std::get<0>(GetParam()), std::get<1>(GetParam()),
                              array_.get(), length_, validdesthi, validdestlo);
  for (size_t i = 0; i < length_; i++) {
    ASSERT_EQF(validdesthi[i], desthi_.get()[i]) << "i = " << i;
    ASSERT_EQF(validdestlo[i], destlo_.get()[i]) << "i = " << i;
  }
}

INSTANTIATE_TEST_CASE_P(
    DaubechiesAndSymlets, WaveletTest,
    ::testing::Combine(
        ::testing::Values(WAVELET_TYPE_DAUBECHIES, WAVELET_TYPE_SYMLET),
        ::testing::Values(4, 6, 8, 12, 16)
    ));

INSTANTIATE_TEST_CASE_P(
    Coiflets, WaveletTest,
    ::testing::Combine(
        ::testing::Values(WAVELET_TYPE_COIFLET),
        ::testing::Values(6, 12)
    ));

INSTANTIATE_TEST_CASE_P(
    DaubechiesAndSymlets, StationaryWaveletTest,
    ::testing::Combine(
        ::testing::Values(WAVELET_TYPE_DAUBECHIES, WAVELET_TYPE_SYMLET),
        ::testing::Values(4, 6, 8, 12, 16)
    ));

INSTANTIATE_TEST_CASE_P(
    Coiflets, StationaryWaveletTest,
    ::testing::Combine(
        ::testing::Values(WAVELET_TYPE_COIFLET),
        ::testing::Values(6, 12)
    ));

#ifdef BENCHMARK
#ifdef SIMD
TEST(Wavelet, SIMDSpeedup) {
  float array[512];
  const int length = sizeof(array) / sizeof(float);  // NOLINT(*)
  for (int i = 0; i < length; i++) {
    array[i] = i;
  }

  std::vector<int> orders { 4, 6, 8, 12, 16 };  // NOLINT(*)

  for (int order : orders) {
    auto prep = wavelet_prepare_array(order, array, length);
    auto desthi = wavelet_allocate_destination(order, length);
    auto destlo = wavelet_allocate_destination(order, length);

    auto checkPointStart = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < BENCHMARK_LENGTH; i++) {
      wavelet_apply(WAVELET_TYPE_DAUBECHIES, order, prep, length, desthi, destlo);  // NOLINT(*)
    }

    auto checkPointFinish = std::chrono::high_resolution_clock::now();
    auto delta1 = checkPointFinish - checkPointStart;
    checkPointStart = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < BENCHMARK_LENGTH; i++) {
      wavelet_apply_na(WAVELET_TYPE_DAUBECHIES, order, array, length,
                       desthi, destlo);
    }

    checkPointFinish = std::chrono::high_resolution_clock::now();
    auto delta2 = checkPointFinish - checkPointStart;
    float ratio = (delta1.count() + 0.f) / delta2.count();
    float speedup = (delta2.count() - delta1.count() + 0.f) / delta2.count();
    printf("[order %i] SIMD version took %i%% of original time. "
        "Speedup is %i%%.\n", order,
        static_cast<int>(ratio * 100), static_cast<int>(speedup * 100));

    free(desthi);
    free(destlo);
    free(prep);
  }
}
#endif
#endif

#include "tests/google/src/gtest_main.cc"
