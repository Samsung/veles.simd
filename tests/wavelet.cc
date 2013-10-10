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
  const float valid_destlo[] = { 1.42184071797210, 4.25026784271829,
      7.07869496746448, 9.90712209221067, 12.7355492169569, 15.5639763417030,
      18.3924034664492, 21.2208305911954, 24.0492577159416, 26.8776848406878,
      29.7061119654340, 32.5345390901802, 35.3629662149264, 37.4782538234490,
      45.3048707044478, 28.8405938767906 };

  const float valid_desthi[] = { -9.91075277401166e-13, -9.90367510222967e-13,
      -9.90194037875369e-13, -9.91873250200115e-13, -9.91456916565880e-13,
      -9.91096094082877e-13, -9.90263426814408e-13, -9.89069937062936e-13,
      -9.91706716746421e-13, -9.92234072683118e-13, -9.92872450922278e-13,
      -9.91484672141496e-13, -9.88431558823777e-13, -15.5030002317990,
      5.58066496329142, -1.39137323046436 };

  for (int i = 0; i < length / 2; i++) {
     ASSERT_NEAR(desthi[i], valid_desthi[i], 1e-5) << i;
     ASSERT_FLOAT_EQ(destlo[i], valid_destlo[i]) << i;
   }
}

TEST(Wavelet, stationary_wavelet_apply_na) {
  int length = 32;
  float array[length], desthi1[length], desthi2[length], destlo[length];
  for (int i = 0; i < length; i++) {
    array[i] = i;
  }
  const float valid_desthi1[] = { -9.91075277401166e-13,
      -9.90107301701571e-13, -9.90367510222967e-13, -9.90624249297412e-13,
      -9.90194037875369e-13, -9.91373649839034e-13, -9.91873250200115e-13,
      -9.91193238597532e-13, -9.91456916565880e-13, -9.89944237694829e-13,
      -9.91096094082877e-13, -9.90901805053568e-13, -9.90263426814408e-13,
      -9.91484672141496e-13, -9.89069937062936e-13, -9.91901005775731e-13,
      -9.91706716746421e-13, -9.88847892458011e-13, -9.92234072683118e-13,
      -9.91595694443959e-13, -9.92872450922278e-13, -9.94343496429906e-13,
      -9.91484672141496e-13, -9.91318138687802e-13, -9.88431558823777e-13,
      7.37209002588238, -15.5030002317990, 4.68518434194794, 5.58066496329142,
      -0.404449011712775, -1.39137323046436, -0.339116857120903 };

  const float valid_desthi2[] = { -2.80091227988777e-12, -2.79960776783383e-12,
      -2.80357681514687e-12, -2.80355599846516e-12, -2.80095391325119e-12,
      -2.79949674553137e-12, -2.79951062331918e-12, -2.80001022368026e-12,
      -2.80267475893936e-12, -2.79856693374825e-12, -2.80492296056423e-12,
      -0.0781250000022623, 0.164291522328916, 0.634073488075181,
      -1.49696584171718, -2.62270640553024, 6.97048991951669,
      13.4936761845669, -2.98585954495631, -19.8119363515072,
      -12.7098068594040, 1.52245837263813, 7.82528131630407,
      8.59130932663576, 5.24090543738087, 1.01894438076528,
      -1.16818198731391, -1.89266864772546, -1.51961243979140,
      -0.776900347899835, -0.320541522330983, -0.0781250000022604 };

  const float valid_destlo[] = { 6.03235928067132, 8.03235928067132,
      10.0323592806713, 12.0323592806713, 14.0323592806713, 16.0323592806713,
      18.0323592806713, 20.0323592806713, 22.0323592806713, 24.0323592806713,
      26.0323592806713, 28.0287655230843, 30.0399167066535, 32.0615267227001,
      33.9634987065767, 35.9320147305194, 38.3103125658258, 40.4883104236778,
      42.2839848729069, 43.7345002903498, 43.7794736932925, 45.1480484137191,
      49.8652419127137, 55.7384062022009, 62.7058766150960, 65.2835749751486,
      58.7895581326311, 46.7708694321525, 31.0673425771182, 16.9214616227404,
      9.00063853315767, 5.73072526035035 };

  stationary_wavelet_apply_na(WAVELET_TYPE_DAUBECHIES, 8, 1, array, length,
                              desthi1, destlo);
  memcpy(array, destlo, sizeof(array));
  stationary_wavelet_apply_na(WAVELET_TYPE_DAUBECHIES, 8, 2, array, length,
                              desthi2, destlo);
   for (int i = 0; i < length; i++) {
     ASSERT_NEAR(desthi1[i], valid_desthi1[i], 1e-5) << i;
     ASSERT_NEAR(desthi2[i], valid_desthi2[i], 1e-5) << i;
     ASSERT_FLOAT_EQ(destlo[i], valid_destlo[i]) << i;
   }
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

class StationaryWaveletTest : public ::testing::TestWithParam<
    std::tuple<WaveletType, int, int>> {
public:
 typedef std::unique_ptr<float, decltype(&std::free)> FloatPtr;

 protected:
  StationaryWaveletTest() : length_(512), array_(nullptr, std::free),
      prep_(nullptr, std::free), desthi_(nullptr, std::free),
      destlo_(nullptr, std::free) {
  }

  virtual void SetUp() override {
    array_ = std::uniquify(mallocf(length_), std::free);
    for (size_t i = 0; i < length_; i++) {
      array_.get()[i] = i;
    }
    desthi_ = std::uniquify(mallocf(length_), std::free);
    destlo_ = std::uniquify(mallocf(length_), std::free);
  }

  size_t length_;
  FloatPtr array_;
  FloatPtr prep_;
  FloatPtr desthi_;
  FloatPtr destlo_;
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
                            std::get<2>(GetParam()),
                           array_.get(), length_, desthi_.get(), destlo_.get());
  float validdesthi[length_], validdestlo[length_];
  stationary_wavelet_apply_na(std::get<0>(GetParam()), std::get<1>(GetParam()),
                              std::get<2>(GetParam()),
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
        ::testing::Values(2, 4, 6, 8, 12, 16)
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
        ::testing::Values(2, 4, 6, 8, 12, 16),
        ::testing::Values(1, 2, 3, 4)
    ));

INSTANTIATE_TEST_CASE_P(
    Coiflets, StationaryWaveletTest,
    ::testing::Combine(
        ::testing::Values(WAVELET_TYPE_COIFLET),
        ::testing::Values(6, 12),
        ::testing::Values(1, 2)
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
