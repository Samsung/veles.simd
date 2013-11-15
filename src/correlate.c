/*! @file correlate.c
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

#ifndef NO_FFTF
#define LIBSIMD_IMPLEMENTATION
#include "inc/simd/correlate.h"
#include "inc/simd/convolve.h"
#include "inc/simd/arithmetic-inl.h"

CrossCorrelationFFTHandle cross_correlate_fft_initialize(size_t xLength,
                                                         size_t hLength) {
  CrossCorrelationFFTHandle handle = convolve_fft_initialize(xLength, hLength);
  handle.reverse = 1;
  return handle;
}

void cross_correlate_fft(CrossCorrelationFFTHandle handle,
                         const float *x, const float *h,
                         float *result) {
  convolve_fft(handle, x, h, result);
}

void cross_correlate_fft_finalize(CrossCorrelationFFTHandle handle) {
  convolve_fft_finalize(handle);
}

CrossCorrelationOverlapSaveHandle cross_correlate_overlap_save_initialize(
    size_t xLength, size_t hLength) {
  CrossCorrelationOverlapSaveHandle handle =
      convolve_overlap_save_initialize(xLength, hLength);
  handle.reverse = 1;
  return handle;
}

void cross_correlate_overlap_save(CrossCorrelationOverlapSaveHandle handle,
                            const float *__restrict x,
                            const float *__restrict h,
                            float *result) {
  convolve_overlap_save(handle, x, h, result);
}

void cross_correlate_overlap_save_finalize(
    CrossCorrelationOverlapSaveHandle handle) {
  convolve_overlap_save_finalize(handle);
}

void cross_correlate_simd(int simd,
                          const float *__restrict x, size_t xLength,
                          const float *__restrict h, size_t hLength,
                          float *__restrict result) {
  for (int n = hLength - 1; n > -(int)xLength; n--) {
    float sum = 0.f;
    int beg = n <= 0? -n : 0;
    int end = -n + hLength;
    if (end > (int)xLength) {
      end = (int)xLength;
    }
    if (simd) {
#ifdef __AVX__
      int simdEnd = beg + ((end - beg) & ~7);
      __m256 accum = _mm256_setzero_ps();
      for (int m = beg; m < simdEnd; m += 8) {
        __m256 xvec = _mm256_loadu_ps(x + m);
        __m256 hvec = _mm256_loadu_ps(h + n + m);
        __m256 mulres = _mm256_mul_ps(xvec, hvec);
        accum = _mm256_add_ps(accum, mulres);
      }
      accum = _mm256_hadd_ps(accum, accum);
      accum = _mm256_hadd_ps(accum, accum);
      sum = accum[0] + accum[4];
      for (int m = simdEnd; m < end; m++) {
        sum += x[m] * h[n + m];
      }
    } else {
#elif defined(__ARM_NEON__)
      int simdEnd = beg + ((end - beg) & ~3);
      float32x4_t accum = vdupq_n_f32(0.f);
      for (int m = beg; m < simdEnd; m += 4) {
        float32x4_t xvec = vld1q_f32(x + m);
        float32x4_t hvec = vld1q_f32(h + n + m);
        accum = vmlaq_f32(accum, xvec, hvec);
      }
      float32x2_t accum2 = vpadd_f32(vget_high_f32(accum),
                                     vget_low_f32(accum));
      sum = vget_lane_f32(accum2, 0) + vget_lane_f32(accum2, 1);
      for (int m = simdEnd; m < end; m++) {
        sum += x[m] * h[n + m];
      }
    } else {
#else
    } {
#endif
      for (int m = beg; m < end; m++) {
        sum += x[m] * h[n + m];
      }
    }
    result[-n + hLength - 1] = sum;
  }
}

CrossCorrelationHandle cross_correlate_initialize(size_t xLength,
                                                size_t hLength) {
  CrossCorrelationHandle handle = convolve_initialize(xLength, hLength);
  switch (handle.algorithm) {
    case kConvolutionAlgorithmFFT:
      handle.handle.fft.reverse = 1;
      break;
    case kConvolutionAlgorithmOverlapSave:
      handle.handle.os.reverse = 1;
      break;
    case kConvolutionAlgorithmBruteForce:
      break;
  }
  return handle;
}

void cross_correlate(CrossCorrelationHandle handle,
                     const float *__restrict x, const float *__restrict h,
                     float *__restrict result) {
  switch (handle.algorithm) {
    case kConvolutionAlgorithmFFT:
    case kConvolutionAlgorithmOverlapSave:
      convolve(handle, x, h, result);
      break;
    case kConvolutionAlgorithmBruteForce:
      cross_correlate_simd(1, x, handle.x_length, h, handle.h_length, result);
      break;
  }
}

void cross_correlate_finalize(CrossCorrelationHandle handle) {
  convolve_finalize(handle);
}
#endif
