/*! @file convolute.h
 *  @brief Calculates the linear convolution of two signals.
 *  @author Markovtsev Vadim <v.markovtsev@samsung.com>
 *  @version 1.0
 *
 *  @section Notes
 *  This code partially conforms to <a href="http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml">Google C++ Style Guide</a>.
 *
 *  @section Copyright
 *  Copyright 2013 Samsung R&D Institute Russia
 */

#ifndef INC_SIMD_CONVOLUTE_H_
#define INC_SIMD_CONVOLUTE_H_

#include <stddef.h>
#include <simd/common.h>
#include <simd/attributes.h>
#include <simd/convolute_structs.h>

SIMD_API_BEGIN

typedef struct ConvoluteFFTHandle ConvoluteFFTHandle;

/// @brief Prepares for the calculation of linear convolution of two signals
/// using the FFT method.
/// @param xLength The length of the first array in float-s.
/// @param hLength The length of the second array in float-s.
/// @return The handle for convolute_fft().
ConvoluteFFTHandle convolute_fft_initialize(size_t xLength, size_t hLength);

/// @brief Calculates the linear convolution of two signals using
/// the FFT method.
/// @param handle The structure obtained from convolute_fft_initialize().
/// @param x The first signal (long one).
/// @param h The second signal (short one).
/// @param result The resulting signal of length xLength.
/// @note result and x may be the same arrays.
void convolute_fft(ConvoluteFFTHandle handle,
                   const float *x, const float *h,
                   float *result) NOTNULL(2, 3, 4);

/// @brief Frees any resources allocated by convolute_fft_initialize().
/// @param handle The structure obtained from convolute_fft_initialize().
void convolute_fft_finalize(ConvoluteFFTHandle handle);

typedef struct ConvoluteOverlapSaveHandle ConvoluteOverlapSaveHandle;

/// @brief Prepares for the calculation of linear convolution of two signals
/// using the overlap-save method.
/// @param xLength The length of the first array in float-s.
/// @param hLength The length of the second array in float-s.
/// @return The handle for convolute_overlap_save().
ConvoluteOverlapSaveHandle convolute_overlap_save_initialize(
    size_t xLength, size_t hLength);

/// @brief Calculates the linear convolution of two signals using
/// the overlap-save method.
/// @param handle The structure obtained from convolute_overlap_save_initialize().
/// @param x The first signal (long one).
/// @param h The second signal (short one).
/// @param result The resulting signal of length xLength.
/// @note result and x may be the same arrays.
void convolute_overlap_save(ConvoluteOverlapSaveHandle handle,
                            const float *__restrict x,
                            const float *__restrict h,
                            float *result) NOTNULL(2, 3, 4);

/// @brief Frees any resources allocated by convolute_overlap_save_initialize().
/// @param handle The structure obtained from convolute_overlap_save_initialize().
void convolute_overlap_save_finalize(ConvoluteOverlapSaveHandle handle);

/// @brief Calculates the linear convolution of two signals using
/// the "brute force" method.
/// @param simd Value indicating whether to use SIMD acceleration or not.
/// @param x The first signal (long one).
/// @param xLength The length of the first array in float-s.
/// @param h The second signal (short one).
/// @param hLength The length of the second array in float-s.
/// @param result The resulting signal of length xLength.
/// @note result and x may be the same arrays.
void convolute_simd(int simd,
                    const float *x, size_t xLength,
                    const float *h, size_t hLength,
                    float *result) NOTNULL(2, 4, 6);

typedef struct ConvoluteHandle ConvoluteHandle;

/// @brief Prepares for the calculation of linear convolution of two signals
/// using the best method.
/// @param xLength The length of the first array in float-s.
/// @param hLength The length of the second array in float-s.
/// @return The handle for convolute().
ConvoluteHandle convolute_initialize(size_t xLength, size_t hLength);

/// @brief Calculates the linear convolution of two signals using
/// the best method.
/// @param handle The structure obtained from convolute_initialize().
/// @param x The first signal (long one).
/// @param h The second signal (short one).
/// @param result The resulting signal of length xLength.
/// @note result and x may be the same arrays.
void convolute(ConvoluteHandle handle,
               const float *x, const float *h,
               float *result) NOTNULL(2, 3, 4);

/// @brief Frees any resources allocated by convolute_overlap_initialize().
/// @param handle The structure obtained from convolute_overlap_initialize().
void convolute_finalize(ConvoluteHandle handle);

SIMD_API_END

#endif  // INC_SIMD_CONVOLUTE_H_
