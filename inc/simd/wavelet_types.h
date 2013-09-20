/*! @file wavelet_types.h
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

#ifndef INC_SIMD_WAVELET_TYPES_H_
#define INC_SIMD_WAVELET_TYPES_H_

#include <simd/common.h>

SIMD_API_BEGIN

typedef enum {
  WAVELET_TYPE_DAUBECHIES,
  WAVELET_TYPE_COIFLET,
  WAVELET_TYPE_SYMLET
} WaveletType;

SIMD_API_END


#endif  // INC_SIMD_WAVELET_TYPES_H_
