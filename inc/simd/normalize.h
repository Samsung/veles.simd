/*! @file normalize.h
 *  @brief Image plane normalization routines.
 *  @author Vadim Markovtsev <v.markovtsev@samsung.com>
 *  @version 1.0
 *
 *  @section Notes
 *  This code partially conforms to <a href="http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml">Google C++ Style Guide</a>.
 *
 *  @section Copyright
 *  Copyright 2013 Samsung R&D Institute Russia
 */

#ifndef INC_SIMD_NORMALIZE_H_
#define INC_SIMD_NORMALIZE_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/// @brief Performs the plane normalization [0, 255] -> [-1, 1].
/// @param simd Value indicating whether to use available SIMD acceleration.
/// @param src The source byte array, stored in row-major format.
/// @param src_stride The stride (the actual width) of the plane.
/// @param width The width of the plane.
/// @param height The height of the plane.
/// @param dst The resulting floating point array.
/// @param dst_stride The stride of dst.
void normalize2D(int simd, const uint8_t* src, int src_stride,
                 int width, int height, float* dst, int dst_stride);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // INC_SIMD_NORMALIZE_H_
