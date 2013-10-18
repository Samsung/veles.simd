/*! @file detect_peaks.h
 *  @brief Find extrema in 1D signal.
 *  @author Vadim Markovtsev <v.markovtsev@samsung.com>
 *  @version 1.0
 *
 *  @section Notes
 *  This code partially conforms to <a href="http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml">Google C++ Style Guide</a>.
 *
 *  @section Copyright
 *  Copyright 2013 Samsung R&D Institute Russia
 */

#ifndef INC_SIMD_DETECT_PEAKS_H_
#define INC_SIMD_DETECT_PEAKS_H_

#include <stddef.h>
#include <simd/attributes.h>
#include <simd/common.h>

SIMD_API_BEGIN

typedef enum {
  kExtremumTypeMaximum = 1,
  kExtremumTypeMinimum,
  kExtremumTypeBoth
} ExtremumType;

typedef struct {
  int position;
  float value;
} ExtremumPoint;

/// @brief Extract maximums and minimums from the series of floating point
/// numbers.
/// @param simd Value indicating whether to use SIMD acceleration.
/// @param data The array of floating point numbers representing the signal.
/// @param size The length of the array (in float-s, not in bytes).
/// @param type The type of the extracted extrema.
/// @param results The pointer to the array of ExtremumPoint-s. That array
/// will be allocated with malloc(), so it should be disposed with free()
/// after it's been used. If no points are found, it is set to NULL.
/// @param resultsLength The number of found extremum points.
void detect_peaks(int simd, const float *data, size_t size, ExtremumType type,
                  ExtremumPoint **results, size_t *resultsLength)
    NOTNULL(2, 5, 6);

SIMD_API_END

#endif  // INC_SIMD_DETECT_PEAKS_H_
