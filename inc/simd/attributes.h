/*! @file attributes.h
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

#ifndef INC_SIMD_ATTRIBUTES_H_
#define INC_SIMD_ATTRIBUTES_H_

#ifndef INLINE
/* Set modifiers to always inline the code */
#define INLINE static __attribute__((always_inline)) inline
#endif

#ifndef MALLOC
/* malloc() function attribute */
#define MALLOC __attribute__ ((__malloc__))
#endif

#ifndef LIBSIMD_IMPLEMENTATION

#ifndef NOTNULL
/* Mark pointer parameters which must not be NULL */
#define NOTNULL(...) __attribute__ ((__nonnull__ (__VA_ARGS__)))
#endif

#else
#define NOTNULL(...)
#endif

#ifndef UNUSED
/* Mark parameters as unused to avoid warnings */
#define UNUSED __attribute__((unused))
#endif

#ifndef WARN_UNUSED_RESULT
/* warn about unused result function attribute */
#define WARN_UNUSED_RESULT __attribute__ ((__warn_unused_result__))
#endif

#endif  // INC_SIMD_ATTRIBUTES_H_
