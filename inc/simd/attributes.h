/*! @file attributes.h
 *  @brief New file description.
 *  @author Markovtsev Vadim <v.markovtsev@samsung.com>
 *  @version 1.0
 *
 *  @section Notes
 *  This code partially conforms to <a href="http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml">Google C++ Style Guide</a>.
 *
 *  @section Copyright
 *  Copyright © 2013 Samsung R&D Institute Russia
 *
 *  @section License
 *  Licensed to the Apache Software Foundation (ASF) under one
 *  or more contributor license agreements.  See the NOTICE file
 *  distributed with this work for additional information
 *  regarding copyright ownership.  The ASF licenses this file
 *  to you under the Apache License, Version 2.0 (the
 *  "License"); you may not use this file except in compliance
 *  with the License.  You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an
 *  "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 *  KIND, either express or implied.  See the License for the
 *  specific language governing permissions and limitations
 *  under the License.
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

#ifndef NOTNULL
#ifndef LIBSIMD_IMPLEMENTATION
/* Mark pointer parameters which must not be NULL */
#define NOTNULL(...) __attribute__ ((__nonnull__ (__VA_ARGS__)))
#else
#define NOTNULL(...)
#endif
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
