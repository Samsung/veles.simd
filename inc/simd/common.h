/*! @file common.h
 *  @brief Common declarations
 *  @author Ernesto Sanches <ernestosanches@gmail.com>
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

#ifndef INC_SIMD_COMMON_H_
#define INC_SIMD_COMMON_H_

#ifdef __cplusplus
#define EXTERN_C_BEGIN extern "C" {
#define EXTERN_C_END }
#else
#define EXTERN_C_BEGIN
#define EXTERN_C_END
#endif

#if __GNUC__ >= 4
#define VISIBILITY_DEFAULT _Pragma("GCC visibility push(default)")
#define VISIBILITY_POP _Pragma("GCC visibility pop")
#else
#define VISIBILITY_DEFAULT
#define VISIBILITY_POP
#endif

#define SIMD_API_BEGIN EXTERN_C_BEGIN VISIBILITY_DEFAULT
#define SIMD_API_END EXTERN_C_END VISIBILITY_POP

#ifndef IMPLEMENTATION
#define NOTNULL(...) __attribute__ ((__nonnull__ (__VA_ARGS__)))
#else
#define NOTNULL(...)
#endif

#endif  // INC_SIMD_COMMON_H_
