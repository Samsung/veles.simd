#  This file is a part of SEAPT, Samsung Extended Autotools Project Template

#  Copyright (c) 2012,2013 Samsung R&D Institute Russia
#  All rights reserved.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met: 
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer. 
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
#  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

AC_INIT(Simd, 1.0)
m4_include([common.ac])

# Check whether to build FFTF dependent code
AC_ARG_ENABLE([enable-simd-fftf],
    AS_HELP_STRING([--disable-simd-fftf], [Disable FFTF dependent code])
)
AS_IF([test "x$enable_simd_fftf" == "xno"], [    
    CPPFLAGS="$CPPFLAGS -DNO_FFTF"
])

# Add an optional feature
AC_ARG_WITH([built-in-fftf],
    AS_HELP_STRING([--with-built-in-fftf],
    [use FFTF build directory from built-in-fftf file]),
    [FFTF_LIBS="\$(abs_top_builddir)/$(cat built-in-fftf)/src/libfftf.la"]
)

AS_IF([test "x$with_built_in_fftf" != "xyes" -a "x$enable_simd_fftf" != "xno"], [
	# Check for installed FFTF
	PKG_CHECK_MODULES([FFTF], [fftf >= 1.0])
])

# Check whether to conduct test benchmarks
AC_ARG_ENABLE([benchmarks],
    AS_HELP_STRING([--enable-benchmarks], [execute SIMD speedup benchmarks during tests evaluation])
)
AS_IF([test "x$enable_benchmarks" == "xyes"], [
    CPPFLAGS="$CPPFLAGS -DBENCHMARK"
])

AC_OUTPUT
