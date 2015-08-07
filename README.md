# veles.simd
Various mathematical routines with SIMD acceleration (SSE/AVX/NEON) in the form of a compact C library.

### Implemented features

*  Conversion between int16_t, int32_t and float
*  Real and complex vector multiplication, addition
*  1D convolution and correlation with best approach detection (naive, overlap-save, FFT)
*  1D peak detection
*  sin, cos, log, exp (delegated to [AVX mathfun](http://software-lisc.fbk.eu/avx_mathfun/) and [NEON mathfun](http://gruntthepeon.free.fr/ssemath/neon_mathfun.html))
*  1D and 2D normalization
*  1D decimated and stationary (undecimated) wavelets

### Building
```
./autogen.sh
mkdir build && cd build
../configure
make -j$(getconf _NPROCESSORS_ONLN)
make instal DESTDIR=...
```

By default, this library makes use of [FFTF](https://github.com/Samsung/FFTF).
You can pass ``--disable-simd-fftf`` to ``configure`` to skip building dependent features.

### Copyright
Copyright © 2013 Samsung R&D Institute Russia

### License
[Apache 2.0](http://www.apache.org/licenses/LICENSE-2.0).
