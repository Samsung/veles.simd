# veles.simd
Various mathematical routines with SIMD acceleration (SSE/AVX/NEON) in the form of a compact C library.

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
