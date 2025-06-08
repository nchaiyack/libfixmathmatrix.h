# `libfixmathmatrix_final.h`

A (nearly) header-only amalgamation of `libfixmath` and `libfixmatrix` fixed-point math libraries as of June 2025. These original libraries were created by Petteri Aimonen and friends (see `AUTHORS`), who implemented all of the numerics.  Both of these are fixed-point math libraries. The `libfixmath` component implements the core type `fix16_t`, Q16.16 fixed-point number (i.e. 32-bit width with 16-bit integer and 16-bit fractional parts, all twos' complement). The `libfixmathmatrix` component handles matrices up to row and column length `FIXMATRIX_MAX_SIZE` (allocating the largest possible size each time.)

The capabilities of each library include:

1. `libfixmath`: 
    • **Mathematical operations** including basic arithmetic (add, subtract, multiply, divide, modulo), trigonometric functions (sin, cos, tan, asin, acos, atan, atan2), and advanced functions (sqrt, exp, log, interpolation) with optional overflow detection and lookup table optimizations.
    • **String conversion utilities** for converting between fix16_t values and human-readable decimal representations.
    • **Signal processing capabilities** including a real-input Fast Fourier Transform (FFT) implementation and auxiliary data types (fract32_t for fractions, utility functions for 64-bit operations).
    • Of separate interest, an implementation of **64-bit integers** for platforms that don't have a C integral type for them (i.e. `long long`.) 

2. `libfixmathmatrix':
    • **Matrix operations** including multiplication, addition, subtraction, transpose, scalar operations, and advanced linear algebra (QR decomposition, Cholesky decomposition, matrix inversion, solving linear systems).
    • **Vector arithmetic** for both 2D and 3D vectors with basic operations (add, subtract, multiply/divide by scalar), norms, normalization, dot products, cross products (3D), and rotation (2D).
    • **Quaternion support** for 3D rotations with operations including multiplication, conjugation, normalization, axis-angle conversion, rotation matrix conversion, and vector rotation.
    • **Error handling framework** with automatic propagation of overflow, dimension mismatch, usage errors, and singular matrix detection across all operations.

**Caveat.** Numerics are hard. I take the importance of their correctness *very, very* seriously. Details of the test suite are given below. I will be transparent about my goals for `libfixmathmatrix.h`: it is designed to do graphics on platforms without FPUs, as well as to use its types as a binary interchange format (up to endianness.) Do *not* expect the behavior of computations to be consistent between
plaforms. Do *not* put it on an MCU and do mission-critical things with it. Per the MIT license: `THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED.`

**Usage.** 

1. **Including and instantiating `libfixmathmatrix_final.h`** There are multiple ways you can include `libfixmathmatrix_final.h` in a file. First, you must make sure that, in addition to `libfixmathmatrix_final.h`, you make accessible its its sub-headers (split out for readability) `implementation_group_[a-g].h` -- still all headers. Then you must make a choice about how to instantiate `libfixmathmatrix_final.h`:

    - *To share one instantiation* of `libfixmathmatrix_final.h` between compilation units, simply `#include <libfixmathmatrix_final.h>` as you wish (e.g. in a shared header); this is good for shared header files. However, you must then choose a compilation unit in which to instantiate `libfixmathmatrix_final.h`, by doing a  `#define LIBFIXMATHMATRIX_IMPLEMENTATION` before a second `#include <libfixmathmatrix_final.h>`.
    
    - *To statically instantiate* `libfixmathmatrix_final.h` in a compilation unit, `#define LIBFIXMATHMATRIX_IMPLEMENTATION` and `#define LIBFIXMATHMATRIX_STATIC` before `#include <libfixmathmatrix_final.h>`. Note that you can therefore override configuration options (see below.)

2. **Important: lookup tables you may optionally choose to compile.** *By default, but entirely optionally,* `FIXMATH_NO_CACHE` is *not* defined, and `FIXMATH_SIN_LUT` *is* defined. As a result, you'll also need to compile the source files `libfixmathmatrix_cache.c` and `libfixmathmatrix_lut.h`. You don't 

**Verification, quality, and testing.**

Here are tests that are currently available in the `tests` subfolder (although you'll currently have to compile them yourself.)

    1. Original tests from 'libfixmath' and 'libfixmathmatrix.h', which test all of the functionalities of the separate files.

    2. `tests/multi-unit-configuration-option-tests`, a test of the different instantiation modes (static vs shared).

    3. `tests/configuration-options`, multiple tests that exercise other compile options, including optimizations and 64-bit emulation;
    <tests_cat2></tests_cat2>

    4. `test/mathematical-properties`, more extensive tests of mathematical identities and properties, some using randomized examples.

**Caveat**. Cursor (i.e. `claude-4-sonnet`) was used to help untangle function dependencies for both `libfixmath` and `libfixmathmatrix`, and thereby "linearize" definitions in all of the files to amalgamate them in one file. Cursor wrote tests to specifications; but "trivial" tests of algebraic identities were entirely written by Cursor. Finally, Cursor helped write descriptions in this README that required looking at all of the functionality implemented by a subset of the codebase.

## Original readmes for `libfixmath` and `libfixmatrix`

These are copied verbatim from their repositories as of pull date; don't worry about any build instructions.

### libfixmath

This is a mirror of the libfixmath's original SVN repository on Google Code.

**Not actively maintained, pull requests welcome.**

Libfixmath implements Q16.16 format fixed point operations in C.

License: <a href="http://www.opensource.org/licenses/mit-license.php">MIT</a>

##### Options

Configuration options are compile definitions that are checked by the preprocessor with `#ifdef` and `#ifndef`.  All of these are undefined by default.

###### `FIXMATH_FAST_SIN`

- `#ifndef`: Most accurate version, accurate to ~2.1%.
- `#ifdef`: Fast implementation, runs at 159% the speed of above 'accurate' version with a slightly lower accuracy of ~2.3%.

###### `FIXMATH_NO_64BIT`

- `#ifndef`: For compilers/platforms that have `uint64_t`.
- `#ifdef`: For compilers/platforms that do not have `uint64_t`.

###### `FIXMATH_NO_CACHE`

- `#ifndef`: Use static memory caches for exponents (32KB) and trigonometry (80KB). 
- `#ifdef`: Do not use caches.

###### `FIXMATH_NO_HARD_DIVISION`

Note: will be automatically defined if `FIXMATH_OPTIMIZE_8BIT` is defined.

- `#ifndef`: For platforms that have hardware integer division.
- `#ifdef`: For platforms that do not have hardware integer division.

###### `FIXMATH_NO_OVERFLOW`

- `#ifndef`: Check for overflow and return the overflow constants. 
- `#ifdef`: Do not check for overflow.

###### `FIXMATH_NO_ROUNDING`

- `#ifndef`: Use rounding. 
- `#ifdef`: Do not use rounding.

###### `FIXMATH_OPTIMIZE_8BIT`

- `#ifndef`: Do not optimize for processors with 8-bit multiplication like Atmel AVR. 
- `#ifdef`: Optimize for processors like Atmel AVR.  Also defines `FIXMATH_NO_HARD_DIVISION` automatically in `fix16.h`.

##### Include the `libfixmath` library in your CMake Project

The simplest way to use `libfixmath` as a dependency is with CMake's [FetchContent API](https://cmake.org/cmake/help/latest/module/FetchContent.html).

```cmake
include(FetchContent)
FetchContent_Declare(
    libfixmath
    GIT_REPOSITORY https://github.com/PetteriAimonen/libfixmath.git
    GIT_TAG <the long git hash of the version you want>
)
FetchContent_MakeAvailable(libfixmath)

target_compile_definitions(libfixmath PRIVATE
    # FIXMATH_FAST_SIN
    # FIXMATH_NO_64BIT
    # FIXMATH_NO_CACHE
    # FIXMATH_NO_HARD_DIVISION
    # FIXMATH_NO_OVERFLOW
    # FIXMATH_NO_ROUNDING
    # FIXMATH_OPTIMIZE_8BIT
)

target_link_libraries(my_cmake_project PRIVATE libfixmath)
```

### Fixed point matrix library

Libfixmatrix is a matrix computation library for microcontrollers.
It is based on the libfixmath_ library, which uses 16.16 bit fixed point values.
The main focus is processors without an FPU, such as ARM Cortex-M3.
The compiled size of the library is less than 5 kB, depending on optimization settings and processor.

The library includes all basic matrix operations, such as multiplication, addition and transposition.
Matrix equation solving (and matrix inversion) is implemented through `QR decomposition`_.
Also `Cholesky decomposition`_ is included. See `function reference`_ for details.

To avoid complexity and dynamic memory allocations, all matrices are allocated a buffer with constant size, specified with parameter 
`FIXMATRIX_MAX_SIZE`. This wastes some memory with matrices smaller than the maximum size, but allows more predictable memory usage.

Libfixmatrix is suited well for tasks involving small matrices (often less than 10x10):
`Kalman filters`_, `transformation matrices`_ and solving `systems of linear equations`_.

.. _libfixmath: http://code.google.com/p/libfixmath/
.. _QR decomposition: http://en.wikipedia.org/wiki/QR_decomposition
.. _Cholesky decomposition: http://en.wikipedia.org/wiki/Cholesky_decomposition
.. _function reference: https://github.com/PetteriAimonen/libfixmatrix/blob/master/FUNCTIONS.rst
.. _Kalman filters: http://en.wikipedia.org/wiki/Kalman_filter
.. _transformation matrices: http://en.wikipedia.org/wiki/Transformation_matrix
.. _systems of linear equations: http://en.wikipedia.org/wiki/System_of_linear_equations
