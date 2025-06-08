/*
 * libfixmathmatrix.h - Combined fixed-point math and matrix library
 * 
 * This is a header-only library combining libfixmath and libfixmatrix
 * into a single include file for easier integration.
 * 
 * USAGE:
 * Include this header file in your project:
 *   #include "libfixmathmatrix.h"
 * 
 * To include implementations (in ONE source file only):
 *   #define LIBFIXMATHMATRIX_IMPLEMENTATION
 *   #include "libfixmathmatrix.h"
 * 
 * CONFIGURATION:
 * This library preserves all original FIXMATH_* configuration options:
 *   - FIXMATH_NO_OVERFLOW: Disable overflow checking
 *   - FIXMATH_NO_64BIT: Use 32-bit implementations
 *   - FIXMATH_NO_ROUNDING: Disable rounding for better performance
 *   - FIXMATH_OPTIMIZE_8BIT: Optimize for 8-bit platforms
 *   - FIXMATH_NO_CACHE: Disable caching in trigonometric functions
 *   - FIXMATH_NO_HARD_DIVISION: Use soft division algorithms
 *   - FIXMATH_FAST_SIN: Use fast sine implementation
 *   - FIXMATH_SIN_LUT: Use lookup table for sine
 *   - FIXMATH_NO_CTYPE: Disable ctype.h usage
 *   - FIXMATRIX_MAX_SIZE: Maximum matrix size (default 8)
 *   - FIXSTRING_NO_STDIO: Disable stdio functions
 * 
 * EXTERNAL DEPENDENCIES:
 * Link with libfixmathmatrix_lut.c for trigonometric lookup tables.
 * Link with libfixmathmatrix_cache.c for performance caches.
 */

#ifndef LIBFIXMATHMATRIX_H
#define LIBFIXMATHMATRIX_H

#ifdef __cplusplus
extern "C"
{
#endif

/* Standard library includes */
#ifdef __KERNEL__
#include <linux/types.h>
#else
#include <stdint.h>
#endif
#include <stdbool.h>
#include <string.h>
#ifndef FIXSTRING_NO_STDIO
#include <stdio.h>
#endif
#if defined(FIXMATH_NO_CTYPE) || defined(__KERNEL__)
/* Provide local implementations when ctype.h is not available */
#else
#include <ctype.h>
#endif

/* Function attributes for optimization */
#ifndef FIXMATH_FUNC_ATTRS
	#if defined(LIBFIXMATHMATRIX_STATIC)
		#define FIXMATH_FUNC_ATTRS static
	#else
		#define FIXMATH_FUNC_ATTRS extern
	#endif
#endif

/* Automatically define FIXMATH_NO_HARD_DIVISION to maintain backwards
 * compatibility with usage of FIXMATH_OPTIMIZE_8BIT.
 */
#if defined(FIXMATH_OPTIMIZE_8BIT)
#  define FIXMATH_NO_HARD_DIVISION
#endif

/*===========================================================================*/
/* TYPE DEFINITIONS                                                          */
/*===========================================================================*/

/* Basic fixed-point types */
typedef int32_t fix16_t;
typedef uint32_t fract32_t;

/* 64-bit integer handling */
#ifndef FIXMATH_NO_64BIT
static inline  int64_t int64_const(int32_t hi, uint32_t lo) { return (((int64_t)hi << 32) | lo); }
static inline  int64_t int64_from_int32(int32_t x) { return (int64_t)x; }
static inline  int32_t int64_hi(int64_t x) { return (x >> 32); }
static inline uint32_t int64_lo(int64_t x) { return (x & ((1ULL << 32) - 1)); }

static inline int64_t int64_add(int64_t x, int64_t y)   { return (x + y);  }
static inline int64_t int64_neg(int64_t x)              { return (-x);     }
static inline int64_t int64_sub(int64_t x, int64_t y)   { return (x - y);  }
static inline int64_t int64_shift(int64_t x, int8_t y)  { return (y < 0 ? (x >> -y) : (x << y)); }

static inline int64_t int64_mul_i32_i32(int32_t x, int32_t y) { return ((int64_t)x * y);  }
static inline int64_t int64_mul_i64_i32(int64_t x, int32_t y) { return (x * y);  }

static inline int64_t int64_div_i64_i32(int64_t x, int32_t y) { return (x / y);  }

static inline int int64_cmp_eq(int64_t x, int64_t y) { return (x == y); }
static inline int int64_cmp_ne(int64_t x, int64_t y) { return (x != y); }
static inline int int64_cmp_gt(int64_t x, int64_t y) { return (x >  y); }
static inline int int64_cmp_ge(int64_t x, int64_t y) { return (x >= y); }
static inline int int64_cmp_lt(int64_t x, int64_t y) { return (x <  y); }
static inline int int64_cmp_le(int64_t x, int64_t y) { return (x <= y); }
#else
/* 32-bit fallback for platforms without 64-bit support */
typedef struct {
	 int32_t hi;
	uint32_t lo;
} _int64_t;

static inline _int64_t int64_const(int32_t hi, uint32_t lo) { return (_int64_t){ hi, lo }; }
static inline _int64_t int64_from_int32(int32_t x) { return (_int64_t){ (x < 0 ? -1 : 0), x }; }
static inline   int32_t int64_hi(_int64_t x) { return x.hi; }
static inline  uint32_t int64_lo(_int64_t x) { return x.lo; }

static inline int int64_cmp_eq(_int64_t x, _int64_t y) { return ((x.hi == y.hi) && (x.lo == y.lo)); }
static inline int int64_cmp_ne(_int64_t x, _int64_t y) { return ((x.hi != y.hi) || (x.lo != y.lo)); }
static inline int int64_cmp_gt(_int64_t x, _int64_t y) { return ((x.hi > y.hi) || ((x.hi == y.hi) && (x.lo >  y.lo))); }
static inline int int64_cmp_ge(_int64_t x, _int64_t y) { return ((x.hi > y.hi) || ((x.hi == y.hi) && (x.lo >= y.lo))); }
static inline int int64_cmp_lt(_int64_t x, _int64_t y) { return ((x.hi < y.hi) || ((x.hi == y.hi) && (x.lo <  y.lo))); }
static inline int int64_cmp_le(_int64_t x, _int64_t y) { return ((x.hi < y.hi) || ((x.hi == y.hi) && (x.lo <= y.lo))); }

static inline _int64_t int64_add(_int64_t x, _int64_t y) {
	_int64_t ret;
	ret.hi = x.hi + y.hi;
	ret.lo = x.lo + y.lo;
	if((ret.lo < x.lo) || (ret.lo < y.lo))
		ret.hi++;
	return ret;
}

static inline _int64_t int64_neg(_int64_t x) {
	_int64_t ret;
	ret.hi = ~x.hi;
	ret.lo = ~x.lo + 1;
	if(ret.lo == 0)
		ret.hi++;
	return ret;
}

static inline _int64_t int64_sub(_int64_t x, _int64_t y) {
	return int64_add(x, int64_neg(y));
}

static inline _int64_t int64_shift(_int64_t x, int8_t y) {
	_int64_t ret = {0,0};
	if(y >= 64 || y <= -64)
		return (_int64_t){ 0, 0 };
	if(y >= 32) {
		ret.hi = (x.lo << (y - 32));
	} 
	else if(y > 0) {
		ret.hi = (x.hi << y) | (x.lo >> (32 - y));
		ret.lo = (x.lo << y);
	}
	else {
		y = -y;
		if(y >= 32){
			ret.lo = (x.hi >> (y - 32));
			ret.hi = (x.hi < 0) ? -1 : 0;
		} else {
			ret.lo = (x.lo >> y) | (x.hi << (32 - y));
		  ret.hi = (x.hi >> y);
		}
	}
	return ret;
}

static inline _int64_t int64_mul_i32_i32(int32_t x, int32_t y) {
	 int16_t hi[2] = { (x >> 16), (y >> 16) };
	uint16_t lo[2] = { (x & 0xFFFF), (y & 0xFFFF) };

	 int32_t r_hi = hi[0] * hi[1];
	 int32_t r_md = (hi[0] * lo[1]) + (hi[1] * lo[0]);
	uint32_t r_lo = lo[0] * lo[1];

	_int64_t r_hilo64 = (_int64_t){ r_hi, r_lo };
	_int64_t r_md64 = int64_shift(int64_from_int32(r_md), 16);

	return int64_add(r_hilo64, r_md64);
}

static inline _int64_t int64_mul_i64_i32(_int64_t x, int32_t y) {
	int neg = ((x.hi ^ y) < 0);
	if(x.hi < 0)
		x = int64_neg(x);
	uint32_t ypos = (y < 0)? (-y) : (y);

	uint32_t _x[4] = { (x.lo & 0xFFFF), (x.lo >> 16), (x.hi & 0xFFFF), (x.hi >> 16) };
	uint32_t _y[2] = { (ypos & 0xFFFF), (ypos >> 16) };

	uint32_t r[4];
	r[0] = (_x[0] * _y[0]);
	r[1] = (_x[1] * _y[0]);
	uint32_t temp_r1 = r[1];
	r[1] += (_x[0] * _y[1]);
	r[2] = (_x[2] * _y[0]) + (_x[1] * _y[1]);
	r[3] = (_x[3] * _y[0]) + (_x[2] * _y[1]);
	// Detect carry bit in r[1]. r[0] can't carry, and r[2]/r[3] don't matter.
	if(r[1] < temp_r1)
		r[3] ++;

	_int64_t middle = int64_shift(int64_const(0, r[1]), 16);
	_int64_t ret;
	ret.lo = r[0];
	ret.hi = (r[3] << 16) + r[2];
	ret = int64_add(ret, middle);
	return (neg ? int64_neg(ret) : ret);
}

static inline _int64_t int64_div_i64_i32(_int64_t x, int32_t y) {
	int neg = ((x.hi ^ y) < 0);
	if(x.hi < 0)
		x = int64_neg(x);
	if(y < 0)
		y = -y;

	_int64_t ret = { (x.hi / y) , (x.lo / y) };
	x.hi = x.hi % y;
	x.lo = x.lo % y;

	_int64_t _y = int64_from_int32(y);

	_int64_t i;
	for(i = int64_from_int32(1); int64_cmp_lt(_y, x); _y = int64_shift(_y, 1), i = int64_shift(i, 1));

	while(x.hi) {
		_y = int64_shift(_y, -1);
		 i = int64_shift(i, -1);
		if(int64_cmp_ge(x, _y)) {
			x = int64_sub(x, _y);
			ret = int64_add(ret, i);
		}
	}

	ret = int64_add(ret, int64_from_int32(x.lo / y));
	return (neg ? int64_neg(ret) : ret);
}

#define int64_t _int64_t
#endif

/* Vector types */
typedef struct {
    fix16_t x;
    fix16_t y;
} v2d;

typedef struct {
	fix16_t x;
	fix16_t y;
	fix16_t z;
} v3d;

/* Matrix type */
#ifndef FIXMATRIX_MAX_SIZE
#define FIXMATRIX_MAX_SIZE 8
#endif

typedef struct {
    uint8_t rows;
    uint8_t columns;
    uint8_t errors;
    fix16_t data[FIXMATRIX_MAX_SIZE][FIXMATRIX_MAX_SIZE];
} mf16;

/* Quaternion type */
typedef struct {
    fix16_t a; // Real part
    fix16_t b; // i
    fix16_t c; // j
    fix16_t d; // k
} qf16;

/*===========================================================================*/
/* EXTERNAL LOOKUP TABLE DECLARATIONS                                       */
/*===========================================================================*/

/* Trigonometric lookup table - implemented in libfixmathmatrix_lut.c */
extern const uint16_t _fix16_sin_lut[];
extern const uint32_t _fix16_sin_lut_count;

/* Performance caches - implemented in libfixmathmatrix_cache.c */
#ifndef FIXMATH_NO_CACHE
extern fix16_t _fix16_sin_cache_index[4096];
extern fix16_t _fix16_sin_cache_value[4096];
extern fix16_t _fix16_atan_cache_index[2][4096];
extern fix16_t _fix16_atan_cache_value[4096];
extern fix16_t _fix16_exp_cache_index[4096];
extern fix16_t _fix16_exp_cache_value[4096];
#endif

/*===========================================================================*/
/* CONSTANTS AND MACROS                                                      */
/*===========================================================================*/

/* Fixed-point constants */
static const fix16_t FOUR_DIV_PI  = 0x145F3;            /*!< Fix16 value of 4/PI */
static const fix16_t _FOUR_DIV_PI2 = 0xFFFF9840;        /*!< Fix16 value of -4/PI² */
static const fix16_t X4_CORRECTION_COMPONENT = 0x399A; 	/*!< Fix16 value of 0.225 */
static const fix16_t PI_DIV_4 = 0x0000C90F;             /*!< Fix16 value of PI/4 */
static const fix16_t THREE_PI_DIV_4 = 0x00025B2F;       /*!< Fix16 value of 3PI/4 */

static const fix16_t fix16_maximum  = 0x7FFFFFFF; /*!< the maximum value of fix16_t */
static const fix16_t fix16_minimum  = 0x80000000; /*!< the minimum value of fix16_t */
static const fix16_t fix16_overflow = 0x80000000; /*!< the value used to indicate overflows when FIXMATH_NO_OVERFLOW is not specified */

static const fix16_t fix16_pi  = 205887;     /*!< fix16_t value of pi */
static const fix16_t fix16_e   = 178145;     /*!< fix16_t value of e */
static const fix16_t fix16_one = 0x00010000; /*!< fix16_t value of 1 */
static const fix16_t fix16_eps = 1;          /*!< fix16_t epsilon */

static const fix16_t fix16_rad_to_deg_mult = 3754936;
static const fix16_t fix16_deg_to_rad_mult = 1144;

/* Matrix error flags */
#define FIXMATRIX_OVERFLOW 0x01
#define FIXMATRIX_DIMERR   0x02
#define FIXMATRIX_USEERR   0x04
#define FIXMATRIX_SINGULAR 0x08
#define FIXMATRIX_NEGATIVE 0x10

/* Conversion macros */
#define F16(x) ((fix16_t)(((x) >= 0) ? ((x) * 65536.0 + 0.5) : ((x) * 65536.0 - 0.5)))

/* String conversion helper macros for F16C */
#define FIXMATH_TOKLEN(token) ( sizeof( #token ) - 1 )

#define FIXMATH_CONSTANT_POW10(times) ( \
  (times == 0) ? 1ULL \
        : (times == 1) ? 10ULL \
            : (times == 2) ? 100ULL \
                : (times == 3) ? 1000ULL \
                    : (times == 4) ? 10000ULL \
                        : (times == 5) ? 100000ULL \
                            : (times == 6) ? 1000000ULL \
                                : (times == 7) ? 10000000ULL \
                                    : 100000000ULL \
)

/** Helper macro for F16C, the type uint64_t is only used at compile time and
 *  shouldn't be visible in the generated code.
 *
 * @note We do not use fix16_one instead of 65536ULL, because the
 *       "use of a const variable in a constant expression is nonstandard in C".
 */
#define FIXMATH_CONVERT_MANTISSA(m) \
( (unsigned) \
    ( \
        ( \
            ( \
                (uint64_t)( ( ( 1 ## m ## ULL ) - FIXMATH_CONSTANT_POW10(FIXMATH_TOKLEN(m)) ) * FIXMATH_CONSTANT_POW10(5 - FIXMATH_TOKLEN(m)) ) \
                * 100000ULL * 65536ULL \
            ) \
            + 5000000000ULL /* rounding: + 0.5 */ \
        ) \
        / \
        10000000000LL \
    ) \
)

#define FIXMATH_COMBINE_I_M(i, m) \
( \
    ( \
        (    i ) \
        << 16 \
    ) \
    | \
    ( \
        FIXMATH_CONVERT_MANTISSA(m) \
        & 0xFFFF \
    ) \
)

/** Create int16_t (Q16.16) constant from separate integer and mantissa part.
 *
 * Only tested on 32-bit ARM Cortex-M0 / x86 Intel.
 *
 * This macro is needed when compiling with options like "--fpu=none",
 * which forbid all and every use of float and related types and
 * would thus make it impossible to have fix16_t constants.
 *
 * Just replace uses of F16() with F16C() like this:
 *   F16(123.1234) becomes F16C(123,1234)
 *
 * @warning Specification of any value outside the mentioned intervals
 *          WILL result in undefined behavior!
 *
 * @note Regardless of the specified minimum and maximum values for i and m below,
 *       the total value of the number represented by i and m MUST be in the interval
 *       ]-32768.00000:32767.99999[ else usage with this macro will yield undefined behavior.
 *
 * @param i Signed integer constant with a value in the interval ]-32768:32767[.
 * @param m Positive integer constant in the interval ]0:99999[ (fractional part/mantissa).
 */
#define F16C(i, m) \
( (fix16_t) \
    ( \
      (( #i[0] ) == '-') \
        ? -FIXMATH_COMBINE_I_M((unsigned)( ( (i) * -1) ), m) \
        : FIXMATH_COMBINE_I_M((unsigned)i, m) \
    ) \
)

/*===========================================================================*/
/* FUNCTION DECLARATIONS                                                     */
/*===========================================================================*/

/* Group A: Core Types & Basic Arithmetic */
FIXMATH_FUNC_ATTRS fix16_t fix16_from_int(int32_t a);
FIXMATH_FUNC_ATTRS float fix16_to_float(fix16_t a);
FIXMATH_FUNC_ATTRS double fix16_to_dbl(fix16_t a);
FIXMATH_FUNC_ATTRS int32_t fix16_to_int(fix16_t a);
FIXMATH_FUNC_ATTRS fix16_t fix16_from_float(float a);
FIXMATH_FUNC_ATTRS fix16_t fix16_from_dbl(double a);
FIXMATH_FUNC_ATTRS fix16_t fix16_abs(fix16_t x);
FIXMATH_FUNC_ATTRS fix16_t fix16_floor(fix16_t x);
FIXMATH_FUNC_ATTRS fix16_t fix16_ceil(fix16_t x);
FIXMATH_FUNC_ATTRS fix16_t fix16_min(fix16_t x, fix16_t y);
FIXMATH_FUNC_ATTRS fix16_t fix16_max(fix16_t x, fix16_t y);
FIXMATH_FUNC_ATTRS fix16_t fix16_clamp(fix16_t x, fix16_t lo, fix16_t hi);
FIXMATH_FUNC_ATTRS fix16_t fix16_rad_to_deg(fix16_t radians);
FIXMATH_FUNC_ATTRS fix16_t fix16_deg_to_rad(fix16_t degrees);
FIXMATH_FUNC_ATTRS fix16_t fix16_sq(fix16_t x);
FIXMATH_FUNC_ATTRS int32_t fix_abs(int32_t x);
FIXMATH_FUNC_ATTRS qf16 qf16_from_v3d(v3d axis, fix16_t angle);
FIXMATH_FUNC_ATTRS v3d qf16_to_v3d(qf16 qf);
FIXMATH_FUNC_ATTRS uint32_t uint32_log2(uint32_t x);
FIXMATH_FUNC_ATTRS fract32_t fract32_create(uint32_t a, uint32_t b);
FIXMATH_FUNC_ATTRS fract32_t fract32_invert(fract32_t a);
#ifndef FIXMATH_NO_64BIT
FIXMATH_FUNC_ATTRS uint32_t fract32_usmul(uint32_t a, fract32_t b);
FIXMATH_FUNC_ATTRS int32_t fract32_smul(int32_t a, fract32_t b);
#endif
FIXMATH_FUNC_ATTRS fix16_t fix16_add(fix16_t a, fix16_t b);
FIXMATH_FUNC_ATTRS fix16_t fix16_sub(fix16_t a, fix16_t b);
#ifndef FIXMATH_NO_OVERFLOW
FIXMATH_FUNC_ATTRS fix16_t fix16_sadd(fix16_t a, fix16_t b);
FIXMATH_FUNC_ATTRS fix16_t fix16_ssub(fix16_t a, fix16_t b);
FIXMATH_FUNC_ATTRS fix16_t fix16_smul(fix16_t a, fix16_t b);
#endif
FIXMATH_FUNC_ATTRS fix16_t fix16_mul(fix16_t a, fix16_t b);

/* Group B: Advanced Math Functions */
FIXMATH_FUNC_ATTRS fix16_t fix16_div(fix16_t a, fix16_t b);
	#ifndef FIXMATH_NO_OVERFLOW
FIXMATH_FUNC_ATTRS fix16_t fix16_sdiv(fix16_t a, fix16_t b);
	#endif
FIXMATH_FUNC_ATTRS fix16_t fix16_mod(fix16_t x, fix16_t y);
FIXMATH_FUNC_ATTRS fix16_t fix16_sqrt(fix16_t a);
FIXMATH_FUNC_ATTRS fix16_t fix16_lerp8(fix16_t a, fix16_t b, uint8_t t);
FIXMATH_FUNC_ATTRS fix16_t fix16_lerp16(fix16_t a, fix16_t b, uint16_t t);
FIXMATH_FUNC_ATTRS fix16_t fix16_lerp32(fix16_t a, fix16_t b, uint32_t t);
FIXMATH_FUNC_ATTRS fix16_t fix16_exp(fix16_t a);
FIXMATH_FUNC_ATTRS fix16_t fix16_log2(fix16_t x);
FIXMATH_FUNC_ATTRS fix16_t fix16_slog2(fix16_t x);
FIXMATH_FUNC_ATTRS fix16_t fix16_log(fix16_t x);
FIXMATH_FUNC_ATTRS fix16_t fix16_sin_parabola(fix16_t a);
FIXMATH_FUNC_ATTRS fix16_t fix16_sin(fix16_t a);
FIXMATH_FUNC_ATTRS fix16_t fix16_cos(fix16_t a);
FIXMATH_FUNC_ATTRS fix16_t fix16_tan(fix16_t a);
FIXMATH_FUNC_ATTRS fix16_t fix16_atan2(fix16_t y, fix16_t x);
FIXMATH_FUNC_ATTRS fix16_t fix16_atan(fix16_t a);
FIXMATH_FUNC_ATTRS fix16_t fix16_asin(fix16_t a);
FIXMATH_FUNC_ATTRS fix16_t fix16_acos(fix16_t a);

/* Group C: String & Utility Functions */
/*!
* Convert fix16_t value to a string.
* Required buffer length for largest values is 13 bytes.
*/
FIXMATH_FUNC_ATTRS void fix16_to_str(fix16_t value, char *buf, int decimals);
FIXMATH_FUNC_ATTRS fix16_t fix16_from_str(const char *str);

/* Group D: Vector Operations */
FIXMATH_FUNC_ATTRS void fa16_unalias(void *dest, void **a, void **b, void *tmp, unsigned size);
FIXMATH_FUNC_ATTRS fix16_t fa16_dot(const fix16_t *a, uint_fast8_t a_stride,
                 const fix16_t *b, uint_fast8_t b_stride,
                                   uint_fast8_t n);
FIXMATH_FUNC_ATTRS fix16_t fa16_norm(const fix16_t *a, uint_fast8_t a_stride, uint_fast8_t n);
FIXMATH_FUNC_ATTRS void v2d_add(v2d *dest, const v2d *a, const v2d *b);
FIXMATH_FUNC_ATTRS void v2d_sub(v2d *dest, const v2d *a, const v2d *b);
FIXMATH_FUNC_ATTRS void v2d_mul_s(v2d *dest, const v2d *a, fix16_t b);
FIXMATH_FUNC_ATTRS void v2d_div_s(v2d *dest, const v2d *a, fix16_t b);
FIXMATH_FUNC_ATTRS fix16_t v2d_dot(const v2d *a, const v2d *b);
FIXMATH_FUNC_ATTRS fix16_t v2d_norm(const v2d *a);
FIXMATH_FUNC_ATTRS void v2d_normalize(v2d *dest, const v2d *a);
FIXMATH_FUNC_ATTRS void v2d_rotate(v2d *dest, const v2d *a, fix16_t angle);
FIXMATH_FUNC_ATTRS void v3d_add(v3d *dest, const v3d *a, const v3d *b);
FIXMATH_FUNC_ATTRS void v3d_sub(v3d *dest, const v3d *a, const v3d *b);
FIXMATH_FUNC_ATTRS void v3d_mul_s(v3d *dest, const v3d *a, fix16_t b);
FIXMATH_FUNC_ATTRS void v3d_div_s(v3d *dest, const v3d *a, fix16_t b);
FIXMATH_FUNC_ATTRS fix16_t v3d_dot(const v3d *a, const v3d *b);
FIXMATH_FUNC_ATTRS fix16_t v3d_norm(const v3d *a);
FIXMATH_FUNC_ATTRS void v3d_normalize(v3d *dest, const v3d *a);
FIXMATH_FUNC_ATTRS void v3d_cross(v3d *dest, const v3d *a, const v3d *b);

/* Group E: Matrix Operations */
FIXMATH_FUNC_ATTRS void mf16_fill(mf16 *dest, fix16_t value);
FIXMATH_FUNC_ATTRS void mf16_transpose(mf16 *dest, const mf16 *matrix);
FIXMATH_FUNC_ATTRS void mf16_add(mf16 *dest, const mf16 *a, const mf16 *b);
FIXMATH_FUNC_ATTRS void mf16_sub(mf16 *dest, const mf16 *a, const mf16 *b);
FIXMATH_FUNC_ATTRS void mf16_mul(mf16 *dest, const mf16 *a, const mf16 *b);
FIXMATH_FUNC_ATTRS void mf16_mul_s(mf16 *dest, const mf16 *matrix, fix16_t scalar);
FIXMATH_FUNC_ATTRS void mf16_div_s(mf16 *dest, const mf16 *matrix, fix16_t scalar);
FIXMATH_FUNC_ATTRS void mf16_mul_bt(mf16 *dest, const mf16 *a, const mf16 *bt);
FIXMATH_FUNC_ATTRS void mf16_mul_at(mf16 *dest, const mf16 *at, const mf16 *b);
FIXMATH_FUNC_ATTRS void mf16_qr_decomposition(mf16 *q, mf16 *r, const mf16 *matrix, int reorthogonalize);
FIXMATH_FUNC_ATTRS void mf16_solve(mf16 *dest, const mf16 *q, const mf16 *r, const mf16 *matrix);
FIXMATH_FUNC_ATTRS void mf16_cholesky(mf16 *dest, const mf16 *matrix);
FIXMATH_FUNC_ATTRS void mf16_invert_lt(mf16 *dest, const mf16 *matrix);

/* Group F: Quaternion Operations */
FIXMATH_FUNC_ATTRS void qf16_conj(qf16 *dest, const qf16 *q);
FIXMATH_FUNC_ATTRS void qf16_add(qf16 *dest, const qf16 *q, const qf16 *r);
FIXMATH_FUNC_ATTRS void qf16_mul_s(qf16 *dest, const qf16 *q, fix16_t s);
FIXMATH_FUNC_ATTRS void qf16_div_s(qf16 *dest, const qf16 *q, fix16_t s);
FIXMATH_FUNC_ATTRS fix16_t qf16_dot(const qf16 *q, const qf16 *r);
FIXMATH_FUNC_ATTRS fix16_t qf16_norm(const qf16 *q);
FIXMATH_FUNC_ATTRS void qf16_normalize(qf16 *dest, const qf16 *q);
FIXMATH_FUNC_ATTRS void qf16_mul(qf16 *dest, const qf16 *q, const qf16 *r);
FIXMATH_FUNC_ATTRS void qf16_pow(qf16 *dest, const qf16 *q, fix16_t power);
FIXMATH_FUNC_ATTRS void qf16_avg(qf16 *dest, const qf16 *q1, const qf16 *q2, fix16_t weight);
FIXMATH_FUNC_ATTRS void qf16_from_axis_angle(qf16 *dest, const v3d *axis, fix16_t angle);
FIXMATH_FUNC_ATTRS void qf16_to_matrix(mf16 *dest, const qf16 *q);

/* Group G: Pretty Printing Functions */
#ifndef FIXSTRING_NO_STDIO
FIXMATH_FUNC_ATTRS void print_fix16_t(FILE *stream, fix16_t value, uint_fast8_t width, uint_fast8_t decimals);
FIXMATH_FUNC_ATTRS void print_mf16(FILE *stream, const mf16 *matrix);
#endif

/* Group H: FFT Functions */
#ifndef INPUT_TYPE
#define INPUT_TYPE uint8_t
#endif

#ifndef INPUT_CONVERT
#define INPUT_CONVERT(x) ((x) << 8)
#endif

#ifndef INPUT_INDEX
#define INPUT_INDEX(x) (x)
#endif

#ifndef OUTPUT_SCALE
#define OUTPUT_SCALE(transform_size)    (fix16_one * 256 / transform_size)
#endif

FIXMATH_FUNC_ATTRS void fix16_fft(INPUT_TYPE *input, fix16_t *real, fix16_t *imag, unsigned transform_length);

#endif /* LIBFIXMATHMATRIX_H */

/*===========================================================================*/
/* IMPLEMENTATION SECTION                                                    */
/*===========================================================================*/

#ifdef LIBFIXMATHMATRIX_IMPLEMENTATION

#include "implementation_group_a.h"
#include "implementation_group_b.h"
#include "implementation_group_c.h"
#include "implementation_group_d.h"
#include "implementation_group_e.h"
#include "implementation_group_f.h"
#include "implementation_group_g.h"
#include "implementation_group_h.h"

#endif /* LIBFIXMATHMATRIX_IMPLEMENTATION */

#ifdef __cplusplus
}
#endif

