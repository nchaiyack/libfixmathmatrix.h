/*
 * libfixmathmatrix_cache.c - Performance caches for libfixmathmatrix
 * 
 * This file contains the performance caches used by trigonometric and
 * exponential functions. By placing them in a separate compilation unit,
 * we maintain the header-only design while providing optional performance
 * optimizations.
 * 
 * Memory usage when FIXMATH_NO_CACHE is not defined:
 * - Trigonometric caches: ~64KB (sin + atan caches)
 * - Exponential caches: ~32KB (exp cache)
 * - Total: ~96KB
 * 
 * On memory-constrained platforms, you can:
 * - Not link this file (disables all caches)
 * - Define FIXMATH_NO_CACHE (disables cache usage)
 * - Place this file in a specific memory section
 */

#include <stdint.h>

/* Forward declaration of fix16_t */
typedef int32_t fix16_t;

/*===========================================================================*/
/* TRIGONOMETRIC FUNCTION CACHES                                            */
/*===========================================================================*/

#ifndef FIXMATH_NO_CACHE

/* Sine function cache (32KB when enabled) */
fix16_t _fix16_sin_cache_index[4096] = { 0 };
fix16_t _fix16_sin_cache_value[4096] = { 0 };

/* Arctangent function caches (48KB when enabled) */
fix16_t _fix16_atan_cache_index[2][4096] = { { 0 }, { 0 } };
fix16_t _fix16_atan_cache_value[4096] = { 0 };

/* Exponential function cache (32KB when enabled) */
fix16_t _fix16_exp_cache_index[4096] = { 0 };
fix16_t _fix16_exp_cache_value[4096] = { 0 };

#endif /* FIXMATH_NO_CACHE */ 