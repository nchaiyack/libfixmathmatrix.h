/*
 * Trigonometric function variants for configuration testing
 * Each variant represents a different compilation configuration
 */

#ifndef TRIG_VARIANTS_H
#define TRIG_VARIANTS_H

/* Configuration-specific implementations of trigonometric functions */
#include <libfixmathmatrix_final.h>

/* Standard implementation (default: cached, polynomial, standard accuracy) */
fix16_t fix16_sin_standard(fix16_t inAngle);
fix16_t fix16_cos_standard(fix16_t inAngle);
fix16_t fix16_tan_standard(fix16_t inAngle);
fix16_t fix16_sin_parabola_standard(fix16_t inAngle);

/* No-cache implementation (FIXMATH_NO_CACHE) */
fix16_t fix16_sin_no_cache(fix16_t inAngle);
fix16_t fix16_cos_no_cache(fix16_t inAngle);
fix16_t fix16_tan_no_cache(fix16_t inAngle);
fix16_t fix16_sin_parabola_no_cache(fix16_t inAngle);

/* Lookup table implementation (FIXMATH_SIN_LUT) */
fix16_t fix16_sin_lut(fix16_t inAngle);
fix16_t fix16_cos_lut(fix16_t inAngle);
fix16_t fix16_tan_lut(fix16_t inAngle);
fix16_t fix16_sin_parabola_lut(fix16_t inAngle);

/* Fast implementation (FIXMATH_FAST_SIN) */
fix16_t fix16_sin_fast(fix16_t inAngle);
fix16_t fix16_cos_fast(fix16_t inAngle);
fix16_t fix16_tan_fast(fix16_t inAngle);
fix16_t fix16_sin_parabola_fast(fix16_t inAngle);

#endif /* TRIG_VARIANTS_H */
