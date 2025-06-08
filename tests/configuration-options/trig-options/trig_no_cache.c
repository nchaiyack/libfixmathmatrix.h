/*
 * No-cache trigonometric implementation
 * Configures libfixmath to disable caching for precise equivalence testing
 */

#define FIXMATH_NO_CACHE
#define LIBFIXMATHMATRIX_IMPLEMENTATION
#include <libfixmathmatrix_final.h>

/* Export no-cache implementation with prefixed names */
fix16_t fix16_sin_no_cache(fix16_t inAngle) {
    return fix16_sin(inAngle);
}

fix16_t fix16_cos_no_cache(fix16_t inAngle) {
    return fix16_cos(inAngle);
}

fix16_t fix16_tan_no_cache(fix16_t inAngle) {
    return fix16_tan(inAngle);
}

fix16_t fix16_sin_parabola_no_cache(fix16_t inAngle) {
    return fix16_sin_parabola(inAngle);
}
