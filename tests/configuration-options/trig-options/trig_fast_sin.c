/*
 * Fast sine trigonometric implementation
 * FIXMATH_FAST_SIN: Uses faster but less accurate algorithm
 */

#define FIXMATH_FAST_SIN
#define LIBFIXMATHMATRIX_STATIC
#define LIBFIXMATHMATRIX_IMPLEMENTATION
#include <libfixmathmatrix_final.h>

/* Export fast implementation with prefixed names */
fix16_t fix16_sin_fast(fix16_t inAngle) {
    return fix16_sin(inAngle);
}

fix16_t fix16_cos_fast(fix16_t inAngle) {
    return fix16_cos(inAngle);
}

fix16_t fix16_tan_fast(fix16_t inAngle) {
    return fix16_tan(inAngle);
}

fix16_t fix16_sin_parabola_fast(fix16_t inAngle) {
    return fix16_sin_parabola(inAngle);
}
