/*
 * Lookup table trigonometric implementation
 * FIXMATH_SIN_LUT: Uses lookup table instead of polynomial
 */

#define FIXMATH_SIN_LUT
#define LIBFIXMATHMATRIX_STATIC
#define LIBFIXMATHMATRIX_IMPLEMENTATION
#include <libfixmathmatrix_final.h>

/* Export LUT implementation with prefixed names */
fix16_t fix16_sin_lut(fix16_t inAngle) {
    return fix16_sin(inAngle);
}

fix16_t fix16_cos_lut(fix16_t inAngle) {
    return fix16_cos(inAngle);
}

fix16_t fix16_tan_lut(fix16_t inAngle) {
    return fix16_tan(inAngle);
}

fix16_t fix16_sin_parabola_lut(fix16_t inAngle) {
    return fix16_sin_parabola(inAngle);
}
