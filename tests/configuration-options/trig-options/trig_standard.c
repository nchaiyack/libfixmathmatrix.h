/*
 * Standard trigonometric implementation 
 * Default configuration: cached, polynomial, standard accuracy
 */

#define LIBFIXMATHMATRIX_STATIC
#define LIBFIXMATHMATRIX_IMPLEMENTATION
#include <libfixmathmatrix_final.h>

/* Export standard implementation with prefixed names */
fix16_t fix16_sin_standard(fix16_t inAngle) {
    return fix16_sin(inAngle);
}

fix16_t fix16_cos_standard(fix16_t inAngle) {
    return fix16_cos(inAngle);
}

fix16_t fix16_tan_standard(fix16_t inAngle) {
    return fix16_tan(inAngle);
}

fix16_t fix16_sin_parabola_standard(fix16_t inAngle) {
    return fix16_sin_parabola(inAngle);
}
