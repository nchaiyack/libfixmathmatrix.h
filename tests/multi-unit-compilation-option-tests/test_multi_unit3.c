/* Test file 3: Unit that tests trigonometric functions */
#include "libfixmathmatrix_final.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

/* Forward declarations for address comparison functions */
extern void* get_fix16_sin_address(void);

/* Test function addresses to verify single instantiation */
void* get_fix16_sin_address_unit3(void) {
    return (void*)fix16_sin;
}

int test_trigonometry_from_unit3(void) {
    printf("=== Testing Trigonometry (Unit 3) ===\n");
    
    /* Test that we're using the same function instances */
    void* addr1 = get_fix16_sin_address();
    void* addr2 = get_fix16_sin_address_unit3();
    printf("fix16_sin address from unit 1: %p\n", addr1);
    printf("fix16_sin address from unit 3: %p\n", addr2);
    assert(addr1 == addr2);
    printf("✓ Function addresses match - single instantiation confirmed\n");
    
    /* Test trigonometric functions */
    fix16_t zero = fix16_from_int(0);
    fix16_t pi_half = fix16_pi >> 1;  /* PI/2 */
    fix16_t pi = fix16_pi;
    
    /* Test sin */
    fix16_t sin_0 = fix16_sin(zero);
    fix16_t sin_pi_2 = fix16_sin(pi_half);
    fix16_t sin_pi = fix16_sin(pi);
    
    printf("sin(0) = %f (expected: ~0)\n", fix16_to_float(sin_0));
    printf("sin(π/2) = %f (expected: ~1)\n", fix16_to_float(sin_pi_2));
    printf("sin(π) = %f (expected: ~0)\n", fix16_to_float(sin_pi));
    /* <add asserts here> */

    /* Test cos */
    fix16_t cos_0 = fix16_cos(zero);
    fix16_t cos_pi_2 = fix16_cos(pi_half);
    fix16_t cos_pi = fix16_cos(pi);
    
    printf("cos(0) = %f (expected: ~1)\n", fix16_to_float(cos_0));
    printf("cos(π/2) = %f (expected: ~0)\n", fix16_to_float(cos_pi_2));
    printf("cos(π) = %f (expected: ~-1)\n", fix16_to_float(cos_pi));
    /* <add asserts here> */

    /* Test atan2 */
    fix16_t atan2_result = fix16_atan2(fix16_one, fix16_one);
    printf("atan2(1, 1) = %f (expected: ~0.785 = π/4)\n", fix16_to_float(atan2_result));
    /* <add asserts here> */

    /* Test degree/radian conversion */
    fix16_t deg_90 = fix16_from_int(90);
    fix16_t rad_from_deg = fix16_deg_to_rad(deg_90);
    fix16_t deg_back = fix16_rad_to_deg(rad_from_deg);
    
    printf("90° -> radians -> degrees: %f -> %f -> %f\n", 
           fix16_to_float(deg_90), fix16_to_float(rad_from_deg), fix16_to_float(deg_back));
    /* <add asserts here> */

    /* Test exp and log functions */
    fix16_t exp_1 = fix16_exp(fix16_one);
    printf("exp(1) = %f (expected: ~2.718)\n", fix16_to_float(exp_1));
    /* <add asserts here> */

    fix16_t log_e = fix16_log(exp_1);
    printf("log(exp(1)) = %f (expected: ~1)\n", fix16_to_float(log_e));
    /* <add asserts here> */

    printf("Unit 3 tests: PASSED\n\n");
    return 1;
} 