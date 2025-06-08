/* Test file 3: Static instantiation - Unit 3 */
#define LIBFIXMATHMATRIX_STATIC
#define LIBFIXMATHMATRIX_IMPLEMENTATION
#include "libfixmathmatrix_final.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

/* Test function addresses - these should be different from other units */
void* get_static_fix16_sin_address_unit3(void) {
    return (void*)fix16_sin;
}

void* get_static_v2d_add_address_unit3(void) {
    return (void*)v2d_add;
}

int test_static_trigonometry_from_unit3(void) {
    printf("=== Testing Static Trigonometry & Vectors (Unit 3) ===\n");
    
    /* Test trigonometric functions - these are static to this unit */
    fix16_t zero = fix16_from_int(0);
    fix16_t pi_quarter = fix16_pi >> 2;  /* PI/4 */
    fix16_t pi_half = fix16_pi >> 1;     /* PI/2 */
    
    /* Test sin and cos with static functions */
    fix16_t sin_0 = fix16_sin(zero);
    fix16_t sin_pi_4 = fix16_sin(pi_quarter);
    fix16_t sin_pi_2 = fix16_sin(pi_half);
    
    fix16_t cos_0 = fix16_cos(zero);
    fix16_t cos_pi_4 = fix16_cos(pi_quarter);
    fix16_t cos_pi_2 = fix16_cos(pi_half);
    
    printf("sin(0) = %f, sin(π/4) = %f, sin(π/2) = %f\n",
           fix16_to_float(sin_0), fix16_to_float(sin_pi_4), fix16_to_float(sin_pi_2));
    printf("cos(0) = %f, cos(π/4) = %f, cos(π/2) = %f\n",
           fix16_to_float(cos_0), fix16_to_float(cos_pi_4), fix16_to_float(cos_pi_2));
        /* <add asserts here> */
        
    /* Test atan2 and atan */
    fix16_t atan2_result = fix16_atan2(fix16_one, fix16_one);
    fix16_t atan_result = fix16_atan(fix16_one);
    
    printf("atan2(1, 1) = %f, atan(1) = %f\n",
           fix16_to_float(atan2_result), fix16_to_float(atan_result));
       /* <add asserts here> */

    /* Test degree/radian conversion with static functions */
    fix16_t deg_45 = fix16_from_int(45);
    fix16_t rad_from_deg = fix16_deg_to_rad(deg_45);
    fix16_t deg_back = fix16_rad_to_deg(rad_from_deg);
    
    printf("45° -> radians -> degrees: %f -> %f -> %f\n",
           fix16_to_float(deg_45), fix16_to_float(rad_from_deg), fix16_to_float(deg_back));
       /* <add asserts here> */

    /* Test exp/log functions */
    fix16_t exp_2 = fix16_exp(fix16_from_int(2));
    fix16_t log_exp2 = fix16_log(exp_2);
    fix16_t log2_8 = fix16_log2(fix16_from_int(8));
    
    printf("exp(2) = %f, log(exp(2)) = %f, log2(8) = %f\n",
           fix16_to_float(exp_2), fix16_to_float(log_exp2), fix16_to_float(log2_8));
       /* <add asserts here> */

    /* Test 2D vector operations with static functions */
    v2d vec_a = {fix16_from_int(6), fix16_from_int(8)};
    v2d vec_b = {fix16_from_int(3), fix16_from_int(4)};
    v2d vec_result;
    
    v2d_add(&vec_result, &vec_a, &vec_b);
    printf("v2d_add([6,8], [3,4]) = [%d,%d] (expected: [9,12])\n",
           fix16_to_int(vec_result.x), fix16_to_int(vec_result.y));
    assert(fix16_to_int(vec_result.x) == 9);
    assert(fix16_to_int(vec_result.y) == 12);
    
    v2d_mul_s(&vec_result, &vec_a, fix16_from_int(2));
    printf("v2d_mul_s([6,8], 2) = [%d,%d] (expected: [12,16])\n",
           fix16_to_int(vec_result.x), fix16_to_int(vec_result.y));
        /* <add asserts here> */

    fix16_t dot_product = v2d_dot(&vec_a, &vec_b);
    printf("v2d_dot([6,8], [3,4]) = %d (expected: 50)\n", fix16_to_int(dot_product));
    assert(fix16_to_int(dot_product) == 50);
    
    fix16_t norm = v2d_norm(&vec_a);
    printf("v2d_norm([6,8]) = %f (expected: 10.0)\n", fix16_to_float(norm));
    `    /* <add asserts here> */

    /* Test 3D vector cross product */
    v3d vec3_a = {fix16_one, fix16_from_int(2), fix16_from_int(3)};
    v3d vec3_b = {fix16_from_int(2), fix16_one, fix16_from_int(2)};
    v3d vec3_result;
    
    v3d_cross(&vec3_result, &vec3_a, &vec3_b);
    printf("v3d_cross([1,2,3], [2,1,2]) = [%d,%d,%d]\n",
           fix16_to_int(vec3_result.x), fix16_to_int(vec3_result.y), fix16_to_int(vec3_result.z));
       /* <add asserts here> */

    /* Test matrix operations */
    mf16 matrix = {{2, 2, 0, {{fix16_from_int(2), fix16_from_int(3)},
                               {fix16_from_int(1), fix16_from_int(4)}}}};
    mf16 matrix_result;
    
    mf16_transpose(&matrix_result, &matrix);
    printf("Matrix transpose result:\n");
    printf("  [%d %d]\n", fix16_to_int(matrix_result.data[0][0]), fix16_to_int(matrix_result.data[0][1]));
    printf("  [%d %d]\n", fix16_to_int(matrix_result.data[1][0]), fix16_to_int(matrix_result.data[1][1]));
       /* <add asserts here> */

    /* Print function addresses for comparison */
    printf("Unit 3 fix16_sin address: %p (static)\n", get_static_fix16_sin_address_unit3());
    printf("Unit 3 v2d_add address: %p (static)\n", get_static_v2d_add_address_unit3());
       /* <add asserts here> */

    printf("Static Unit 3 tests: PASSED\n\n");
    return 1;
} 