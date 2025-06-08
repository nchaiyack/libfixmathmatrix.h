/*
 * Test program for libfixmathmatrix.h
 * 
 * This tests the basic functionality of the combined header-only library
 * with configurable inlining.
 */

#include <stdio.h>
#include <stdarg.h>  // For enhanced printing functions

/* Test with default configuration */
#define TRIG_FUNCTIONS_AVAILABLE
#define LIBFIXMATHMATRIX_IMPLEMENTATION
#include <libfixmathmatrix_final.h>

/* Include shared test utilities */
#include "../test_utilities.h"

int main()
{
    printf("Testing libfixmathmatrix...\n");
    
    /* Test basic constants */
    printf("fix16_pi = ");
    print_fix16_t(stdout, fix16_pi, 0, 6);
    printf(" (should be ~3.141593)\n");
        assert(fix16_pi == 205887);
    printf("fix16_e = ");
    print_fix16_t(stdout, fix16_e, 0, 6);
    printf(" (should be ~2.718282)\n");
        assert(fix16_e == 178145);
    printf("fix16_one = ");
    print_fix16_t(stdout, fix16_one, 0, 6);
    printf(" (should be 1.000000)\n");
        assert(fix16_one == 65536);
    
    /* Test basic arithmetic */
    fix16_t a = F16(3.14159);
    fix16_t b = F16(2.71828);
    
    printf("a = ");
    print_fix16_t(stdout, a, 0, 6);
    printf(" (F16(3.14159))\n");
        assert(fix16_approx_equal(a, 205887, 10)); // π with small tolerance
    printf("b = ");
    print_fix16_t(stdout, b, 0, 6);
    printf(" (F16(2.71828))\n");
        assert(fix16_approx_equal(b, 178145, 10)); // e with small tolerance
    
    /* Test uint32 function */
    uint32_t log_val = uint32_log2(1024);
    printf("uint32_log2(1024) = %u (should be 10)\n", log_val);
        assert(log_val == 10);
    
    /* Test fract32 functions */
    fract32_t fract = fract32_create(1, 3);
    printf("fract32_create(1, 3) = %u\n", fract);
        assert(fract != 0); // Should create a valid fraction
    
    fract32_t inv_fract = fract32_invert(fract);
    printf("fract32_invert() = %u\n", inv_fract);
        assert(inv_fract != 0); // Should create a valid inverted fraction
    
    /* Test fix16 arithmetic */
    fix16_t sum = fix16_add(a, b);
    printf("fix16_add(pi, e) = ");
    print_fix16_t(stdout, sum, 0, 6);
    printf("\n");
        assert(fix16_approx_equal(sum, F16(5.86), 1000)); // π + e ≈ 5.86
    
    fix16_t product = fix16_mul(a, b);
    printf("fix16_mul(pi, e) = ");
    print_fix16_t(stdout, product, 0, 6);
    printf("\n");
        assert(fix16_approx_equal(product, F16(8.54), 2000)); // π × e ≈ 8.54
    
    /* Test new division function */
    fix16_t quotient = fix16_div(a, b);
    printf("fix16_div(pi, e) = ");
    print_fix16_t(stdout, quotient, 0, 6);
    printf("\n");
        assert(fix16_approx_equal(quotient, F16(1.156), 1000)); // π / e ≈ 1.156
    
    /* Test square root */
    fix16_t sqrt_val = fix16_sqrt(F16(4.0));
    printf("fix16_sqrt(4.0) = ");
    print_fix16_t(stdout, sqrt_val, 0, 6);
    printf(" (should be ~2.000000)\n");
        assert(fix16_approx_equal(sqrt_val, F16(2.0), 100)); // √4 = 2

    /* Test interpolation */
    fix16_t lerp_val = fix16_lerp16(F16(0.0), F16(10.0), 32768); // 50%
    printf("fix16_lerp16(0, 10, 50%%) = ");
    print_fix16_t(stdout, lerp_val, 0, 6);
    printf(" (should be ~5.000000)\n");
        assert(fix16_approx_equal(lerp_val, F16(5.0), 100)); // 50% between 0 and 10 = 5

    /* Test trigonometric functions (when available) */
    #ifdef TRIG_FUNCTIONS_AVAILABLE
    fix16_t sin_val = fix16_sin(fix16_pi / 2); // sin(π/2) = 1
    printf("fix16_sin(π/2) = ");
    print_fix16_t(stdout, sin_val, 0, 6);
    printf(" (should be ~1.000000)\n");
        assert(fix16_approx_equal(sin_val, fix16_one, 1000)); // sin(π/2) = 1

    fix16_t cos_val = fix16_cos(0); // cos(0) = 1  
    printf("fix16_cos(0) = ");
    print_fix16_t(stdout, cos_val, 0, 6);
    printf(" (should be ~1.000000)\n");
        assert(fix16_approx_equal(cos_val, fix16_one, 1000)); // cos(0) = 1
        
    fix16_t atan_val = fix16_atan(fix16_one); // atan(1) = π/4
    printf("fix16_atan(1) = ");
    print_fix16_t(stdout, atan_val, 0, 6);
    printf(" (should be ~0.785398 = π/4)\n");
        assert(fix16_approx_equal(atan_val, fix16_pi / 4, 1000)); // atan(1) = π/4
    #endif
    
    /* Test exponential and logarithmic functions */
    fix16_t exp_val = fix16_exp(fix16_one); // e^1 = e
    printf("fix16_exp(1) = ");
    print_fix16_t(stdout, exp_val, 0, 6);
    printf(" (should be ~2.718282 = e)\n");
        assert(fix16_approx_equal(exp_val, fix16_e, 2000)); // e^1 = e

    fix16_t ln_val = fix16_log(fix16_e); // ln(e) = 1
    printf("fix16_log(e) = ");
    print_fix16_t(stdout, ln_val, 0, 6);
    printf(" (should be ~1.000000)\n");
        assert(fix16_approx_equal(ln_val, fix16_one, 2000)); // ln(e) = 1

    fix16_t log2_val = fix16_log2(F16(8.0)); // log2(8) = 3
    printf("fix16_log2(8) = ");
    print_fix16_t(stdout, log2_val, 0, 6);
    printf(" (should be ~3.000000)\n");
        assert(fix16_approx_equal(log2_val, F16(3.0), 1000)); // log2(8) = 3

    fix16_t slog2_val = fix16_slog2(F16(4.0)); // log2(4) = 2
    printf("fix16_slog2(4) = ");
    print_fix16_t(stdout, slog2_val, 0, 6);
    printf(" (should be ~2.000000)\n");
        assert(fix16_approx_equal(slog2_val, F16(2.0), 1000)); // log2(4) = 2

    /* Test string conversion functions */
    char str_buffer[32];
    fix16_to_str(fix16_pi, str_buffer, 4);
    printf("fix16_to_str(π, 4 decimals) = \"%s\" (should be ~3.1416)\n", str_buffer);
        assert(str_approx_equal(str_buffer, "3.141")); // π ≈ 3.1416

    fix16_to_str(fix16_e, str_buffer, 3);
    printf("fix16_to_str(e, 3 decimals) = \"%s\" (should be ~2.718)\n", str_buffer);
        assert(str_approx_equal(str_buffer, "2.71")); // e ≈ 2.718

    fix16_to_str(F16(-1.5), str_buffer, 2);
    printf("fix16_to_str(-1.5, 2 decimals) = \"%s\" (should be \"-1.50\")\n", str_buffer);
        assert(str_approx_equal(str_buffer, "-1.5")); // -1.5

    fix16_t parsed_val = fix16_from_str("3.14159");
    printf("fix16_from_str(\"3.14159\") = ");
    print_fix16_t(stdout, parsed_val, 0, 6);
    printf(" (should be ~3.141593 = π)\n");
        assert(fix16_approx_equal(parsed_val, fix16_pi, 100)); // Parse π

    fix16_t parsed_neg = fix16_from_str("-2.718");
    printf("fix16_from_str(\"-2.718\") = ");
    print_fix16_t(stdout, parsed_neg, 0, 6);
    printf(" (should be ~-2.718000 = -e)\n");
        assert(fix16_approx_equal(parsed_neg, -fix16_e, 1000)); // Parse -e

    fix16_t parsed_int = fix16_from_str("42");
    printf("fix16_from_str(\"42\") = ");
    print_fix16_t(stdout, parsed_int, 0, 6);
    printf(" (should be ~42.000000)\n");
        assert(fix16_approx_equal(parsed_int, F16(42.0), 10)); // Parse 42

    /* FFT tests have been moved to dedicated test suite */
    print_test_section("FFT Functions");
    printf("FFT tests have been extracted to tests/mathematical-properties/fft-tests.c\n");
    printf("Run that dedicated test suite for comprehensive FFT mathematical property validation.\n");
    printf("The FFT implementation includes:\n");
    printf("  • DC Component Preservation\n");
    printf("  • Impulse Response (Flat Spectrum)\n");
    printf("  • Linearity Property\n");
    printf("  • Conjugate Symmetry (Hermitian Property)\n");
    printf("  • Known Frequency Response\n");
    printf("  • Energy Conservation (Parseval's Theorem)\n");
    printf("  • Multiple Transform Lengths (Power-of-2)\n");
    printf("  • Phase Relationship Validation\n");
    
    /* Test array utility functions */
    print_test_section("Array Utilities Test");
    
    // Test vector dot product
    fix16_t vec_a[4] = {F16(1.0), F16(2.0), F16(3.0), F16(4.0)};
    fix16_t vec_b[4] = {F16(2.0), F16(3.0), F16(4.0), F16(5.0)};
    
    fix16_t dot_result = fa16_dot(vec_a, 1, vec_b, 1, 4);
    printf("fa16_dot([1,2,3,4] • [2,3,4,5]) = ");
    print_fix16_t(stdout, dot_result, 0, 6);
    printf(" (should be ~40.000000)\n");
        assert(fix16_approx_equal(dot_result, F16(40.0), 1000)); // 1×2 + 2×3 + 3×4 + 4×5 = 40

    // Test vector norm
    fix16_t norm_vec[3] = {F16(3.0), F16(4.0), F16(0.0)};
    fix16_t norm_result = fa16_norm(norm_vec, 1, 3);
    printf("fa16_norm([3,4,0]) = ");
    print_fix16_t(stdout, norm_result, 0, 6);
    printf(" (should be ~5.000000)\n");
        assert(fix16_approx_equal(norm_result, F16(5.0), 1000)); // √(3² + 4² + 0²) = 5

    // Test with negative values
    fix16_t neg_vec[2] = {F16(-3.0), F16(4.0)};
    fix16_t neg_norm = fa16_norm(neg_vec, 1, 2);
    printf("fa16_norm([-3,4]) = ");
    print_fix16_t(stdout, neg_norm, 0, 6);
    printf(" (should be ~5.000000)\n");
    assert(fix16_approx_equal(neg_norm, F16(5.0), 1000)); // √((-3)² + 4²) = 5

    /* Test 2D vector operations */
    print_test_section("2D Vector Operations Test");
    
    v2d vec2a = {F16(3.0), F16(4.0)};
    v2d vec2b = {F16(1.0), F16(2.0)};
    v2d vec2_result;
    
    // Test 2D vector addition
    v2d_add(&vec2_result, &vec2a, &vec2b);
    printf("v2d_add([3,4] + [1,2]) = [");
    print_fix16_t(stdout, vec2_result.x, 0, 4);
    printf(", ");
    print_fix16_t(stdout, vec2_result.y, 0, 4);
    printf("] (should be [4,6])\n");
        v2d expected_add = {F16(4.0), F16(6.0)};
        assert(v2d_approx_equal(&vec2_result, &expected_add, 10)); // Using new utility function

    // Test 2D vector subtraction  
    v2d_sub(&vec2_result, &vec2a, &vec2b);
    printf("v2d_sub([3,4] - [1,2]) = [");
    print_fix16_t(stdout, vec2_result.x, 0, 4);
    printf(", ");
    print_fix16_t(stdout, vec2_result.y, 0, 4);
    printf("] (should be [2,2])\n");
        v2d expected_sub = {F16(2.0), F16(2.0)};
        assert(v2d_approx_equal(&vec2_result, &expected_sub, 10));

    // Test 2D scalar multiplication
    v2d_mul_s(&vec2_result, &vec2a, F16(2.0));
    printf("v2d_mul_s([3,4] * 2) = [");
    print_fix16_t(stdout, vec2_result.x, 0, 4);
    printf(", ");
    print_fix16_t(stdout, vec2_result.y, 0, 4);
    printf("] (should be [6,8])\n");
        v2d expected_mul = {F16(6.0), F16(8.0)};
        assert(v2d_approx_equal(&vec2_result, &expected_mul, 10));

    // Test 2D dot product
    fix16_t dot2d = v2d_dot(&vec2a, &vec2b);
    printf("v2d_dot([3,4] • [1,2]) = ");
    print_fix16_t(stdout, dot2d, 0, 6);
    printf(" (should be ~11.000000)\n");
        assert(fix16_approx_equal(dot2d, F16(11.0), 100)); // 3×1 + 4×2 = 11

    // Test 2D norm
    fix16_t norm2d = v2d_norm(&vec2a);
    printf("v2d_norm([3,4]) = ");
    print_fix16_t(stdout, norm2d, 0, 6);
    printf(" (should be ~5.000000)\n");
        assert(fix16_approx_equal(norm2d, F16(5.0), 1000)); // √(3² + 4²) = 5

    /* Test 3D vector operations */
    print_test_section("3D Vector Operations Test");
    
    v3d vec3a = {F16(1.0), F16(2.0), F16(3.0)};
    v3d vec3b = {F16(4.0), F16(5.0), F16(6.0)};
    v3d vec3_result;
    
    // Test 3D vector addition
    v3d_add(&vec3_result, &vec3a, &vec3b);
    printf("v3d_add([1,2,3] + [4,5,6]) = [");
    print_fix16_t(stdout, vec3_result.x, 0, 4);
    printf(", ");
    print_fix16_t(stdout, vec3_result.y, 0, 4);
    printf(", ");
    print_fix16_t(stdout, vec3_result.z, 0, 4);
    printf("] (should be [5,7,9])\n");
        v3d expected_add3d = {F16(5.0), F16(7.0), F16(9.0)};
        assert(v3d_approx_equal(&vec3_result, &expected_add3d, 10)); // Using new utility function

    // Test 3D vector subtraction
    v3d_sub(&vec3_result, &vec3b, &vec3a);
    printf("v3d_sub([4,5,6] - [1,2,3]) = [");
    print_fix16_t(stdout, vec3_result.x, 0, 4);
    printf(", ");
    print_fix16_t(stdout, vec3_result.y, 0, 4);
    printf(", ");
    print_fix16_t(stdout, vec3_result.z, 0, 4);
    printf("] (should be [3,3,3])\n");
        v3d expected_sub3d = {F16(3.0), F16(3.0), F16(3.0)};
        assert(v3d_approx_equal(&vec3_result, &expected_sub3d, 10));

    // Test 3D dot product
    fix16_t dot3d = v3d_dot(&vec3a, &vec3b);
    printf("v3d_dot([1,2,3] • [4,5,6]) = ");
    print_fix16_t(stdout, dot3d, 0, 6);
    printf(" (should be ~32.000000)\n");
        assert(fix16_approx_equal(dot3d, F16(32.0), 100)); // 1×4 + 2×5 + 3×6 = 32

    // Test 3D cross product
    v3d cross_a = {F16(1.0), F16(0.0), F16(0.0)};
    v3d cross_b = {F16(0.0), F16(1.0), F16(0.0)};
    v3d_cross(&vec3_result, &cross_a, &cross_b);
    printf("v3d_cross([1,0,0] × [0,1,0]) = [");
    print_fix16_t(stdout, vec3_result.x, 0, 4);
    printf(", ");
    print_fix16_t(stdout, vec3_result.y, 0, 4);
    printf(", ");
    print_fix16_t(stdout, vec3_result.z, 0, 4);
    printf("] (should be [0,0,1])\n");
        // Using decimal comparison utility
        assert(v3d_approx_equal_decimal(&vec3_result, 0.0, 0.0, 1.0, 0.01));

    // Test 3D norm  
    fix16_t norm3d = v3d_norm(&vec3a);
    printf("v3d_norm([1,2,3]) = ");
    print_fix16_t(stdout, norm3d, 0, 6);
    printf(" (should be ~3.741657)\n");
    assert(fix16_approx_equal_decimal(norm3d, 3.742, 0.01)); // Using decimal comparison

    /* Test conversions */
    printf("fix16_to_float(pi) = %f\n", fix16_to_float(a));
    printf("fix16_to_int(pi) = %d\n", fix16_to_int(a));
        float float_pi = fix16_to_float(a);
        int int_pi = fix16_to_int(a);
        assert(fabs(float_pi - 3.14159) < 0.01); // Should be approximately π
        assert(int_pi == 3); // Integer part of π
    
    /* Test type definitions */
    v2d vec2 = {F16(1.0), F16(2.0)};
    v3d vec3 = {F16(1.0), F16(2.0), F16(3.0)};
    printf("v2d: (");
    print_fix16_t(stdout, vec2.x, 0, 4);
    printf(", ");
    print_fix16_t(stdout, vec2.y, 0, 4);
    printf(")\n");
    printf("v3d: (");
    print_fix16_t(stdout, vec3.x, 0, 4);
    printf(", ");
    print_fix16_t(stdout, vec3.y, 0, 4);
    printf(", ");
    print_fix16_t(stdout, vec3.z, 0, 4);
    printf(")\n");
        assert(fix16_approx_equal(vec2.x, F16(1.0), 10));
        assert(fix16_approx_equal(vec2.y, F16(2.0), 10));
        assert(fix16_approx_equal(vec3.x, F16(1.0), 10));
        assert(fix16_approx_equal(vec3.y, F16(2.0), 10));
        assert(fix16_approx_equal(vec3.z, F16(3.0), 10));

    mf16 matrix;
    matrix.rows = 2;
    matrix.columns = 2;
    matrix.errors = 0;
    printf("Matrix size: %dx%d\n", matrix.rows, matrix.columns);
        assert(matrix.rows == 2);
        assert(matrix.columns == 2);
        assert(mf16_no_errors(&matrix)); // Using new utility function
        assert(mf16_is_square(&matrix));  // Using new utility function

    qf16 quat = {F16(1.0), F16(0.0), F16(0.0), F16(0.0)};
    printf("Quaternion: (");
    print_fix16_t(stdout, quat.a, 0, 4);
    printf(", ");
    print_fix16_t(stdout, quat.b, 0, 4);
    printf(", ");
    print_fix16_t(stdout, quat.c, 0, 4);
    printf(", ");
    print_fix16_t(stdout, quat.d, 0, 4);
    printf(")\n");
        qf16 expected_quat = {F16(1.0), F16(0.0), F16(0.0), F16(0.0)};
        assert(qf16_approx_equal(&quat, &expected_quat, 10)); // Using new utility function
        
    printf("All tests completed successfully!\n");
    return 0;
} 