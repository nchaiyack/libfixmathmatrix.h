/*
 * example_usage_test.c - Example of using shared test utilities
 * 
 * This demonstrates how other test families can leverage the shared
 * test_utilities.h header for consistent testing patterns.
 */

#include <stdio.h>

/* Include the library being tested */
#define LIBFIXMATHMATRIX_IMPLEMENTATION
#include <libfixmathmatrix_final.h>

/* Include shared test utilities */
#include "test_utilities.h"

/* Example test function demonstrating various utility usage patterns */
int test_matrix_operations() {
    int tests_passed = 0;
    int total_tests = 0;
    
    print_test_section("Matrix Operations Example");
    
    // Create test matrices
    mf16 matrix_a, matrix_b, result;
    matrix_a.rows = matrix_a.columns = 2;
    matrix_b.rows = matrix_b.columns = 2;
    matrix_a.errors = matrix_b.errors = 0;
    
    // Fill matrices with test data
    matrix_a.data[0][0] = F16(1.0); matrix_a.data[0][1] = F16(2.0);
    matrix_a.data[1][0] = F16(3.0); matrix_a.data[1][1] = F16(4.0);
    
    matrix_b.data[0][0] = F16(2.0); matrix_b.data[0][1] = F16(0.0);
    matrix_b.data[1][0] = F16(1.0); matrix_b.data[1][1] = F16(2.0);
    
    // Test 1: Matrix addition
    total_tests++;
    mf16_add(&result, &matrix_a, &matrix_b);
    
    // Create expected result
    mf16 expected_add;
    expected_add.rows = expected_add.columns = 2;
    expected_add.errors = 0;
    expected_add.data[0][0] = F16(3.0); expected_add.data[0][1] = F16(2.0);
    expected_add.data[1][0] = F16(4.0); expected_add.data[1][1] = F16(6.0);
    
    if (mf16_approx_equal(&result, &expected_add, 100) && mf16_no_errors(&result)) {
        printf("✅ Matrix addition test passed\n");
        tests_passed++;
    } else {
        printf("❌ Matrix addition test failed\n");
    }
    
    // Test 2: Matrix properties
    total_tests++;
    if (mf16_is_square(&matrix_a) && mf16_same_dimensions(&matrix_a, &matrix_b)) {
        printf("✅ Matrix properties test passed\n");
        tests_passed++;
    } else {
        printf("❌ Matrix properties test failed\n");
    }
    
    return tests_passed == total_tests;
}

/* Example test using enhanced assertion macros */
int test_enhanced_assertions() {
    print_test_section("Enhanced Assertions Example");
    
    // Test various assertion patterns
    fix16_t pi_approx = F16(3.14159);
    
    // Using convenience macros
    ASSERT_FIX16_EQUAL(pi_approx, fix16_pi);
    ASSERT_FIX16_APPROX(fix16_sqrt(F16(16.0)), F16(4.0), 100);
    
    // String testing
    char buffer[32];
    fix16_to_str(fix16_e, buffer, 3);
    ASSERT_STR_PREFIX(buffer, "2.71");
    
    // Vector testing
    v3d vec1 = {F16(1.0), F16(2.0), F16(3.0)};
    v3d vec2 = {F16(1.0), F16(2.0), F16(3.0)};
    ASSERT_V3D_APPROX(&vec1, &vec2, 10);
    
    printf("✅ All enhanced assertion tests passed\n");
    return 1;
}

/* Example test using percentage-based tolerances */
int test_percentage_tolerances() {
    print_test_section("Percentage Tolerance Example");
    
    fix16_t calculated_e = fix16_exp(fix16_one); // e^1 should equal e
    
    // Test with 1% tolerance
    if (fix16_approx_equal_percent(calculated_e, fix16_e, 1.0)) {
        printf("✅ e^1 ≈ e within 1%% tolerance\n");
    } else {
        printf("❌ e^1 calculation outside 1%% tolerance\n");
        return 0;
    }
    
    // Test with decimal comparison
    if (fix16_approx_equal_decimal(fix16_pi, 3.14159, 0.001)) {
        printf("✅ π ≈ 3.14159 within decimal tolerance\n");
    } else {
        printf("❌ π constant outside decimal tolerance\n");
        return 0;
    }
    
    return 1;
}

/* Example test demonstrating vector utilities */
int test_vector_utilities() {
    print_test_section("Vector Utilities Example");
    
    // Test 3D vector operations
    v3d unit_x = {F16(1.0), F16(0.0), F16(0.0)};
    v3d unit_y = {F16(0.0), F16(1.0), F16(0.0)};
    v3d cross_result;
    
    v3d_cross(&cross_result, &unit_x, &unit_y);
    
    // Should get unit_z vector (0, 0, 1)
    if (v3d_approx_equal_decimal(&cross_result, 0.0, 0.0, 1.0, 0.01)) {
        printf("✅ Cross product i × j = k\n");
    } else {
        printf("❌ Cross product calculation failed\n");
        return 0;
    }
    
    // Test 2D vector operations
    v2d vec_a = {F16(3.0), F16(4.0)};
    v2d vec_b = {F16(-4.0), F16(3.0)};
    
    // These should be perpendicular (dot product = 0)
    fix16_t dot_product = v2d_dot(&vec_a, &vec_b);
    if (fix16_approx_equal(dot_product, F16(0.0), 100)) {
        printf("✅ Perpendicular vectors have zero dot product\n");
    } else {
        printf("❌ Perpendicular vector test failed\n");
        return 0;
    }
    
    return 1;
}

/* Main function demonstrating complete test suite */
int main() {
    printf("=== Shared Test Utilities Example ===\n");
    printf("Demonstrating usage patterns for test_utilities.h\n");
    
    int total_suites = 0;
    int passed_suites = 0;
    
    // Run test suites
    total_suites++; if (test_matrix_operations()) passed_suites++;
    total_suites++; if (test_enhanced_assertions()) passed_suites++;
    total_suites++; if (test_percentage_tolerances()) passed_suites++;
    total_suites++; if (test_vector_utilities()) passed_suites++;
    
    // Print final summary
    print_test_summary("Example Test Suite", passed_suites, total_suites);
    
    printf("\n🎯 Key Utilities Demonstrated:\n");
    printf("   • fix16_approx_equal() - Fixed-point comparison with tolerance\n");
    printf("   • fix16_approx_equal_decimal() - Fixed-point to decimal comparison\n");
    printf("   • fix16_approx_equal_percent() - Percentage-based tolerance\n");
    printf("   • v2d_approx_equal(), v3d_approx_equal() - Vector comparisons\n");
    printf("   • v3d_approx_equal_decimal() - Vector to decimal comparison\n");
    printf("   • mf16_approx_equal() - Matrix comparisons\n");
    printf("   • mf16_no_errors(), mf16_is_square() - Matrix property checks\n");
    printf("   • qf16_approx_equal() - Quaternion comparisons\n");
    printf("   • ASSERT_* macros - Enhanced assertions with messages\n");
    printf("   • print_test_section(), print_test_summary() - Formatted output\n");
    
    return (passed_suites == total_suites) ? 0 : 1;
} 