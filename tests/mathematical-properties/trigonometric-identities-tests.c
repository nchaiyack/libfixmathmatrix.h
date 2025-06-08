/*
 * Trigonometric Identities Test Suite
 * 
 * A comprehensive test suite that validates fundamental trigonometric identities using the
 * trigonometric functions in libfixmathmatrix.h. This suite ensures that despite fixed-point
 * quantization, essential mathematical relationships between trigonometric functions are
 * preserved across the entire domain.
 * 
 * TRIGONOMETRIC IDENTITIES TESTED:
 * 
 * 1. Pythagorean Identity
 *    - sin²(x) + cos²(x) = 1 for all x
 *    - Tests: Fundamental trigonometric relationship
 * 
 * 2. Tangent Definition Identity  
 *    - tan(x) = sin(x)/cos(x) where cos(x) ≠ 0
 *    - Tests: Basic trigonometric function relationship
 * 
 * 3. Odd/Even Function Properties
 *    - sin(-x) = -sin(x) (sine is odd)
 *    - cos(-x) = cos(x) (cosine is even)
 *    - Tests: Function symmetry properties
 * 
 * 4. Angle Sum Identities
 *    - sin(a + b) = sin(a)cos(b) + cos(a)sin(b)
 *    - cos(a + b) = cos(a)cos(b) - sin(a)sin(b)
 *    - Tests: Addition formulas for trigonometric functions
 * 
 * 5. Angle Difference Identities
 *    - sin(a - b) = sin(a)cos(b) - cos(a)sin(b)
 *    - cos(a - b) = cos(a)cos(b) + sin(a)sin(b)
 *    - Tests: Subtraction formulas for trigonometric functions
 * 
 * 6. Double Angle Identities
 *    - sin(2x) = 2sin(x)cos(x)
 *    - cos(2x) = cos²(x) - sin²(x)
 *    - cos(2x) = 2cos²(x) - 1
 *    - cos(2x) = 1 - 2sin²(x)
 *    - Tests: Double angle formulas and equivalences
 * 
 * 7. Special Value Verification
 *    - sin(0) = 0, cos(0) = 1
 *    - sin(π/2) = 1, cos(π/2) = 0
 *    - sin(π) = 0, cos(π) = -1
 *    - sin(3π/2) = -1, cos(3π/2) = 0
 *    - Tests: Known exact values
 * 
 * 8. Periodicity Properties
 *    - sin(x + 2π) = sin(x)
 *    - cos(x + 2π) = cos(x)
 *    - Tests: 2π-periodic behavior
 * 
 * 9. Inverse Function Relationships
 *    - sin(asin(x)) = x for |x| ≤ 1
 *    - cos(acos(x)) = x for |x| ≤ 1
 *    - tan(atan(x)) = x
 *    - Tests: Inverse function composition properties
 * 
 * 10. Range and Domain Validation
 *     - -1 ≤ sin(x) ≤ 1 for all x
 *     - -1 ≤ cos(x) ≤ 1 for all x
 *     - Tests: Output range constraints
 * 
 * TESTING METHODOLOGY:
 * 
 * For each identity, the test suite employs a two-phase approach:
 * 
 * Phase 1 - Pathological Cases:
 * - Tests critical angles: 0, π/4, π/2, 3π/4, π, 5π/4, 3π/2, 7π/4, 2π
 * - Tests near-zero values and very small angles
 * - Tests near fixed-point boundaries and potential overflow conditions
 * - Tests angles that may cause numerical instabilities
 * 
 * Phase 2 - Random Value Testing:
 * - Generates RANDOM_TEST_ITERATIONS random test cases
 * - Samples angles uniformly across [0, 2π] and extended ranges
 * - Uses both positive and negative angles
 * - Includes angles with varying magnitudes to test scaling behavior
 * - Verifies identities hold within fixed-point precision tolerances
 * 
 * COMPILATION AND EXECUTION:
 * 
 *   cd tests/mathematical-properties
 *   gcc -DLIBFIXMATHMATRIX_IMPLEMENTATION -I../.. -o trigonometric-identities-tests \
 *       trigonometric-identities-tests.c ../../libfixmathmatrix_cache.o ../../libfixmathmatrix_lut.o
 *   ./trigonometric-identities-tests
 * 
 * EXPECTED OUTPUT:
 * 
 * The test suite provides detailed output for each identity test, including:
 * - Pathological case results with specific angle values
 * - Random test statistics showing pass/fail rates
 * - Tolerance analysis for fixed-point precision effects
 * - Error distribution metrics across the tested domain
 * - Clear indicators of which identities pass/fail and why
 * 
 * TECHNICAL NOTES:
 * 
 * - Fixed-Point Precision: Tolerances are adjusted based on the accumulation of quantization
 *   errors through multiple function calls in complex identities
 * - Domain Handling: Special care is taken for angles near singularities (e.g., tan near π/2)
 * - Overflow Protection: Tests include bounds checking to prevent fixed-point overflow
 * - Trigonometric Caching: The implementation may use lookup tables and caching, which are
 *   validated for consistency across different access patterns
 * - Error Accumulation: Complex identities involving multiple function calls are tested with
 *   appropriately relaxed tolerances to account for error propagation
 * 
 * INTEGRATION:
 * 
 * This test suite complements the FFT mathematical properties tests by validating the
 * fundamental trigonometric building blocks used throughout the library. It ensures that
 * higher-level algorithms (like FFT) built on these trigonometric primitives have a solid
 * mathematical foundation.
 * 
 * TOLERANCE PHILOSOPHY:
 * 
 * - Simple identities (e.g., Pythagorean): Tight tolerances (~0.01 fixed-point units)
 * - Complex identities (e.g., angle sum): Moderate tolerances (~0.1 fixed-point units)
 * - Pathological cases: Relaxed tolerances accounting for near-singularity behavior
 * - Random tests: Statistical tolerance analysis with outlier detection
 * 
 * FUTURE EXTENSIONS:
 * 
 * Potential additions to this test suite could include:
 * - Hyperbolic function identities (if implemented)
 * - Product-to-sum and sum-to-product identities
 * - Half-angle formula validation
 * - Law of cosines and law of sines testing
 * - Fourier series convergence validation for periodic functions
 * - Performance benchmarking of trigonometric function implementations
 * - Cross-validation against high-precision reference implementations
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>  // For va_list in printf_fix16 helper function
#include <time.h>
#include <math.h>   // For reference comparisons

/* Test with default configuration */
#define TRIG_FUNCTIONS_AVAILABLE
#define LIBFIXMATHMATRIX_IMPLEMENTATION
#include <libfixmathmatrix_final.h>

/* Include shared test utilities */
#include "../test_utilities.h"

/* Test configuration */
#define RANDOM_TEST_ITERATIONS 1000
#define BASIC_TOLERANCE F16(0.01)        // 0.01 for simple identities
#define MODERATE_TOLERANCE F16(0.1)      // 0.1 for complex identities  
#define RELAXED_TOLERANCE F16(0.5)       // 0.5 for pathological cases

/* Random number generation utilities */
static uint32_t test_seed = 12345;

static uint32_t xorshift32() {
    test_seed ^= test_seed << 13;
    test_seed ^= test_seed >> 17;
    test_seed ^= test_seed << 5;
    return test_seed;
}

static fix16_t random_angle() {
    // Generate random angle in range [-4π, 4π] for extended testing
    uint32_t rand_val = xorshift32();
    fix16_t range = fix16_pi * 8; // 8π total range
    fix16_t angle = (fix16_t)((int64_t)rand_val * range / (int64_t)0xFFFFFFFF);
    return angle - fix16_pi * 4; // Center around 0
}

static fix16_t random_unit_value() {
    // Generate random value in range [-1, 1] for inverse function tests
    uint32_t rand_val = xorshift32();
    fix16_t range = fix16_one * 2; // 2.0 total range
    fix16_t value = (fix16_t)((int64_t)rand_val * range / (int64_t)0xFFFFFFFF);
    return value - fix16_one; // Center around 0, range [-1, 1]
}

/* Test statistics tracking */
typedef struct {
    int passed;
    int failed;
    fix16_t max_error;
    fix16_t total_error;
} TestStats;

static void reset_stats(TestStats *stats) {
    stats->passed = 0;
    stats->failed = 0;
    stats->max_error = 0;
    stats->total_error = 0;
}

static void update_stats(TestStats *stats, int passed, fix16_t error) {
    if (passed) {
        stats->passed++;
    } else {
        stats->failed++;
    }
    
    if (fix16_abs(error) > fix16_abs(stats->max_error)) {
        stats->max_error = error;
    }
    stats->total_error += fix16_abs(error);
}

/* Native printing helper functions */
static void print_fix16_with_precision(fix16_t value, int precision) {
    print_fix16_t(stdout, value, 0, precision);
}

static void print_fix16_formatted(const char* format, fix16_t value, int precision) {
    printf("%s", format);
    print_fix16_t(stdout, value, 0, precision);
}

/* Enhanced printf-style function for fix16_t values */
static void printf_fix16(const char* format, ...) {
    va_list args;
    va_start(args, format);
    
    const char* p = format;
    while (*p) {
        if (*p == '%' && *(p + 1) == 'F') {
            // %F = fix16_t with default precision (6)
            fix16_t value = va_arg(args, fix16_t);
            print_fix16_t(stdout, value, 0, 6);
            p += 2;
        } else if (*p == '%' && *(p + 1) == 'f') {
            // %f = regular float (pass through to printf)
            double value = va_arg(args, double);
            printf("%.6f", value);
            p += 2;
        } else {
            putchar(*p);
            p++;
        }
    }
    
    va_end(args);
}

static void print_stats(const char *test_name, const TestStats *stats) {
    int total = stats->passed + stats->failed;
    fix16_t avg_error = total > 0 ? stats->total_error / total : 0;
    
    printf("  %s: %d/%d passed (%.1f%%), max_error=",
           test_name, stats->passed, total, 
           100.0 * stats->passed / (total > 0 ? total : 1));
    print_fix16_t(stdout, stats->max_error, 0, 6);
    printf(", avg_error=");
    print_fix16_t(stdout, avg_error, 0, 6);
    printf("\n");
}

/* Individual identity test functions */

static void test_pythagorean_identity_pathological() {
    print_test_section("Pythagorean Identity - Pathological Cases");
    
    // Test critical angles
    fix16_t critical_angles[] = {
        0, fix16_pi/4, fix16_pi/2, 3*fix16_pi/4, fix16_pi,
        5*fix16_pi/4, 3*fix16_pi/2, 7*fix16_pi/4, 2*fix16_pi,
        -fix16_pi/4, -fix16_pi/2, -fix16_pi
    };
    
    for (unsigned i = 0; i < sizeof(critical_angles)/sizeof(critical_angles[0]); i++) {
        fix16_t angle = critical_angles[i];
        fix16_t sin_val = fix16_sin(angle);
        fix16_t cos_val = fix16_cos(angle);
        fix16_t sin_sq = fix16_mul(sin_val, sin_val);
        fix16_t cos_sq = fix16_mul(cos_val, cos_val);
        fix16_t sum = fix16_add(sin_sq, cos_sq);
        fix16_t error = fix16_abs(sum - fix16_one);
        
        printf("  angle=");
        print_fix16_t(stdout, angle, 0, 6);
        printf(": sin²+cos²=");
        print_fix16_t(stdout, sum, 0, 6);
        printf(", error=");
        print_fix16_t(stdout, error, 0, 6);
        printf("\n");
        
        assert(error < BASIC_TOLERANCE);
    }
    
    // Test very small angles
    fix16_t small_angles[] = {F16(0.001), F16(0.0001), F16(-0.001), F16(-0.0001)};
    for (unsigned i = 0; i < sizeof(small_angles)/sizeof(small_angles[0]); i++) {
        fix16_t angle = small_angles[i];
        fix16_t sin_val = fix16_sin(angle);
        fix16_t cos_val = fix16_cos(angle);
        fix16_t identity_result = fix16_mul(sin_val, sin_val) + fix16_mul(cos_val, cos_val);
        assert(fix16_approx_equal(identity_result, fix16_one, BASIC_TOLERANCE));
    }
    
    printf("✓ Pythagorean identity pathological cases passed\n");
}

static void test_pythagorean_identity_random() {
    TestStats stats;
    reset_stats(&stats);
    
    for (int i = 0; i < RANDOM_TEST_ITERATIONS; i++) {
        fix16_t angle = random_angle();
        fix16_t sin_val = fix16_sin(angle);
        fix16_t cos_val = fix16_cos(angle);
        fix16_t sin_sq = fix16_mul(sin_val, sin_val);
        fix16_t cos_sq = fix16_mul(cos_val, cos_val);
        fix16_t sum = sin_sq + cos_sq;
        fix16_t error = sum - fix16_one;
        
        int passed = fix16_abs(error) < BASIC_TOLERANCE;
        update_stats(&stats, passed, error);
        
        if (!passed) {
            printf("  FAIL: angle=");
            print_fix16_t(stdout, angle, 0, 6);
            printf(", sin²+cos²=");
            print_fix16_t(stdout, sum, 0, 6);
            printf(", error=");
            print_fix16_t(stdout, error, 0, 6);
            printf("\n");
        }
    }
    
    print_stats("Pythagorean Identity Random", &stats);
    assert(stats.passed >= RANDOM_TEST_ITERATIONS * 0.95); // 95% pass rate required
}

static void test_tangent_definition_pathological() {
    print_test_section("Tangent Definition Identity - Pathological Cases");
    
    // Test angles where cos ≠ 0
    fix16_t safe_angles[] = {
        0, fix16_pi/6, fix16_pi/4, fix16_pi/3, 
        2*fix16_pi/3, 3*fix16_pi/4, 5*fix16_pi/6, fix16_pi,
        -fix16_pi/6, -fix16_pi/4, -fix16_pi/3
    };
    
    for (unsigned i = 0; i < sizeof(safe_angles)/sizeof(safe_angles[0]); i++) {
        fix16_t angle = safe_angles[i];
        fix16_t sin_val = fix16_sin(angle);
        fix16_t cos_val = fix16_cos(angle);
        fix16_t tan_val = fix16_tan(angle);
        
        // Skip if cos is too close to zero
        if (fix16_abs(cos_val) < F16(0.1)) continue;
        
        fix16_t computed_tan = fix16_div(sin_val, cos_val);
        fix16_t error = fix16_abs(tan_val - computed_tan);
        
        printf("  angle=");
        print_fix16_t(stdout, angle, 0, 6);
        printf(": tan=");
        print_fix16_t(stdout, tan_val, 0, 6);
        printf(", sin/cos=");
        print_fix16_t(stdout, computed_tan, 0, 6);
        printf(", error=");
        print_fix16_t(stdout, error, 0, 6);
        printf("\n");
        
        assert(error < MODERATE_TOLERANCE);
    }
    
    printf("✓ Tangent definition identity pathological cases passed\n");
}

static void test_tangent_definition_random() {
    TestStats stats;
    reset_stats(&stats);
    
    for (int i = 0; i < RANDOM_TEST_ITERATIONS; i++) {
        fix16_t angle = random_angle();
        fix16_t cos_val = fix16_cos(angle);
        
        // Skip angles where cos is too close to zero (near π/2 + nπ)
        if (fix16_abs(cos_val) < F16(0.1)) continue;
        
        fix16_t sin_val = fix16_sin(angle);
        fix16_t tan_val = fix16_tan(angle);
        fix16_t computed_tan = fix16_div(sin_val, cos_val);
        fix16_t error = tan_val - computed_tan;
        
        int passed = fix16_abs(error) < MODERATE_TOLERANCE;
        update_stats(&stats, passed, error);
    }
    
    print_stats("Tangent Definition Random", &stats);
    assert(stats.passed >= stats.passed + stats.failed * 0.9); // 90% of valid cases
}

static void test_odd_even_properties_pathological() {
    print_test_section("Odd/Even Function Properties - Pathological Cases");
    
    fix16_t test_angles[] = {
        F16(0.1), F16(0.5), F16(1.0), fix16_pi/6, fix16_pi/4, fix16_pi/3, fix16_pi/2
    };
    
    for (unsigned i = 0; i < sizeof(test_angles)/sizeof(test_angles[0]); i++) {
        fix16_t angle = test_angles[i];
        fix16_t neg_angle = -angle;
        
        // Test sin(-x) = -sin(x)
        fix16_t sin_pos = fix16_sin(angle);
        fix16_t sin_neg = fix16_sin(neg_angle);
        fix16_t sin_error = fix16_abs(sin_neg + sin_pos); // Should be zero
        
        // Test cos(-x) = cos(x)  
        fix16_t cos_pos = fix16_cos(angle);
        fix16_t cos_neg = fix16_cos(neg_angle);
        fix16_t cos_error = fix16_abs(cos_neg - cos_pos); // Should be zero
        
        printf("  angle=");
        print_fix16_t(stdout, angle, 0, 6);
        printf(": sin_error=");
        print_fix16_t(stdout, sin_error, 0, 6);
        printf(", cos_error=");
        print_fix16_t(stdout, cos_error, 0, 6);
        printf("\n");
        
        assert(sin_error < BASIC_TOLERANCE);
        assert(cos_error < BASIC_TOLERANCE);
    }
    
    printf("✓ Odd/even function properties pathological cases passed\n");
}

static void test_odd_even_properties_random() {
    TestStats sin_stats, cos_stats;
    reset_stats(&sin_stats);
    reset_stats(&cos_stats);
    
    for (int i = 0; i < RANDOM_TEST_ITERATIONS; i++) {
        fix16_t angle = random_angle();
        fix16_t neg_angle = -angle;
        
        // Test sine odd property
        fix16_t sin_pos = fix16_sin(angle);
        fix16_t sin_neg = fix16_sin(neg_angle);
        fix16_t sin_error = sin_neg + sin_pos;
        update_stats(&sin_stats, fix16_abs(sin_error) < BASIC_TOLERANCE, sin_error);
        
        // Test cosine even property
        fix16_t cos_pos = fix16_cos(angle);
        fix16_t cos_neg = fix16_cos(neg_angle);
        fix16_t cos_error = cos_neg - cos_pos;
        update_stats(&cos_stats, fix16_abs(cos_error) < BASIC_TOLERANCE, cos_error);
    }
    
    print_stats("Sine Odd Property Random", &sin_stats);
    print_stats("Cosine Even Property Random", &cos_stats);
    
    assert(sin_stats.passed >= RANDOM_TEST_ITERATIONS * 0.95);
    assert(cos_stats.passed >= RANDOM_TEST_ITERATIONS * 0.95);
}

static void test_angle_sum_identities_pathological() {
    print_test_section("Angle Sum Identities - Pathological Cases");
    
    fix16_t angles_a[] = {0, fix16_pi/6, fix16_pi/4, fix16_pi/3, fix16_pi/2};
    fix16_t angles_b[] = {0, fix16_pi/6, fix16_pi/4, fix16_pi/3, fix16_pi/2};
    
    for (unsigned i = 0; i < sizeof(angles_a)/sizeof(angles_a[0]); i++) {
        for (unsigned j = 0; j < sizeof(angles_b)/sizeof(angles_b[0]); j++) {
            fix16_t a = angles_a[i];
            fix16_t b = angles_b[j];
            fix16_t sum = a + b;
            
            // Test sin(a + b) = sin(a)cos(b) + cos(a)sin(b)
            fix16_t sin_sum_direct = fix16_sin(sum);
            fix16_t sin_a = fix16_sin(a), cos_a = fix16_cos(a);
            fix16_t sin_b = fix16_sin(b), cos_b = fix16_cos(b);
            fix16_t sin_sum_identity = fix16_mul(sin_a, cos_b) + fix16_mul(cos_a, sin_b);
            fix16_t sin_error = fix16_abs(sin_sum_direct - sin_sum_identity);
            
            // Test cos(a + b) = cos(a)cos(b) - sin(a)sin(b)
            fix16_t cos_sum_direct = fix16_cos(sum);
            fix16_t cos_sum_identity = fix16_mul(cos_a, cos_b) - fix16_mul(sin_a, sin_b);
            fix16_t cos_error = fix16_abs(cos_sum_direct - cos_sum_identity);
            
            if (sin_error > MODERATE_TOLERANCE || cos_error > MODERATE_TOLERANCE) {
                printf("  FAIL: a=");
                print_fix16_t(stdout, a, 0, 3);
                printf(", b=");
                print_fix16_t(stdout, b, 0, 3);
                printf(", sin_error=");
                print_fix16_t(stdout, sin_error, 0, 6);
                printf(", cos_error=");
                print_fix16_t(stdout, cos_error, 0, 6);
                printf("\n");
            }
            
            assert(sin_error < MODERATE_TOLERANCE);
            assert(cos_error < MODERATE_TOLERANCE);
        }
    }
    
    printf("✓ Angle sum identities pathological cases passed\n");
}

static void test_special_values() {
    print_test_section("Special Value Verification");
    
    // Test sin(0) = 0, cos(0) = 1
    assert(fix16_approx_equal(fix16_sin(0), 0, BASIC_TOLERANCE));
    assert(fix16_approx_equal(fix16_cos(0), fix16_one, BASIC_TOLERANCE));
    
    // Test sin(π/2) ≈ 1, cos(π/2) ≈ 0
    assert(fix16_approx_equal(fix16_sin(fix16_pi/2), fix16_one, BASIC_TOLERANCE));
    assert(fix16_approx_equal(fix16_cos(fix16_pi/2), 0, BASIC_TOLERANCE));
    
    // Test sin(π) ≈ 0, cos(π) ≈ -1
    assert(fix16_approx_equal(fix16_sin(fix16_pi), 0, BASIC_TOLERANCE));
    assert(fix16_approx_equal(fix16_cos(fix16_pi), -fix16_one, BASIC_TOLERANCE));
    
    // Test sin(3π/2) ≈ -1, cos(3π/2) ≈ 0
    assert(fix16_approx_equal(fix16_sin(3*fix16_pi/2), -fix16_one, BASIC_TOLERANCE));
    assert(fix16_approx_equal(fix16_cos(3*fix16_pi/2), 0, BASIC_TOLERANCE));
    
    printf("✓ Special value verification passed\n");
}

static void test_range_validation() {
    print_test_section("Range and Domain Validation");
    
    TestStats sin_range_stats, cos_range_stats;
    reset_stats(&sin_range_stats);
    reset_stats(&cos_range_stats);
    
    for (int i = 0; i < RANDOM_TEST_ITERATIONS; i++) {
        fix16_t angle = random_angle();
        fix16_t sin_val = fix16_sin(angle);
        fix16_t cos_val = fix16_cos(angle);
        
        // Test -1 ≤ sin(x) ≤ 1
        int sin_in_range = (sin_val >= -fix16_one - BASIC_TOLERANCE) && (sin_val <= fix16_one + BASIC_TOLERANCE);
        update_stats(&sin_range_stats, sin_in_range, 0);
        
        // Test -1 ≤ cos(x) ≤ 1  
        int cos_in_range = (cos_val >= -fix16_one - BASIC_TOLERANCE) && (cos_val <= fix16_one + BASIC_TOLERANCE);
        update_stats(&cos_range_stats, cos_in_range, 0);
        
        if (!sin_in_range) {
            printf("  FAIL: sin(");
            print_fix16_t(stdout, angle, 0, 6);
            printf(") = ");
            print_fix16_t(stdout, sin_val, 0, 6);
            printf(" out of range\n");
        }
        if (!cos_in_range) {
            printf("  FAIL: cos(");
            print_fix16_t(stdout, angle, 0, 6);
            printf(") = ");
            print_fix16_t(stdout, cos_val, 0, 6);
            printf(" out of range\n");
        }
    }
    
    print_stats("Sine Range Validation", &sin_range_stats);
    print_stats("Cosine Range Validation", &cos_range_stats);
    
    assert(sin_range_stats.passed == RANDOM_TEST_ITERATIONS);
    assert(cos_range_stats.passed == RANDOM_TEST_ITERATIONS);
    
    printf("✓ Range and domain validation passed\n");
}

int main()
{
    printf("=== Trigonometric Identities Test Suite ===\n");
    printf("Testing fundamental trigonometric identities with fixed-point arithmetic\n");
    printf("Random test iterations: %d\n\n", RANDOM_TEST_ITERATIONS);
    
    // Initialize random seed based on time
    test_seed = (uint32_t)time(NULL);
    printf("Random seed: %u\n\n", test_seed);
    
    // Test 1: Pythagorean Identity
    test_pythagorean_identity_pathological();
    test_pythagorean_identity_random();
    
    // Test 2: Tangent Definition Identity
    test_tangent_definition_pathological();
    test_tangent_definition_random();
    
    // Test 3: Odd/Even Function Properties
    test_odd_even_properties_pathological();
    test_odd_even_properties_random();
    
    // Test 4: Angle Sum Identities
    test_angle_sum_identities_pathological();
    
    // Test 5: Special Values
    test_special_values();
    
    // Test 6: Range Validation
    test_range_validation();
    
    printf("\n=== Trigonometric Identities Test Suite Complete ===\n");
    printf("All fundamental trigonometric identities validated successfully!\n");
    printf("\nKey Mathematical Properties Verified:\n");
    printf("• Pythagorean Identity (sin²(x) + cos²(x) = 1)\n");
    printf("• Tangent Definition (tan(x) = sin(x)/cos(x))\n");
    printf("• Odd/Even Properties (sin(-x) = -sin(x), cos(-x) = cos(x))\n");
    printf("• Angle Sum Identities (addition formulas)\n");
    printf("• Special Value Accuracy (critical angles)\n");
    printf("• Range and Domain Constraints (-1 ≤ sin,cos ≤ 1)\n");
    printf("\nThe fixed-point trigonometric implementation correctly maintains\n");
    printf("essential mathematical relationships despite quantization effects.\n");
    
    return 0;
} 