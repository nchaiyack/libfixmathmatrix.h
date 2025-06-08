# Test Guidelines for 'libfixmathmatrix' Contributors New to Numerics

## Purpose and Scope

This document establishes testing standards for all test suites written for libfixmathmatrix. As a numerical library intended for graphics rendering and binary interchange formats, **we aspire to have good mathematical properties**, but aren't experts. Here are some guidelines that might help us write good tests, especially for those that don't touch numerics/the nitty gritty of numerical types often.

## Summary of Guidelines

1. **Implement comprehensive edge case testing** for all arithmetic operations
2. **Create reference test vector databases** for cross-platform validation
3. **Develop mathematical property verification** for all operations
4. **Establish error tolerance calculation framework** based on theoretical analysis
5. **Build performance validation suite** comparing optimized vs reference implementations

### Stretch Goals
1. Implement property-based testing framework
2. Create automated test generation tools
3. Develop numerical stability analysis tools
4. Build cross-compiler validation system
5. Implement statistical error analysis reporting

## Core Testing Principles

### Principle 1: Full Domain Coverage
**REQUIREMENT**: Test functions across their complete input domains, not artificially restricted subsets.

**MANDATORY PRACTICES**:
- Use full 32-bit ranges for int32_t inputs: `INT32_MIN` to `INT32_MAX`
- Test complete Q16.16 fixed-point range: `fix16_minimum` to `fix16_maximum`
- Never artificially restrict ranges "to avoid overflow" - test overflow behavior explicitly

**IMPLEMENTATION STANDARD**:
```c
// REQUIRED: Full range testing
int32_t test_value = random_int32();  // Full range

// FORBIDDEN: Artificial range restriction
int32_t test_value = random_int32() >> 8;  // Don't do this
```

### Principle 2: Mathematically Rigorous Tolerances
**REQUIREMENT**: Calculate error tolerances based on theoretical analysis, not arbitrary "close enough" values.

**TOLERANCE CALCULATION STANDARD**:
```c
// Base tolerance on theoretical precision limits
#define Q16_16_EPSILON (1)  // Smallest representable difference
#define SINGLE_OP_ERROR (Q16_16_EPSILON)
#define MULTI_OP_ERROR(n) ((n) * Q16_16_EPSILON * 2)  // 2x safety factor

// Example usage
assert(fix16_approx_equal(result, expected, SINGLE_OP_ERROR));
```

**FORBIDDEN PRACTICES**:
- Tolerances larger than 100 units without mathematical justification
- Using identical tolerances for operations with different error characteristics
- "Trial and error" tolerance adjustment until tests pass

### Principle 3: Comprehensive Edge Case Testing
**REQUIREMENT**: Every test suite MUST include systematic edge case coverage.

**MANDATORY EDGE CASES**:
1. **Boundary Values**: MIN, MAX, zero, ±1
2. **Sign Transitions**: Positive to negative, zero crossings
3. **Overflow/Underflow**: Deliberate out-of-range operations
4. **Error Conditions**: Division by zero, domain violations
5. **Precision Limits**: Values near representational boundaries

**IMPLEMENTATION TEMPLATE**:
```c
static void test_edge_cases_FUNCTION_NAME(void) {
    // Boundary values
    assert_function_behavior(fix16_minimum);
    assert_function_behavior(fix16_maximum);
    assert_function_behavior(0);
    assert_function_behavior(fix16_one);
    assert_function_behavior(-fix16_one);
    
    // Overflow testing
    fix16_t overflow_result = function_under_test(fix16_maximum, fix16_one);
    assert(overflow_result == fix16_overflow);  // Must handle overflow properly
    
    // Error condition testing
    fix16_t error_result = fix16_div(fix16_one, 0);
    assert(error_result == fix16_overflow);  // Must handle division by zero
}
```

### Principle 4: Mathematical Property Verification
**REQUIREMENT**: Test that mathematical identities and properties hold within precision limits.

**MANDATORY PROPERTY TESTS**:
1. **Algebraic Properties**: Commutativity, associativity, distributivity
2. **Identity Elements**: additive/multiplicative identities
3. **Inverse Relationships**: operation reversibility where applicable
4. **Mathematical Constants**: π, e, trigonometric identities to full precision

**PROPERTY TESTING TEMPLATE**:
```c
static void test_mathematical_properties(void) {
    fix16_t a = generate_test_value();
    fix16_t b = generate_test_value();
    fix16_t c = generate_test_value();
    
    // Commutativity: a + b = b + a
    fix16_t ab = fix16_add(a, b);
    fix16_t ba = fix16_add(b, a);
    assert(fix16_cmp_eq(ab, ba));
    
    // Associativity: (a + b) + c = a + (b + c)
    fix16_t left = fix16_add(fix16_add(a, b), c);
    fix16_t right = fix16_add(a, fix16_add(b, c));
    assert(fix16_approx_equal(left, right, MULTI_OP_ERROR(3)));
    
    // Inverse property: a - a = 0
    fix16_t zero_result = fix16_sub(a, a);
    assert(fix16_cmp_eq(zero_result, 0));
}
```

### Principle 5: Cross-Platform Consistency
**REQUIREMENT**: Ensure bit-exact reproducible results across all target platforms.

**CROSS-PLATFORM TESTING STANDARD**:
```c
typedef struct {
    int32_t input_a;
    int32_t input_b;
    int32_t expected_result;
    const char* operation_name;
} ReferenceTestVector;

// REQUIRED: Comprehensive reference test vectors
static const ReferenceTestVector reference_vectors[] = {
    {0x10000, 0x20000, 0x30000, "add_basic"},
    {0x7FFFFFFF, 0x1, fix16_overflow, "add_overflow"},
    {fix16_pi, fix16_e, 0x56A52, "mul_constants"},
    // ... Must include edge cases and typical values
};

static void test_cross_platform_consistency(void) {
    for (size_t i = 0; i < ARRAY_SIZE(reference_vectors); i++) {
        const ReferenceTestVector* tv = &reference_vectors[i];
        fix16_t result = operation_under_test(tv->input_a, tv->input_b);
        if (result != tv->expected_result) {
            printf("FAIL: %s - Expected 0x%X, got 0x%X\n", 
                   tv->operation_name, tv->expected_result, result);
            assert(0);
        }
    }
}
```

## General waffle about advice

### Input Generation Standards
**REQUIREMENT**: Use systematic and reproducible input generation.

```c
// REQUIRED: Seeded, reproducible random generation
static uint32_t test_seed = 0;

static void initialize_test_environment(void) {
    test_seed = (uint32_t)time(NULL);  // Or use fixed seed for reproducible tests
    printf("Test seed: %u (save this for test reproduction)\n", test_seed);
    srand(test_seed);
}

// REQUIRED: Systematic input space coverage
static void test_input_space_coverage(void) {
    // Test specific critical values
    test_critical_values();
    
    // Test systematic ranges
    for (int exp = -15; exp <= 15; exp++) {
        fix16_t base_value = fix16_from_int(1) << exp;
        test_function_at_value(base_value);
    }
    
    // Test random sampling
    for (int i = 0; i < REQUIRED_RANDOM_ITERATIONS; i++) {
        fix16_t random_value = generate_full_range_random();
        test_function_at_value(random_value);
    }
}
```

### Error Handling Verification
**REQUIREMENT**: Every test suite MUST verify proper error handling.

```c
static void test_error_handling(void) {
    // Division by zero
    assert(fix16_div(fix16_one, 0) == fix16_overflow);
    
    // Overflow conditions
    assert(fix16_add(fix16_maximum, fix16_one) == fix16_overflow);
    assert(fix16_sub(fix16_minimum, fix16_one) == fix16_overflow);
    
    // Domain violations (where applicable)
    assert(fix16_sqrt(-fix16_one) == fix16_overflow);
    
    // Result validation after error conditions
    // Library state should remain consistent after errors
}
```

### Performance Validation
**REQUIREMENT**: Verify that optimizations don't compromise correctness.

```c
static void test_performance_vs_correctness(void) {
    // Compare optimized vs reference implementations
    for (int i = 0; i < PERFORMANCE_TEST_ITERATIONS; i++) {
        fix16_t input = generate_test_value();
        
        fix16_t fast_result = optimized_function(input);
        fix16_t reference_result = reference_function(input);
        
        assert(fix16_approx_equal(fast_result, reference_result, 
                                 OPTIMIZATION_ERROR_TOLERANCE));
    }
}
```

### Documentation Requirements
**REQUIREMENT**: Every test function MUST include comprehensive documentation.

```c
/**
 * Test: fix16_add edge cases and mathematical properties
 * 
 * Coverage:
 * - Boundary values (MIN, MAX, zero, ±1)
 * - Overflow conditions
 * - Mathematical properties (commutativity, associativity)
 * - Sign handling
 * - Error propagation
 * 
 * Expected behavior:
 * - Commutative: a + b = b + a
 * - Associative within error bounds
 * - Proper overflow detection and handling
 * - Identity: a + 0 = a
 * 
 * Error tolerance: SINGLE_OP_ERROR (1 unit)
 * Rationale: Addition should be exact in fixed-point arithmetic
 */
static void test_fix16_add_comprehensive(void) {
    // Implementation here
}
```

## Test Suite Structure Requirements

### Mandatory Test Categories
Every function MUST have test suites covering:

1. **Basic Functionality**: Typical use cases with known results
2. **Edge Cases**: Boundary conditions and extreme values
3. **Error Conditions**: Invalid inputs and overflow scenarios  
4. **Mathematical Properties**: Algebraic identities and relationships
5. **Cross-Platform Consistency**: Reference test vectors
6. **Performance Validation**: Optimized vs reference comparison

### Test Organization Standard
```c
// File structure: tests/CATEGORY/FUNCTION_NAME_tests.c

int main(void) {
    initialize_test_environment();
    
    test_FUNCTION_basic_functionality();
    test_FUNCTION_edge_cases();
    test_FUNCTION_error_conditions();
    test_FUNCTION_mathematical_properties();
    test_FUNCTION_cross_platform_consistency();
    test_FUNCTION_performance_validation();
    
    print_test_summary();
    return get_test_failure_count();
}
```

### Test Metrics Requirements
**REQUIREMENT**: Every test suite MUST report quantitative metrics.

```c
typedef struct {
    unsigned int tests_executed;
    unsigned int tests_passed;
    unsigned int tests_failed;
    double max_observed_error;
    double avg_observed_error;
    unsigned int edge_cases_covered;
    unsigned int performance_samples;
} TestMetrics;

static void report_test_metrics(const TestMetrics* metrics) {
    printf("=== Test Metrics ===\n");
    printf("Tests executed: %u\n", metrics->tests_executed);
    printf("Pass rate: %.2f%%\n", 100.0 * metrics->tests_passed / metrics->tests_executed);
    printf("Max error observed: %.6f units\n", metrics->max_observed_error);
    printf("Edge cases covered: %u\n", metrics->edge_cases_covered);
    
    // REQUIREMENT: All metrics must meet minimum standards
    assert(metrics->tests_passed == metrics->tests_executed);  // 100% pass rate required
    assert(metrics->edge_cases_covered >= MINIMUM_EDGE_CASES);
    assert(metrics->max_observed_error <= MAXIMUM_ACCEPTABLE_ERROR);
}
```