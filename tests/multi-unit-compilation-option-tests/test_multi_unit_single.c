/* Test file 1: Main compilation unit with implementation */
#define LIBFIXMATHMATRIX_IMPLEMENTATION
#include "libfixmathmatrix_final.h"
#include <stdio.h>
#include <assert.h>

/* Forward declarations for functions from other units */
extern int test_basic_math_from_unit2(void);
extern int test_trigonometry_from_unit3(void);
extern int test_vectors_from_unit4(void);

/* Test function addresses to verify single instantiation */
void* get_fix16_add_address(void) {
    return (void*)fix16_add;
}

void* get_fix16_mul_address(void) {
    return (void*)fix16_mul;
}

void* get_fix16_sin_address(void) {
    return (void*)fix16_sin;
}

void* get_v2d_add_address(void) {
    return (void*)v2d_add;
}

int test_core_functions_unit1(void) {
    printf("=== Testing Core Functions (Unit 1) ===\n");
    
    /* Test basic arithmetic */
    fix16_t a = fix16_from_int(10);
    fix16_t b = fix16_from_int(5);
    fix16_t result = fix16_add(a, b);
    
    printf("fix16_add(10, 5) = %d (expected: %d)\n", 
           fix16_to_int(result), 15);
    assert(fix16_to_int(result) == 15);
    
    result = fix16_mul(a, b);
    printf("fix16_mul(10, 5) = %d (expected: %d)\n", 
           fix16_to_int(result), 50);
    assert(fix16_to_int(result) == 50);
    
    /* Test conversion functions */
    fix16_t pi_approx = fix16_from_float(3.14159f);
    float back_to_float = fix16_to_float(pi_approx);
    printf("float->fix16->float: 3.14159 -> %f\n", back_to_float);
    assert(back_to_float > 3.14f && back_to_float < 3.15f);
    
    printf("Unit 1 tests: PASSED\n\n");
    return 1;
}

int main(void) {
    printf("=========================================\n");
    printf("MULTI-UNIT TEST: Single Instantiation\n");
    printf("=========================================\n");
    
    /* Test this unit */
    int result1 = test_core_functions_unit1();
    
    /* Test other units */
    int result2 = test_basic_math_from_unit2();
    int result3 = test_trigonometry_from_unit3();
    int result4 = test_vectors_from_unit4();
    
    if (result1 && result2 && result3 && result4) {
        printf("========================================\n");
        printf("ALL TESTS PASSED!\n");
        printf("Single instantiation test successful.\n");
        printf("========================================\n");
        return 0;
    } else {
        printf("========================================\n");
        printf("SOME TESTS FAILED!\n");
        printf("========================================\n");
        return 1;
    }
} 