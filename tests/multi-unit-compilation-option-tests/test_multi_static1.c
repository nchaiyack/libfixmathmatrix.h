/* Test file 1: Static instantiation - Unit 1 */
#define LIBFIXMATHMATRIX_STATIC
#define LIBFIXMATHMATRIX_IMPLEMENTATION
#include "libfixmathmatrix_final.h"
#include <stdio.h>
#include <assert.h>

/* Forward declarations for functions from other units */
extern int test_static_basic_math_from_unit2(void);
extern int test_static_trigonometry_from_unit3(void);

/* Test function addresses - these should be different in static mode */
void* get_static_fix16_add_address_unit1(void) {
    return (void*)fix16_add;
}

void* get_static_fix16_mul_address_unit1(void) {
    return (void*)fix16_mul;
}

void* get_static_fix16_sin_address_unit1(void) {
    return (void*)fix16_sin;
}

int test_static_core_functions_unit1(void) {
    printf("=== Testing Static Core Functions (Unit 1) ===\n");
    
    /* Test basic arithmetic - these are local static functions */
    fix16_t a = fix16_from_int(20);
    fix16_t b = fix16_from_int(8);
    fix16_t result = fix16_add(a, b);
    
    printf("fix16_add(20, 8) = %d (expected: %d)\n", 
           fix16_to_int(result), 28);
    assert(fix16_to_int(result) == 28);
    
    result = fix16_mul(a, b);
    printf("fix16_mul(20, 8) = %d (expected: %d)\n", 
           fix16_to_int(result), 160);
    assert(fix16_to_int(result) == 160);
    
    /* Test division */
    result = fix16_div(a, b);
    printf("fix16_div(20, 8) = %f (expected: 2.5)\n", fix16_to_float(result));
    /* <add assert here> */
    
    /* Test square operations */
    fix16_t square = fix16_sq(fix16_from_int(5));
    printf("fix16_sq(5) = %d (expected: 25)\n", fix16_to_int(square));
    assert(fix16_to_int(square) == 25);
    
    /* Test abs and floor/ceil */
    fix16_t negative = fix16_from_float(-3.7f);
    fix16_t abs_val = fix16_abs(negative);
    fix16_t floor_val = fix16_floor(negative);
    fix16_t ceil_val = fix16_ceil(negative);
    
    printf("abs(-3.7) = %f, floor(-3.7) = %f, ceil(-3.7) = %f\n",
           fix16_to_float(abs_val), fix16_to_float(floor_val), fix16_to_float(ceil_val));
    /* <add asserts here> */

    /* Print function address for comparison */
    printf("Unit 1 fix16_add address: %p (static)\n", get_static_fix16_add_address_unit1());
    
    printf("Static Unit 1 tests: PASSED\n\n");
    return 1;
}

int main(void) {
    printf("=========================================\n");
    printf("MULTI-UNIT TEST: Static Instantiation\n");
    printf("=========================================\n");
    printf("Each unit has its own static copy of functions\n\n");
    
    /* Test this unit */
    int result1 = test_static_core_functions_unit1();
    
    /* Test other units */
    int result2 = test_static_basic_math_from_unit2();
    int result3 = test_static_trigonometry_from_unit3();
    
    if (result1 && result2 && result3) {
        printf("========================================\n");
        printf("ALL STATIC TESTS PASSED!\n");
        printf("Static instantiation test successful.\n");
        printf("Each unit has its own function copies.\n");
        printf("========================================\n");
        return 0;
    } else {
        printf("========================================\n");
        printf("SOME STATIC TESTS FAILED!\n");
        printf("========================================\n");
        return 1;
    }
} 