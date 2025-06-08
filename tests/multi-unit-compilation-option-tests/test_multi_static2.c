/* Test file 2: Static instantiation - Unit 2 */
#define LIBFIXMATHMATRIX_STATIC
#define LIBFIXMATHMATRIX_IMPLEMENTATION
#include "libfixmathmatrix_final.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>

/* Test function addresses - these should be different from other units */
void* get_static_fix16_add_address_unit2(void) {
    return (void*)fix16_add;
}

void* get_static_fix16_mul_address_unit2(void) {
    return (void*)fix16_mul;
}

int test_static_basic_math_from_unit2(void) {
    printf("=== Testing Static Basic Math (Unit 2) ===\n");
    
    /* These functions are static to this compilation unit */
    fix16_t a = fix16_from_int(15);
    fix16_t b = fix16_from_int(4);
    
    /* Test that static functions work correctly */
    fix16_t quotient = fix16_div(a, b);
    fix16_t remainder = fix16_mod(a, b);
    
    printf("fix16_div(15, 4) = %f\n", fix16_to_float(quotient));
    /* <add asserts here> */
    printf("fix16_mod(15, 4) = %f\n", fix16_to_float(remainder));
    /* <add asserts here> */
    
    /* Test square root with local static function */
    fix16_t sqrt_val = fix16_sqrt(fix16_from_int(25));
    printf("fix16_sqrt(25) = %d (expected: 5)\n", fix16_to_int(sqrt_val));
    assert(fix16_to_int(sqrt_val) == 5);
    
    /* Test min/max/clamp with static functions */
    fix16_t min_val = fix16_min(fix16_from_int(7), fix16_from_int(12));
    fix16_t max_val = fix16_max(fix16_from_int(7), fix16_from_int(12));
    fix16_t clamped = fix16_clamp(fix16_from_int(20), fix16_from_int(5), fix16_from_int(15));
    
    printf("min(7, 12) = %d, max(7, 12) = %d, clamp(20, 5, 15) = %d\n",
           fix16_to_int(min_val), fix16_to_int(max_val), fix16_to_int(clamped));
    
    assert(fix16_to_int(min_val) == 7);
    assert(fix16_to_int(max_val) == 12);
    assert(fix16_to_int(clamped) == 15);
    
    /* Test lerp functions */
    fix16_t lerp_result = fix16_lerp16(fix16_from_int(0), fix16_from_int(10), 32768); /* 50% */
    printf("lerp16(0, 10, 50%%) = %d (expected: 5)\n", fix16_to_int(lerp_result));
        /* <add asserts here> */

    /* Print function address for comparison with other units */
    printf("Unit 2 fix16_add address: %p (static)\n", get_static_fix16_add_address_unit2());
        /* <add asserts here> */

    /* Test string conversion */
    char buffer[32];
    fix16_to_str(fix16_from_float(3.14159f), buffer, 3);
    int len = strlen(buffer);
    printf("fix16_to_str(3.14159, 3) = \"%s\" (length: %d)\n", buffer, len);
    /* <add asserts here> */

    fix16_t from_str = fix16_from_str("2.718");
    printf("fix16_from_str(\"2.718\") = %f\n", fix16_to_float(from_str));
    /* <add asserts here> */

    printf("Static Unit 2 tests: PASSED\n\n");
    return 1;
} 