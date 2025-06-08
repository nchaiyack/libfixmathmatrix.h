/* Test file 2: Unit that only uses declarations */
#include "libfixmathmatrix_final.h"
#include <stdio.h>
#include <assert.h>

/* Forward declarations for address comparison functions */
extern void* get_fix16_add_address(void);
extern void* get_fix16_mul_address(void);

/* Test function addresses to verify single instantiation */
void* get_fix16_add_address_unit2(void) {
    return (void*)fix16_add;
}

void* get_fix16_mul_address_unit2(void) {
    return (void*)fix16_mul;
}

int test_basic_math_from_unit2(void) {
    printf("=== Testing Basic Math (Unit 2) ===\n");
    
    /* Test that we're using the same function instances */
    void* addr1 = get_fix16_add_address();
    void* addr2 = get_fix16_add_address_unit2();
    printf("fix16_add address from unit 1: %p\n", addr1);
    printf("fix16_add address from unit 2: %p\n", addr2);
    assert(addr1 == addr2);
    printf("✓ Function addresses match - single instantiation confirmed\n");
    
    /* Test division and modulo */
    fix16_t a = fix16_from_int(100);
    fix16_t b = fix16_from_int(3);
    fix16_t quotient = fix16_div(a, b);
    fix16_t remainder = fix16_mod(a, b);
    
    printf("fix16_div(100, 3) = %f\n", fix16_to_float(quotient));
    printf("fix16_mod(100, 3) = %f\n", fix16_to_float(remainder));
     /* <add asserts here> */

    /* Test square root */
    fix16_t sqrt_val = fix16_sqrt(fix16_from_int(16));
    printf("fix16_sqrt(16) = %d (expected: 4)\n", fix16_to_int(sqrt_val));
    assert(fix16_to_int(sqrt_val) == 4);
    
    /* Test min/max/clamp */
    fix16_t min_val = fix16_min(fix16_from_int(10), fix16_from_int(5));
    fix16_t max_val = fix16_max(fix16_from_int(10), fix16_from_int(5));
    fix16_t clamped = fix16_clamp(fix16_from_int(15), fix16_from_int(5), fix16_from_int(10));
    
    assert(fix16_to_int(min_val) == 5);
    assert(fix16_to_int(max_val) == 10);
    assert(fix16_to_int(clamped) == 10);
    
    printf("Unit 2 tests: PASSED\n\n");
    return 1;
} 