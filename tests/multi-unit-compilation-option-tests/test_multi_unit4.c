/* Test file 4: Unit that tests vector and matrix operations */
#include "libfixmathmatrix_final.h"
#include <stdio.h>
#include <assert.h>

/* Forward declarations for address comparison functions */
extern void* get_v2d_add_address(void);

/* Test function addresses to verify single instantiation */
void* get_v2d_add_address_unit4(void) {
    return (void*)v2d_add;
}

int test_vectors_from_unit4(void) {
    printf("=== Testing Vectors & Matrices (Unit 4) ===\n");
    
    /* Test that we're using the same function instances */
    void* addr1 = get_v2d_add_address();
    void* addr2 = get_v2d_add_address_unit4();
    printf("v2d_add address from unit 1: %p\n", addr1);
    printf("v2d_add address from unit 4: %p\n", addr2);
    assert(addr1 == addr2);
    printf("✓ Function addresses match - single instantiation confirmed\n");
    
    /* Test 2D vector operations */
    v2d vec_a = {fix16_from_int(3), fix16_from_int(4)};
    v2d vec_b = {fix16_from_int(1), fix16_from_int(2)};
    v2d vec_result;
    
    v2d_add(&vec_result, &vec_a, &vec_b);
    printf("v2d_add([3,4], [1,2]) = [%d,%d] (expected: [4,6])\n", 
           fix16_to_int(vec_result.x), fix16_to_int(vec_result.y));
    assert(fix16_to_int(vec_result.x) == 4);
    assert(fix16_to_int(vec_result.y) == 6);
    
    v2d_sub(&vec_result, &vec_a, &vec_b);
    printf("v2d_sub([3,4], [1,2]) = [%d,%d] (expected: [2,2])\n", 
           fix16_to_int(vec_result.x), fix16_to_int(vec_result.y));
    assert(fix16_to_int(vec_result.x) == 2);
    assert(fix16_to_int(vec_result.y) == 2);
    
    fix16_t dot_product = v2d_dot(&vec_a, &vec_b);
    printf("v2d_dot([3,4], [1,2]) = %d (expected: 11)\n", fix16_to_int(dot_product));
    assert(fix16_to_int(dot_product) == 11);
    
    fix16_t norm = v2d_norm(&vec_a);
    printf("v2d_norm([3,4]) = %f (expected: 5.0)\n", fix16_to_float(norm));
       /* <add asserts here> */

    /* Test 3D vector operations */
    v3d vec3_a = {fix16_from_int(1), fix16_from_int(2), fix16_from_int(3)};
    v3d vec3_b = {fix16_from_int(4), fix16_from_int(5), fix16_from_int(6)};
    v3d vec3_result;
    
    
    v3d_cross(&vec3_result, &vec3_a, &vec3_b);
    printf("v3d_cross([1,2,3], [4,5,6]) = [%d,%d,%d] (expected: [-3,6,-3])\n", 
           fix16_to_int(vec3_result.x), fix16_to_int(vec3_result.y), fix16_to_int(vec3_result.z));
 /* <add asserts here> */

    /* Test matrix operations */
    mf16 matrix_a = {{2, 2, 0, {{fix16_from_int(1), fix16_from_int(2)}, 
                                 {fix16_from_int(3), fix16_from_int(4)}}}};
    mf16 matrix_b = {{2, 2, 0, {{fix16_from_int(2), fix16_from_int(0)}, 
                                 {fix16_from_int(1), fix16_from_int(2)}}}};
    mf16 matrix_result;
    
    mf16_add(&matrix_result, &matrix_a, &matrix_b);
    printf("Matrix addition result:\n");
    printf("  [%d %d]\n", fix16_to_int(matrix_result.data[0][0]), fix16_to_int(matrix_result.data[0][1]));
    printf("  [%d %d]\n", fix16_to_int(matrix_result.data[1][0]), fix16_to_int(matrix_result.data[1][1]));
    printf("Expected: [3 2; 4 6]\n");
        /* <add asserts here> */
    
    mf16_transpose(&matrix_result, &matrix_a);
    printf("Matrix transpose result:\n");
    printf("  [%d %d]\n", fix16_to_int(matrix_result.data[0][0]), fix16_to_int(matrix_result.data[0][1]));
    printf("  [%d %d]\n", fix16_to_int(matrix_result.data[1][0]), fix16_to_int(matrix_result.data[1][1]));
    printf("Expected: [1 3; 2 4]\n");
        /* <add asserts here> */
    
    /* Test quaternion operations */
    qf16 quat_a = {fix16_one, 0, 0, 0};  /* Identity quaternion */
    qf16 quat_b = {fix16_one >> 1, fix16_one >> 1, 0, 0};  /* Rotation quaternion */
    qf16 quat_result;
    
    qf16_mul(&quat_result, &quat_a, &quat_b);
    printf("Quaternion multiplication result: [%f, %f, %f, %f]\n",
           fix16_to_float(quat_result.a), fix16_to_float(quat_result.b),
           fix16_to_float(quat_result.c), fix16_to_float(quat_result.d));
        /* <add asserts here> */

    fix16_t quat_norm = qf16_norm(&quat_b);
    printf("Quaternion norm: %f\n", fix16_to_float(quat_norm));
        /* <add asserts here> */

    printf("Unit 4 tests: PASSED\n\n");
    return 1;
} 