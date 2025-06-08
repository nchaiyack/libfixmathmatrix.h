/*
 * test_utilities.h - Shared test utilities for libfixmathmatrix test suite
 * 
 * This header provides common testing utilities and helper functions
 * that can be used across different test families in the libfixmathmatrix
 * test suite.
 * 
 * USAGE:
 *   #include "test_utilities.h"
 * 
 * DEPENDENCIES:
 *   - assert.h
 *   - string.h 
 *   - math.h
 *   - libfixmathmatrix_final.h (or equivalent)
 */

#ifndef TEST_UTILITIES_H
#define TEST_UTILITIES_H

#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/* FIXED-POINT TESTING UTILITIES                                            */
/*===========================================================================*/

/**
 * @brief Test if two fix16_t values are approximately equal within tolerance
 * 
 * @param a First fixed-point value
 * @param b Second fixed-point value  
 * @param tolerance Maximum allowed difference (in fix16_t units)
 * @return 1 if values are within tolerance, 0 otherwise
 */
static inline int fix16_approx_equal(fix16_t a, fix16_t b, fix16_t tolerance) {
    fix16_t diff = fix16_abs(fix16_sub(a, b));
    return diff <= tolerance;
}

/**
 * @brief Test if a fix16_t value is approximately equal to a decimal value
 * 
 * @param actual The fix16_t value to test
 * @param expected The expected decimal value
 * @param tolerance Maximum allowed difference (in decimal units)
 * @return 1 if values are within tolerance, 0 otherwise
 */
static inline int fix16_approx_equal_decimal(fix16_t actual, double expected, double tolerance) {
    double actual_decimal = fix16_to_dbl(actual);
    return fabs(actual_decimal - expected) <= tolerance;
}

/**
 * @brief Test if a fix16_t value is exactly equal to expected
 * 
 * @param actual The fix16_t value to test
 * @param expected The expected fix16_t value
 * @return 1 if exactly equal, 0 otherwise
 */
static inline int fix16_exact_equal(fix16_t actual, fix16_t expected) {
    return actual == expected;
}

/**
 * @brief Test if a fix16_t value is within a percentage of expected value
 * 
 * @param actual The fix16_t value to test
 * @param expected The expected fix16_t value
 * @param percent_tolerance Tolerance as percentage (e.g., 5.0 for 5%)
 * @return 1 if within percentage tolerance, 0 otherwise
 */
static inline int fix16_approx_equal_percent(fix16_t actual, fix16_t expected, double percent_tolerance) {
    if (expected == 0) {
        return fix16_abs(actual) <= F16(percent_tolerance / 100.0);
    }
    fix16_t tolerance = fix16_abs(fix16_mul(expected, F16(percent_tolerance / 100.0)));
    return fix16_approx_equal(actual, expected, tolerance);
}

/*===========================================================================*/
/* STRING TESTING UTILITIES                                                 */
/*===========================================================================*/

/**
 * @brief Test if a string starts with expected prefix
 * 
 * @param actual The actual string
 * @param expected_prefix The expected prefix
 * @return 1 if string starts with prefix, 0 otherwise
 */
static inline int str_approx_equal(const char *actual, const char *expected_prefix) {
    return strncmp(actual, expected_prefix, strlen(expected_prefix)) == 0;
}

/**
 * @brief Test if two strings are exactly equal
 * 
 * @param actual The actual string
 * @param expected The expected string
 * @return 1 if strings are identical, 0 otherwise
 */
static inline int str_exact_equal(const char *actual, const char *expected) {
    return strcmp(actual, expected) == 0;
}

/**
 * @brief Test if a string contains a substring
 * 
 * @param haystack The string to search in
 * @param needle The substring to find
 * @return 1 if substring found, 0 otherwise
 */
static inline int str_contains(const char *haystack, const char *needle) {
    return strstr(haystack, needle) != NULL;
}

/*===========================================================================*/
/* VECTOR TESTING UTILITIES                                                 */
/*===========================================================================*/

/**
 * @brief Test if two 2D vectors are approximately equal
 * 
 * @param actual The actual v2d vector
 * @param expected The expected v2d vector
 * @param tolerance Maximum allowed difference per component
 * @return 1 if vectors are within tolerance, 0 otherwise
 */
static inline int v2d_approx_equal(const v2d *actual, const v2d *expected, fix16_t tolerance) {
    return fix16_approx_equal(actual->x, expected->x, tolerance) &&
           fix16_approx_equal(actual->y, expected->y, tolerance);
}

/**
 * @brief Test if two 3D vectors are approximately equal
 * 
 * @param actual The actual v3d vector
 * @param expected The expected v3d vector
 * @param tolerance Maximum allowed difference per component
 * @return 1 if vectors are within tolerance, 0 otherwise
 */
static inline int v3d_approx_equal(const v3d *actual, const v3d *expected, fix16_t tolerance) {
    return fix16_approx_equal(actual->x, expected->x, tolerance) &&
           fix16_approx_equal(actual->y, expected->y, tolerance) &&
           fix16_approx_equal(actual->z, expected->z, tolerance);
}

/**
 * @brief Test if a 3D vector has expected components
 * 
 * @param actual The actual v3d vector
 * @param x Expected x component (decimal)
 * @param y Expected y component (decimal)
 * @param z Expected z component (decimal)
 * @param tolerance Maximum allowed difference per component (decimal)
 * @return 1 if all components match within tolerance, 0 otherwise
 */
static inline int v3d_approx_equal_decimal(const v3d *actual, double x, double y, double z, double tolerance) {
    return fix16_approx_equal_decimal(actual->x, x, tolerance) &&
           fix16_approx_equal_decimal(actual->y, y, tolerance) &&
           fix16_approx_equal_decimal(actual->z, z, tolerance);
}

/*===========================================================================*/
/* MATRIX TESTING UTILITIES                                                 */
/*===========================================================================*/

/**
 * @brief Test if two matrices have the same dimensions
 * 
 * @param a First matrix
 * @param b Second matrix
 * @return 1 if dimensions match, 0 otherwise
 */
static inline int mf16_same_dimensions(const mf16 *a, const mf16 *b) {
    return (a->rows == b->rows) && (a->columns == b->columns);
}

/**
 * @brief Test if two matrices are approximately equal
 * 
 * @param actual The actual matrix
 * @param expected The expected matrix
 * @param tolerance Maximum allowed difference per element
 * @return 1 if matrices are within tolerance, 0 otherwise
 */
static inline int mf16_approx_equal(const mf16 *actual, const mf16 *expected, fix16_t tolerance) {
    if (!mf16_same_dimensions(actual, expected)) {
        return 0;
    }
    
    for (int row = 0; row < actual->rows; row++) {
        for (int col = 0; col < actual->columns; col++) {
            if (!fix16_approx_equal(actual->data[row][col], expected->data[row][col], tolerance)) {
                return 0;
            }
        }
    }
    return 1;
}

/**
 * @brief Test if a matrix has no errors
 * 
 * @param matrix The matrix to test
 * @return 1 if no errors, 0 if errors detected
 */
static inline int mf16_no_errors(const mf16 *matrix) {
    return matrix->errors == 0;
}

/**
 * @brief Test if a matrix is square
 * 
 * @param matrix The matrix to test
 * @return 1 if square, 0 otherwise
 */
static inline int mf16_is_square(const mf16 *matrix) {
    return matrix->rows == matrix->columns;
}

/*===========================================================================*/
/* QUATERNION TESTING UTILITIES                                             */
/*===========================================================================*/

/**
 * @brief Test if two quaternions are approximately equal
 * 
 * @param actual The actual quaternion
 * @param expected The expected quaternion
 * @param tolerance Maximum allowed difference per component
 * @return 1 if quaternions are within tolerance, 0 otherwise
 */
static inline int qf16_approx_equal(const qf16 *actual, const qf16 *expected, fix16_t tolerance) {
    return fix16_approx_equal(actual->a, expected->a, tolerance) &&
           fix16_approx_equal(actual->b, expected->b, tolerance) &&
           fix16_approx_equal(actual->c, expected->c, tolerance) &&
           fix16_approx_equal(actual->d, expected->d, tolerance);
}

/**
 * @brief Test if a quaternion is approximately normalized (|q| ≈ 1)
 * 
 * @param q The quaternion to test
 * @param tolerance Maximum allowed deviation from unit length
 * @return 1 if approximately normalized, 0 otherwise
 */
static inline int qf16_is_normalized(const qf16 *q, fix16_t tolerance) {
    fix16_t norm = qf16_norm(q);
    return fix16_approx_equal(norm, fix16_one, tolerance);
}

/*===========================================================================*/
/* ENHANCED ASSERTION MACROS                                                */
/*===========================================================================*/

/**
 * @brief Enhanced assertion with custom message
 */
#define ASSERT_MSG(condition, message) \
    do { \
        if (!(condition)) { \
            printf("ASSERTION FAILED: %s\n", message); \
            printf("  File: %s, Line: %d\n", __FILE__, __LINE__); \
            assert(condition); \
        } \
    } while(0)

/**
 * @brief Assert fix16_t values are approximately equal with automatic tolerance
 */
#define ASSERT_FIX16_EQUAL(actual, expected) \
    ASSERT_MSG(fix16_approx_equal((actual), (expected), 100), \
               "Fix16 values not approximately equal")

/**
 * @brief Assert fix16_t values are approximately equal with custom tolerance
 */
#define ASSERT_FIX16_APPROX(actual, expected, tolerance) \
    ASSERT_MSG(fix16_approx_equal((actual), (expected), (tolerance)), \
               "Fix16 values not within specified tolerance")

/**
 * @brief Assert string starts with expected prefix
 */
#define ASSERT_STR_PREFIX(actual, prefix) \
    ASSERT_MSG(str_approx_equal((actual), (prefix)), \
               "String does not start with expected prefix")

/**
 * @brief Assert vectors are approximately equal
 */
#define ASSERT_V3D_APPROX(actual, expected, tolerance) \
    ASSERT_MSG(v3d_approx_equal((actual), (expected), (tolerance)), \
               "3D vectors not approximately equal")

/*===========================================================================*/
/* TEST REPORTING UTILITIES                                                 */
/*===========================================================================*/

/**
 * @brief Print a test section header
 * 
 * @param section_name Name of the test section
 */
static inline void print_test_section(const char *section_name) {
    printf("\n=== %s ===\n", section_name);
}

/**
 * @brief Print test results summary
 * 
 * @param test_name Name of the test
 * @param passed Number of tests passed
 * @param total Total number of tests
 */
static inline void print_test_summary(const char *test_name, int passed, int total) {
    printf("\n%s: %d/%d tests passed", test_name, passed, total);
    if (passed == total) {
        printf(" ✅\n");
    } else {
        printf(" ❌\n");
    }
}

#ifdef __cplusplus
}
#endif

#endif /* TEST_UTILITIES_H */ 