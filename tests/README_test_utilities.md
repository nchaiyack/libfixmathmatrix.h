# Test Utilities for libfixmathmatrix

This directory contains shared testing utilities that provide consistent, reusable testing patterns across all test families in the libfixmathmatrix test suite.

## Overview

The `test_utilities.h` header provides:
- **Fixed-point comparison functions** with configurable tolerance
- **Vector and matrix testing utilities** for multi-dimensional data
- **Enhanced assertion macros** with informative error messages  
- **String testing helpers** for text output validation
- **Test reporting functions** for consistent output formatting

## Quick Start

```c
#include <stdio.h>
#define LIBFIXMATHMATRIX_IMPLEMENTATION
#include <libfixmathmatrix_final.h>
#include "test_utilities.h"

int main() {
    // Test fixed-point arithmetic with tolerance
    assert(fix16_approx_equal(fix16_sqrt(F16(4.0)), F16(2.0), 100));
    
    // Test vectors
    v3d vec1 = {F16(1.0), F16(2.0), F16(3.0)};
    v3d vec2 = {F16(1.0), F16(2.0), F16(3.0)};
    assert(v3d_approx_equal(&vec1, &vec2, 10));
    
    // Enhanced assertions with messages
    ASSERT_FIX16_EQUAL(fix16_pi, F16(3.14159));
    
    return 0;
}
```

## Core Functions

### Fixed-Point Testing

| Function | Description | Example |
|----------|-------------|---------|
| `fix16_approx_equal(a, b, tolerance)` | Compare fix16_t values within tolerance | `fix16_approx_equal(result, F16(2.0), 100)` |
| `fix16_approx_equal_decimal(fix16, decimal, tol)` | Compare fix16_t to decimal value | `fix16_approx_equal_decimal(fix16_pi, 3.14159, 0.001)` |
| `fix16_approx_equal_percent(a, b, percent)` | Compare with percentage tolerance | `fix16_approx_equal_percent(calc, expected, 5.0)` |
| `fix16_exact_equal(a, b)` | Exact equality check | `fix16_exact_equal(fix16_one, 65536)` |

### Vector Testing

| Function | Description | Example |
|----------|-------------|---------|
| `v2d_approx_equal(v1, v2, tolerance)` | Compare 2D vectors | `v2d_approx_equal(&result, &expected, 100)` |
| `v3d_approx_equal(v1, v2, tolerance)` | Compare 3D vectors | `v3d_approx_equal(&result, &expected, 100)` |
| `v3d_approx_equal_decimal(v, x, y, z, tol)` | Compare 3D vector to decimal components | `v3d_approx_equal_decimal(&vec, 1.0, 2.0, 3.0, 0.01)` |

### Matrix Testing  

| Function | Description | Example |
|----------|-------------|---------|
| `mf16_approx_equal(m1, m2, tolerance)` | Compare matrices element-wise | `mf16_approx_equal(&result, &expected, 100)` |
| `mf16_same_dimensions(m1, m2)` | Check dimension compatibility | `mf16_same_dimensions(&A, &B)` |
| `mf16_no_errors(matrix)` | Check error flags | `mf16_no_errors(&result)` |
| `mf16_is_square(matrix)` | Check if matrix is square | `mf16_is_square(&matrix)` |

### String Testing

| Function | Description | Example |
|----------|-------------|---------|
| `str_approx_equal(actual, prefix)` | Check string starts with prefix | `str_approx_equal(buffer, "3.14")` |
| `str_exact_equal(actual, expected)` | Exact string equality | `str_exact_equal(buffer, "2.718")` |
| `str_contains(haystack, needle)` | Check substring existence | `str_contains(output, "PASS")` |

### Enhanced Assertion Macros

| Macro | Description | Example |
|-------|-------------|---------|
| `ASSERT_MSG(condition, message)` | Assert with custom message | `ASSERT_MSG(x > 0, "Value must be positive")` |
| `ASSERT_FIX16_EQUAL(actual, expected)` | Fix16 equality with default tolerance | `ASSERT_FIX16_EQUAL(result, F16(2.0))` |
| `ASSERT_FIX16_APPROX(actual, expected, tol)` | Fix16 equality with custom tolerance | `ASSERT_FIX16_APPROX(result, expected, 500)` |
| `ASSERT_STR_PREFIX(actual, prefix)` | String prefix assertion | `ASSERT_STR_PREFIX(buffer, "3.14")` |
| `ASSERT_V3D_APPROX(actual, expected, tol)` | 3D vector assertion | `ASSERT_V3D_APPROX(&result, &expected, 100)` |

### Test Reporting

| Function | Description | Example |
|----------|-------------|---------|
| `print_test_section(name)` | Print formatted section header | `print_test_section("Matrix Tests")` |
| `print_test_summary(name, passed, total)` | Print test results summary | `print_test_summary("Math Tests", 15, 16)` |

## Integration Examples

### Basic Usage Pattern
```c
void test_arithmetic() {
    print_test_section("Arithmetic Tests");
    
    fix16_t result = fix16_add(fix16_pi, fix16_e);
    assert(fix16_approx_equal(result, F16(5.86), 1000));
    
    printf("✅ Addition test passed\n");
}
```

### Matrix Testing Pattern
```c
void test_matrix_multiply() {
    mf16 A, B, result;
    // ... initialize matrices ...
    
    mf16_mul(&result, &A, &B);
    
    // Verify no computational errors
    assert(mf16_no_errors(&result));
    
    // Check dimensions
    assert(result.rows == A.rows && result.columns == B.columns);
    
    // Compare with expected result
    assert(mf16_approx_equal(&result, &expected, 100));
}
```

### Vector Testing Pattern
```c
void test_cross_product() {
    v3d i = {F16(1.0), F16(0.0), F16(0.0)};
    v3d j = {F16(0.0), F16(1.0), F16(0.0)};
    v3d result;
    
    v3d_cross(&result, &i, &j);
    
    // Should equal k = (0, 0, 1)
    assert(v3d_approx_equal_decimal(&result, 0.0, 0.0, 1.0, 0.01));
}
```

## Tolerance Guidelines

### Recommended Tolerances

| Operation Type | Suggested Tolerance | Rationale |
|----------------|-------------------|-----------|
| **Constants** | `10` | Fixed constants should be exact |
| **Basic arithmetic** | `100-1000` | Rounding errors in Q16.16 format |
| **Trigonometric** | `1000` | Approximation algorithms introduce error |
| **Exponential/Log** | `2000` | Complex mathematical functions |
| **Vector norms** | `1000` | Square root operations |
| **Matrix operations** | `100-1000` | Depends on matrix size and condition |

### Decimal Tolerances
- **Constants**: `0.001` (3 decimal places)
- **General calculations**: `0.01` (2 decimal places)  
- **Complex operations**: `0.1` (1 decimal place)

## Files

- **`test_utilities.h`** - Main header with all utility functions
- **`example_usage_test.c`** - Comprehensive usage examples
- **`autogenerated-diagnostics/basic_diagnostics.c`** - Real-world usage example

## Contributing

When adding new utilities:

1. **Follow naming conventions**: `type_operation_qualifier()` format
2. **Add comprehensive documentation** with `@brief`, `@param`, `@return`
3. **Use `static inline`** for header-only functions
4. **Include usage examples** in comments
5. **Add corresponding assertion macros** for convenience
6. **Test thoroughly** with the example usage file

## License

These utilities are part of the libfixmathmatrix test suite and follow the same licensing terms as the main library. 