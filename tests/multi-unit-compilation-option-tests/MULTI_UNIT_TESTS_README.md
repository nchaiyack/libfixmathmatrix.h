# Multi-Unit Compilation Tests for libfixmathmatrix_final.h

This directory contains comprehensive tests to validate that the refactored `libfixmathmatrix_final.h` header-only library works correctly with multiple compilation units in two different usage patterns.

## Test Scenarios

### 1. Single Instantiation Pattern

**Files:** `test_multi_unit_single.c`, `test_multi_unit2.c`, `test_multi_unit3.c`, `test_multi_unit4.c`

**Pattern:**
- Multiple `.c` files include `libfixmathmatrix_final.h`
- **Only ONE file** defines `LIBFIXMATHMATRIX_IMPLEMENTATION`
- All units share the same function instances
- Function addresses should be identical across units

**Verification:**
- Tests that all compilation units can use the library functions
- Validates that function addresses are identical across units (single instantiation)
- Ensures no multiple definition errors during linking
- Smaller binary size due to single function instances

### 2. Static Instantiation Pattern

**Files:** `test_multi_static1.c`, `test_multi_static2.c`, `test_multi_static3.c`

**Pattern:**
- Multiple `.c` files include `libfixmathmatrix_final.h`
- **All files** define both `LIBFIXMATHMATRIX_STATIC` and `LIBFIXMATHMATRIX_IMPLEMENTATION`
- Each unit gets its own static copy of functions
- Function addresses should be different across units

**Verification:**
- Tests that static instantiation works without linking conflicts
- Validates that each unit has its own function copies
- Ensures no symbol conflicts during linking
- Larger binary size due to code duplication

## Test Coverage

Both test scenarios comprehensively test all functional groups:

### Group A: Core Types & Basic Arithmetic
- `fix16_from_int`, `fix16_to_float`, `fix16_add`, `fix16_mul`
- Type conversions and basic operations

### Group B: Advanced Math Functions  
- `fix16_div`, `fix16_sqrt`, `fix16_sin`, `fix16_exp`
- Mathematical functions and trigonometry

### Group C: String & Utility Functions
- `fix16_to_str`, `fix16_from_str`
- String conversion utilities

### Group D: Vector Operations
- `v2d_add`, `v2d_dot`, `v3d_cross`
- 2D and 3D vector mathematics

### Group E: Matrix Operations
- `mf16_add`, `mf16_transpose`
- Matrix manipulation functions

### Group F: Quaternion Operations
- `qf16_mul`, `qf16_normalize`
- Quaternion mathematics

### Group G: Pretty Printing Functions
- `print_mf16` (if stdio available)
- Formatted output functions

## Running the Tests

### Quick Start

```bash
# Run all tests with comprehensive validation
./run_multi_unit_tests.sh
```

### Manual Testing

```bash
# Build both test executables
make all

# Run individual tests
make test-single
make test-static

# Run both tests
make test
```

### Advanced Analysis

```bash
# Analyze symbol tables
make analyze-symbols

# Compare binary sizes
make size-comparison

# Verify linking behavior
make verify-single
make verify-static

# Show test information
make test-info
```

## Expected Results

### Single Instantiation Test
- **✓ PASS**: Function addresses identical across units
- **✓ PASS**: All library functions work correctly
- **✓ PASS**: No multiple definition errors
- **✓ PASS**: Smaller binary size

### Static Instantiation Test  
- **✓ PASS**: Function addresses different across units
- **✓ PASS**: All library functions work correctly
- **✓ PASS**: No symbol conflicts
- **✓ PASS**: Larger binary size (expected)

## Key Validations

### 1. Function Address Consistency
The tests verify that:
- Single instantiation: Same function addresses across units
- Static instantiation: Different function addresses across units

### 2. Linking Behavior
- Single instantiation: One definition rule satisfied
- Static instantiation: Static linkage prevents conflicts

### 3. Functional Correctness
- All arithmetic operations produce correct results
- All mathematical functions work properly
- All vector/matrix operations behave correctly

### 4. Configuration Options
- `LIBFIXMATHMATRIX_STATIC` correctly affects function attributes
- `LIBFIXMATHMATRIX_IMPLEMENTATION` properly includes implementations
- `FIXMATH_FUNC_ATTRS` adapts to static vs extern linkage

## File Structure

```
.
├── Makefile                        # Build system
├── run_multi_unit_tests.sh        # Comprehensive test runner
├── MULTI_UNIT_TESTS_README.md     # This documentation
│
├── test_multi_unit_single.c       # Single instantiation: main unit
├── test_multi_unit2.c             # Single instantiation: unit 2
├── test_multi_unit3.c             # Single instantiation: unit 3
├── test_multi_unit4.c             # Single instantiation: unit 4
│
├── test_multi_static1.c           # Static instantiation: unit 1
├── test_multi_static2.c           # Static instantiation: unit 2
└── test_multi_static3.c           # Static instantiation: unit 3
```

## Build Requirements

- GCC compiler (or compatible)
- Make utility
- `libfixmathmatrix_final.h`
- Implementation group files (`implementation_group_*.h`)
- Math library (`-lm`)

## Troubleshooting

### Build Failures
- Ensure `libfixmathmatrix_final.h` is present
- Check that implementation group files exist
- Verify GCC is available and working

### Test Failures
- Review function address comparisons in output
- Check for unexpected linking errors
- Validate that correct macros are defined

### Symbol Analysis Issues
- Use `nm` and `objdump` to inspect symbols manually
- Check for unexpected global symbols in static build
- Verify function attributes are applied correctly

## What This Proves

These tests demonstrate that `libfixmathmatrix_final.h` successfully implements a proper header-only library pattern:

1. **Single Instantiation**: Works like a traditional library with shared functions
2. **Static Instantiation**: Works like a pure header-only library with private functions
3. **No Conflicts**: Both patterns compile and link without errors
4. **Full Functionality**: All library features work in both patterns
5. **Proper Attributes**: Function attributes adapt correctly to usage pattern

This validates that the refactoring from a monolithic implementation to a proper header-only library was successful and maintains all intended functionality while providing flexible usage patterns. 