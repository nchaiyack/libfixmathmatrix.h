/* Implementation - Group D: Vector Operations */

/* Forward declarations for functions that may call each other */
FIXMATH_FUNC_ATTRS fix16_t fa16_dot(const fix16_t *a, uint_fast8_t a_stride,
                                   const fix16_t *b, uint_fast8_t b_stride,
                                   uint_fast8_t n);
FIXMATH_FUNC_ATTRS fix16_t fa16_norm(const fix16_t *a, uint_fast8_t a_stride, uint_fast8_t n);

/* Count leading zeros - processor optimized when available (duplicate detection) */
#ifndef CLZ_ALREADY_DEFINED
#define CLZ_ALREADY_DEFINED
#ifdef __GNUC__
#define clz(x) (__builtin_clzl(x) - (8 * sizeof(long) - 32))
#else
static uint8_t clz(uint32_t x)
{
    uint8_t result = 0;
    if (x == 0) return 32;
    while (!(x & 0xF0000000)) { result += 4; x <<= 4; }
    while (!(x & 0x80000000)) { result += 1; x <<= 1; }
    return result;
}
#endif
#endif

/* Scale a value by a power of 2 with overflow detection */
static fix16_t scale_value(fix16_t value, int_fast8_t scale)
{
    if (scale > 0)
    {
        fix16_t temp = value << scale;
        if (temp >> scale != value)
            return fix16_overflow;
        else
            return temp;
    }
    else if (scale < 0)
    {
        return value >> -scale;
    }
    else
    {
        return value;
    }
}

/* Integer log2 for array size calculations (32-bit version) */
#ifdef FIXMATH_NO_64BIT
static uint_fast8_t ilog2(uint_fast8_t v)
{
    uint_fast8_t result = 0;
    if (v & 0xF0) { result += 4; v >>= 4; }
    while (v) { result++; v >>= 1; }
    return result;
}
#endif

/*---------------------------------------------------------------------------*/
/* ARRAY UTILITY FUNCTIONS                                                   */
/*---------------------------------------------------------------------------*/

/* Memory aliasing helper - used when dest overlaps with source arrays */
FIXMATH_FUNC_ATTRS void fa16_unalias(void *dest, void **a, void **b, void *tmp, unsigned size)
{
    if (dest == *a)
    {
        memcpy(tmp, *a, size);
        *a = tmp;
        
        if (dest == *b)
            *b = tmp;
    }
    else if (dest == *b)
    {
        memcpy(tmp, *b, size);
        *b = tmp;
    }
}

/* Vector dot product - conditional compilation for 32-bit vs 64-bit arithmetic */
#ifdef FIXMATH_NO_64BIT

/* Vector dot product - 32-bit version without 64-bit arithmetic */
FIXMATH_FUNC_ATTRS fix16_t fa16_dot(const fix16_t *a, uint_fast8_t a_stride,
                                   const fix16_t *b, uint_fast8_t b_stride,
                                   uint_fast8_t n)
{
    fix16_t sum = 0;
    
    while (n--)
    {
        // Compute result
        if (*a != 0 && *b != 0)
        {
            fix16_t product = fix16_mul(*a, *b);
            sum = fix16_add(sum, product);
            
            if (sum == fix16_overflow || product == fix16_overflow)
                return fix16_overflow;
        }
        
        // Go to next item
        a += a_stride;
        b += b_stride;
    }
    
    return sum;
}

/* Vector norm calculation - 32-bit version */
FIXMATH_FUNC_ATTRS fix16_t fa16_norm(const fix16_t *a, uint_fast8_t a_stride, uint_fast8_t n)
{
    fix16_t sum = 0;
    fix16_t max = 0;
    
    // Calculate inclusive OR of all values to find out the maximum.
    {
        uint_fast8_t i;
        const fix16_t *p = a;
        for (i = 0; i < n; i++, p += a_stride)
        {
            max |= fix16_abs(*p);
        }
    }
    
    // To avoid overflows, the values before squaring can be max 128.0,
    // i.e. v & 0xFF800000 must be 0. Also, to avoid overflow in sum,
    // we need additional log2(n) bits of space.
    int_fast8_t scale = clz(max) - 9 - ilog2(n) / 2;
    
    while (n--)
    {
        fix16_t val = scale_value(*a, scale);
        fix16_t product = fix16_mul(val, val);
        sum = fix16_add(sum, product);
        
        a += a_stride;
    }
    
    if (sum == fix16_overflow)
        return sum;
    
    fix16_t result = fix16_sqrt(sum);
    return scale_value(result, -scale);
}

#else

/* Vector dot product - optimized 64-bit version for ARM SMLAL instruction */
FIXMATH_FUNC_ATTRS fix16_t fa16_dot(const fix16_t *a, uint_fast8_t a_stride,
                                   const fix16_t *b, uint_fast8_t b_stride,
                                   uint_fast8_t n)
{
    int64_t sum = 0;
    
    while (n--)
    {
        if (*a != 0 && *b != 0)
        {
            sum += (int64_t)(*a) * (*b);
        }
        
        // Go to next item
        a += a_stride;
        b += b_stride;
    }
    
    // The upper 17 bits should all be the same (the sign).
    uint32_t upper = sum >> 47;
    if (sum < 0)
    {
        upper = ~upper;
        
        #ifndef FIXMATH_NO_ROUNDING
        // This adjustment is required in order to round -1/2 correctly
        sum--;
        #endif
    }
    
    #ifndef FIXMATH_NO_OVERFLOW
    if (upper)
        return fix16_overflow;
    #endif
    
    fix16_t result = sum >> 16;

    #ifndef FIXMATH_NO_ROUNDING
    result += (sum & 0x8000) >> 15;
    #endif
    
    return result;
}

/* Vector norm calculation - optimized 64-bit version */
FIXMATH_FUNC_ATTRS fix16_t fa16_norm(const fix16_t *a, uint_fast8_t a_stride, uint_fast8_t n)
{
    int64_t sum = 0;
    
    while (n--)
    {
        if (*a != 0)
        {
            sum += (int64_t)(*a) * (*a);
        }
        
        a += a_stride;
    }
    
    int_fast8_t scale = 0;
    uint32_t highpart = (uint32_t)(sum >> 32);
    uint32_t lowpart = (uint32_t)sum;
    if (highpart)
        scale = 33 - clz(highpart);
    else if (lowpart & 0x80000000)
        scale = 1;
    
    if (scale & 1) scale++;
    
    fix16_t result = fix16_sqrt((uint32_t)(sum >> scale));
    result = scale_value(result, scale / 2 - 8);
    
    return result;
}

#endif

/*---------------------------------------------------------------------------*/
/* 2D VECTOR OPERATIONS                                                      */
/*---------------------------------------------------------------------------*/

FIXMATH_FUNC_ATTRS void v2d_add(v2d *dest, const v2d *a, const v2d *b)
{
    dest->x = fix16_add(a->x, b->x);
    dest->y = fix16_add(a->y, b->y);
}

FIXMATH_FUNC_ATTRS void v2d_sub(v2d *dest, const v2d *a, const v2d *b)
{
    dest->x = fix16_sub(a->x, b->x);
    dest->y = fix16_sub(a->y, b->y);
}

FIXMATH_FUNC_ATTRS void v2d_mul_s(v2d *dest, const v2d *a, fix16_t b)
{
    dest->x = fix16_mul(a->x, b);
    dest->y = fix16_mul(a->y, b);
}

FIXMATH_FUNC_ATTRS void v2d_div_s(v2d *dest, const v2d *a, fix16_t b)
{
    dest->x = fix16_div(a->x, b);
    dest->y = fix16_div(a->y, b);
}

FIXMATH_FUNC_ATTRS fix16_t v2d_dot(const v2d *a, const v2d *b)
{
    return fix16_add(fix16_mul(a->x, b->x), fix16_mul(a->y, b->y));
}

FIXMATH_FUNC_ATTRS fix16_t v2d_norm(const v2d *a)
{
    return fa16_norm(&a->x, &a->y - &a->x, 2);
}

FIXMATH_FUNC_ATTRS void v2d_normalize(v2d *dest, const v2d *a)
{
    v2d_div_s(dest, a, v2d_norm(a));
}

FIXMATH_FUNC_ATTRS void v2d_rotate(v2d *dest, const v2d *a, fix16_t angle)
{
    fix16_t c = fix16_cos(angle);
    fix16_t s = fix16_sin(angle);
    
    dest->x = fix16_add(fix16_mul(c, a->x), fix16_mul(-s, a->y));
    dest->y = fix16_add(fix16_mul(s, a->x), fix16_mul(c, a->y));
}

/*---------------------------------------------------------------------------*/
/* 3D VECTOR OPERATIONS                                                      */
/*---------------------------------------------------------------------------*/

FIXMATH_FUNC_ATTRS void v3d_add(v3d *dest, const v3d *a, const v3d *b)
{
    dest->x = fix16_add(a->x, b->x);
    dest->y = fix16_add(a->y, b->y);
    dest->z = fix16_add(a->z, b->z);
}

FIXMATH_FUNC_ATTRS void v3d_sub(v3d *dest, const v3d *a, const v3d *b)
{
    dest->x = fix16_sub(a->x, b->x);
    dest->y = fix16_sub(a->y, b->y);
    dest->z = fix16_sub(a->z, b->z);
}

FIXMATH_FUNC_ATTRS void v3d_mul_s(v3d *dest, const v3d *a, fix16_t b)
{
    dest->x = fix16_mul(a->x, b);
    dest->y = fix16_mul(a->y, b);
    dest->z = fix16_mul(a->z, b);
}

FIXMATH_FUNC_ATTRS void v3d_div_s(v3d *dest, const v3d *a, fix16_t b)
{
    dest->x = fix16_div(a->x, b);
    dest->y = fix16_div(a->y, b);
    dest->z = fix16_div(a->z, b);
}

FIXMATH_FUNC_ATTRS fix16_t v3d_dot(const v3d *a, const v3d *b)
{
    return fa16_dot(&a->x, &a->y - &a->x, &b->x, &b->y - &b->x, 3);
}

FIXMATH_FUNC_ATTRS fix16_t v3d_norm(const v3d *a)
{
    return fa16_norm(&a->x, &a->y - &a->x, 3);
}

FIXMATH_FUNC_ATTRS void v3d_normalize(v3d *dest, const v3d *a)
{
    v3d_div_s(dest, a, v3d_norm(a));
}

FIXMATH_FUNC_ATTRS void v3d_cross(v3d *dest, const v3d *a, const v3d *b)
{
    v3d tmp;
    fa16_unalias(dest, (void**)&a, (void**)&b, &tmp, sizeof(tmp));
    
    dest->x = fix16_sub(fix16_mul(a->y, b->z), fix16_mul(a->z, b->y));
    dest->y = fix16_sub(fix16_mul(a->z, b->x), fix16_mul(a->x, b->z));
    dest->z = fix16_sub(fix16_mul(a->x, b->y), fix16_mul(a->y, b->x));
} 