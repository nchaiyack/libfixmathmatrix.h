/* Implementation - Group A: Core Types & Basic Arithmetic */

/* Forward declaration for fix16_mul (defined later in this group) */
FIXMATH_FUNC_ATTRS fix16_t fix16_mul(fix16_t inArg0, fix16_t inArg1);

/*---------------------------------------------------------------------------*/
/* CONVERSION FUNCTIONS                                                      */
/*---------------------------------------------------------------------------*/

FIXMATH_FUNC_ATTRS fix16_t fix16_from_int(int32_t a) { 
    return a * fix16_one; 
}

FIXMATH_FUNC_ATTRS float fix16_to_float(fix16_t a) { 
    return (float)a / fix16_one; 
}

FIXMATH_FUNC_ATTRS double fix16_to_dbl(fix16_t a) { 
    return (double)a / fix16_one; 
}

FIXMATH_FUNC_ATTRS int32_t fix16_to_int(fix16_t a)
{
#ifdef FIXMATH_NO_ROUNDING
    return (a >> 16);
#else
    if (a >= 0)
        return (a + (fix16_one >> 1)) / fix16_one;
    return (a - (fix16_one >> 1)) / fix16_one;
#endif
}

FIXMATH_FUNC_ATTRS fix16_t fix16_from_float(float a)
{
    float temp = a * fix16_one;
#ifndef FIXMATH_NO_ROUNDING
    temp += (temp >= 0) ? 0.5f : -0.5f;
#endif
    return (fix16_t)temp;
}

FIXMATH_FUNC_ATTRS fix16_t fix16_from_dbl(double a)
{
    double temp = a * fix16_one;
#ifndef FIXMATH_NO_ROUNDING
    temp += (double)((temp >= 0) ? 0.5f : -0.5f);
#endif
    return (fix16_t)temp;
}

/*---------------------------------------------------------------------------*/
/* BASIC OPERATIONS                                                          */
/*---------------------------------------------------------------------------*/

FIXMATH_FUNC_ATTRS fix16_t fix16_abs(fix16_t x) { 
    return (fix16_t)(x < 0 ? -(uint32_t)x : (uint32_t)x); 
}

FIXMATH_FUNC_ATTRS fix16_t fix16_floor(fix16_t x) { 
    return (x & 0xFFFF0000UL); 
}

FIXMATH_FUNC_ATTRS fix16_t fix16_ceil(fix16_t x) { 
    return (x & 0xFFFF0000UL) + (x & 0x0000FFFFUL ? fix16_one : 0); 
}

FIXMATH_FUNC_ATTRS fix16_t fix16_min(fix16_t x, fix16_t y) { 
    return (x < y ? x : y); 
}

FIXMATH_FUNC_ATTRS fix16_t fix16_max(fix16_t x, fix16_t y) { 
    return (x > y ? x : y); 
}

FIXMATH_FUNC_ATTRS fix16_t fix16_clamp(fix16_t x, fix16_t lo, fix16_t hi) { 
    return fix16_min(fix16_max(x, lo), hi); 
}

/*---------------------------------------------------------------------------*/
/* ANGLE AND SQUARE OPERATIONS                                               */
/*---------------------------------------------------------------------------*/

FIXMATH_FUNC_ATTRS fix16_t fix16_rad_to_deg(fix16_t radians) { 
    return fix16_mul(radians, fix16_rad_to_deg_mult); 
}

FIXMATH_FUNC_ATTRS fix16_t fix16_deg_to_rad(fix16_t degrees) { 
    return fix16_mul(degrees, fix16_deg_to_rad_mult); 
}

FIXMATH_FUNC_ATTRS fix16_t fix16_sq(fix16_t x) { 
    return fix16_mul(x, x); 
}

FIXMATH_FUNC_ATTRS int32_t fix_abs(int32_t in)
{
    if(in == fix16_minimum)
    {
        // minimum value can't be negated
        return 0x80000000;
    }
    else
    {
        return (in >= 0 ? in : -in);
    }
}

/*---------------------------------------------------------------------------*/
/* QUATERNION CONVERSION                                                     */
/*---------------------------------------------------------------------------*/

FIXMATH_FUNC_ATTRS qf16 qf16_from_v3d(v3d axis, fix16_t angle)
{
    qf16 q;
    q.a = angle;
    q.b = axis.x;
    q.c = axis.y;
    q.d = axis.z;
    return q;
}

FIXMATH_FUNC_ATTRS v3d qf16_to_v3d(qf16 qf)
{
    v3d v;
    v.x = qf.b;
    v.y = qf.c;
    v.z = qf.d;
    return v;
}

/*---------------------------------------------------------------------------*/
/* UINT32/FRACT32 FUNCTIONS                                                  */
/*---------------------------------------------------------------------------*/

FIXMATH_FUNC_ATTRS uint32_t uint32_log2(uint32_t inVal)
{
    if(inVal == 0)
        return 0;
    uint32_t tempOut = 0;
    if(inVal >= (1UL << 16)) { inVal >>= 16; tempOut += 16; }
    if(inVal >= (1 <<  8)) { inVal >>=  8; tempOut +=  8; }
    if(inVal >= (1 <<  4)) { inVal >>=  4; tempOut +=  4; }
    if(inVal >= (1 <<  2)) { inVal >>=  2; tempOut +=  2; }
    if(inVal >= (1 <<  1)) {               tempOut +=  1; }
    return tempOut;
}

FIXMATH_FUNC_ATTRS fract32_t fract32_create(uint32_t inNumerator, uint32_t inDenominator)
{
    if(inDenominator <= inNumerator)
        return 0xFFFFFFFF;
    uint32_t tempMod = (inNumerator % inDenominator);
    uint32_t tempDiv = (0xFFFFFFFF / (inDenominator - 1));
    return (tempMod * tempDiv);
}

FIXMATH_FUNC_ATTRS fract32_t fract32_invert(fract32_t inFract)
{
    return (0xFFFFFFFF - inFract);
}

#ifndef FIXMATH_NO_64BIT
FIXMATH_FUNC_ATTRS uint32_t fract32_usmul(uint32_t inVal, fract32_t inFract)
{
    return (uint32_t)(((uint64_t)inVal * (uint64_t)inFract) >> 32);
}

FIXMATH_FUNC_ATTRS int32_t fract32_smul(int32_t inVal, fract32_t inFract)
{
    if(inVal < 0)
        return -(int32_t)fract32_usmul(-inVal, inFract);
    return fract32_usmul(inVal, inFract);
}
#endif

/*---------------------------------------------------------------------------*/
/* BASIC ARITHMETIC OPERATIONS                                               */
/*---------------------------------------------------------------------------*/

/* Subtraction and addition with overflow detection.
 * The versions without overflow detection are inlined in the header.
 */
#ifndef FIXMATH_NO_OVERFLOW
FIXMATH_FUNC_ATTRS fix16_t fix16_add(fix16_t a, fix16_t b)
{
    // Use unsigned integers because overflow with signed integers is
    // an undefined operation (http://www.airs.com/blog/archives/120).
    uint32_t _a = a;
    uint32_t _b = b;
    uint32_t sum = _a + _b;

    // Overflow can only happen if sign of a == sign of b, and then
    // it causes sign of sum != sign of a.
    if (!((_a ^ _b) & 0x80000000) && ((_a ^ sum) & 0x80000000))
        return fix16_overflow;
    
    return sum;
}

FIXMATH_FUNC_ATTRS fix16_t fix16_sub(fix16_t a, fix16_t b)
{
    uint32_t _a = a;
    uint32_t _b = b;
    uint32_t diff = _a - _b;

    // Overflow can only happen if sign of a != sign of b, and then
    // it causes sign of diff != sign of a.
    if (((_a ^ _b) & 0x80000000) && ((_a ^ diff) & 0x80000000))
        return fix16_overflow;
    
    return diff;
}

/* Saturating arithmetic */
FIXMATH_FUNC_ATTRS fix16_t fix16_sadd(fix16_t a, fix16_t b)
{
    fix16_t result = fix16_add(a, b);

    if (result == fix16_overflow)
        return (a >= 0) ? fix16_maximum : fix16_minimum;

    return result;
}    

FIXMATH_FUNC_ATTRS fix16_t fix16_ssub(fix16_t a, fix16_t b)
{
    fix16_t result = fix16_sub(a, b);

    if (result == fix16_overflow)
        return (a >= 0) ? fix16_maximum : fix16_minimum;

    return result;
}

FIXMATH_FUNC_ATTRS fix16_t fix16_smul(fix16_t inArg0, fix16_t inArg1)
{
    fix16_t result = fix16_mul(inArg0, inArg1);

    if (result == fix16_overflow)
    {
        if ((inArg0 >= 0) == (inArg1 >= 0))
            return fix16_maximum;
        else
            return fix16_minimum;
    }

    return result;
}
#else

FIXMATH_FUNC_ATTRS fix16_t fix16_add(fix16_t inArg0, fix16_t inArg1) { 
    return (inArg0 + inArg1); 
}

FIXMATH_FUNC_ATTRS fix16_t fix16_sub(fix16_t inArg0, fix16_t inArg1) { 
    return (inArg0 - inArg1); 
}

#endif

/* 64-bit implementation for fix16_mul. Fastest version for e.g. ARM Cortex M3.
 * Performs a 32*32 -> 64bit multiplication. The middle 32 bits are the result,
 * bottom 16 bits are used for rounding, and upper 16 bits are used for overflow
 * detection.
 */
 
#if !defined(FIXMATH_NO_64BIT) && !defined(FIXMATH_OPTIMIZE_8BIT)
FIXMATH_FUNC_ATTRS fix16_t fix16_mul(fix16_t inArg0, fix16_t inArg1)
{
    int64_t product = (int64_t)inArg0 * inArg1;
    
    #ifndef FIXMATH_NO_OVERFLOW
    // The upper 17 bits should all be the same (the sign).
    uint32_t upper = (product >> 47);
    #endif
    
    if (product < 0)
    {
        #ifndef FIXMATH_NO_OVERFLOW
        if (~upper)
                return fix16_overflow;
        #endif
        
        #ifndef FIXMATH_NO_ROUNDING
        // This adjustment is required in order to round -1/2 correctly
        product--;
        #endif
    }
    else
    {
        #ifndef FIXMATH_NO_OVERFLOW
        if (upper)
                return fix16_overflow;
        #endif
    }
    
    #ifdef FIXMATH_NO_ROUNDING
    return product >> 16;
    #else
    fix16_t result = product >> 16;
    result += (product & 0x8000) >> 15;
    
    return result;
    #endif
}
#endif

/* 32-bit implementation of fix16_mul. Potentially fast on 16-bit processors,
 * and this is a relatively good compromise for compilers that do not support
 * uint64_t. Uses 16*16->32bit multiplications.
 */
#if defined(FIXMATH_NO_64BIT) && !defined(FIXMATH_OPTIMIZE_8BIT)
FIXMATH_FUNC_ATTRS fix16_t fix16_mul(fix16_t inArg0, fix16_t inArg1)
{
    // Each argument is divided to 16-bit parts.
    //                    AB
    //            *     CD
    // -----------
    //                    BD    16 * 16 -> 32 bit products
    //                 CB
    //                 AD
    //                AC
    //             |----| 64 bit product
    int32_t A = (inArg0 >> 16), C = (inArg1 >> 16);
    uint32_t B = (inArg0 & 0xFFFF), D = (inArg1 & 0xFFFF);
    
    int32_t AC = A*C;
    int32_t AD_CB = A*D + C*B;
    uint32_t BD = B*D;
    
    int32_t product_hi = AC + (AD_CB >> 16);
    
    // Handle carry from lower 32 bits to upper part of result.
    uint32_t ad_cb_temp = AD_CB << 16;
    uint32_t product_lo = BD + ad_cb_temp;
    if (product_lo < BD)
        product_hi++;
    
#ifndef FIXMATH_NO_OVERFLOW
    // The upper 17 bits should all be the same (the sign).
    if (product_hi >> 31 != product_hi >> 15)
        return fix16_overflow;
#endif
    
#ifdef FIXMATH_NO_ROUNDING
    return (product_hi << 16) | (product_lo >> 16);
#else
    // Subtracting 0x8000 (= 0.5) and then using signed right shift
    // achieves proper rounding to result-1, except in the corner
    // case of negative numbers and lowest word = 0x8000.
    // To handle that, we also have to subtract 1 for negative numbers.
    uint32_t product_lo_tmp = product_lo;
    product_lo -= 0x8000;
    product_lo -= (uint32_t)product_hi >> 31;
    if (product_lo > product_lo_tmp)
        product_hi--;
    
    // Discard the lowest 16 bits. Note that this is not exactly the same
    // as dividing by 0x10000. For example if product = -1, result will
    // also be -1 and not 0. This is compensated by adding +1 to the result
    // and compensating this in turn in the rounding above.
    fix16_t result = (product_hi << 16) | (product_lo >> 16);
    result += 1;
    return result;
#endif
}
#endif

/* 8-bit implementation of fix16_mul. Fastest on e.g. Atmel AVR.
 * Uses 8*8->16bit multiplications, and also skips any bytes that
 * are zero.
 */
#if defined(FIXMATH_OPTIMIZE_8BIT)
FIXMATH_FUNC_ATTRS fix16_t fix16_mul(fix16_t inArg0, fix16_t inArg1)
{
    uint32_t _a = fix_abs(inArg0);
    uint32_t _b = fix_abs(inArg1);
    
    uint8_t va[4] = {_a, (_a >> 8), (_a >> 16), (_a >> 24)};
    uint8_t vb[4] = {_b, (_b >> 8), (_b >> 16), (_b >> 24)};
    
    uint32_t low = 0;
    uint32_t mid = 0;
    
    // Result column i depends on va[0..i] and vb[i..0]

    #ifndef FIXMATH_NO_OVERFLOW
    // i = 6
    if (va[3] && vb[3]) return fix16_overflow;
    #endif
    
    // i = 5
    if (va[2] && vb[3]) mid += (uint16_t)va[2] * vb[3];
    if (va[3] && vb[2]) mid += (uint16_t)va[3] * vb[2];
    mid <<= 8;
    
    // i = 4
    if (va[1] && vb[3]) mid += (uint16_t)va[1] * vb[3];
    if (va[2] && vb[2]) mid += (uint16_t)va[2] * vb[2];
    if (va[3] && vb[1]) mid += (uint16_t)va[3] * vb[1];
    
    #ifndef FIXMATH_NO_OVERFLOW
    if (mid & 0xFF000000) return fix16_overflow;
    #endif
    mid <<= 8;
    
    // i = 3
    if (va[0] && vb[3]) mid += (uint16_t)va[0] * vb[3];
    if (va[1] && vb[2]) mid += (uint16_t)va[1] * vb[2];
    if (va[2] && vb[1]) mid += (uint16_t)va[2] * vb[1];
    if (va[3] && vb[0]) mid += (uint16_t)va[3] * vb[0];
    
    #ifndef FIXMATH_NO_ROUNDING
    low += 0x8000;
    #endif
    mid += (low >> 16);
    
    #ifndef FIXMATH_NO_OVERFLOW
    if (mid & 0x80000000) return fix16_overflow;
    #endif
    
    fix16_t result = mid;
    
    /* Figure out the sign of result */
    if ((inArg0 ^ inArg1) & 0x80000000)
    {
        #ifndef FIXMATH_NO_OVERFLOW
        if (result == fix16_minimum)
            return fix16_overflow;
        #endif
        
        result = -result;
    }
    
    return result;
}
#endif 