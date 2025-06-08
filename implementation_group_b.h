/* Implementation - Group B: Advanced Math Functions */

/* Forward declarations for functions that may call each other */
FIXMATH_FUNC_ATTRS fix16_t fix16_div(fix16_t a, fix16_t b);
FIXMATH_FUNC_ATTRS fix16_t fix16_sin(fix16_t inAngle);
FIXMATH_FUNC_ATTRS fix16_t fix16_cos(fix16_t inAngle);

/* Utility function for bit counting leading zeros */
#if defined(__GNUC__) && !defined(FIXMATH_NO_HARD_DIVISION)
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

/* Static helper function for right shift with rounding */
static inline fix16_t fix16_rs(fix16_t x)
{
    #ifdef FIXMATH_NO_ROUNDING
        return (x >> 1);
    #else
        fix16_t y = (x >> 1) + (x & 1);
        return y;
    #endif
}

/* Static helper function for logarithms */
static fix16_t fix16__log2_inner(fix16_t x)
{
    fix16_t result = 0;
    
    while(x >= fix16_from_int(2))
    {
        result++;
        x = fix16_rs(x);
    }

    if(x == 0) return (result << 16);

    uint_fast8_t i;
    for(i = 16; i > 0; i--)
    {
        x = fix16_mul(x, x);
        result <<= 1;
        if(x >= fix16_from_int(2))
        {
            result |= 1;
            x = fix16_rs(x);
        }
    }
    #ifndef FIXMATH_NO_ROUNDING
        x = fix16_mul(x, x);
        if(x >= fix16_from_int(2)) result++;
    #endif
    
    return result;
}

/*---------------------------------------------------------------------------*/
/* DIVISION AND MODULO FUNCTIONS                                             */
/*---------------------------------------------------------------------------*/

/* Hardware division implementation (faster when available) */
#if !defined(FIXMATH_NO_HARD_DIVISION)
FIXMATH_FUNC_ATTRS fix16_t fix16_div(fix16_t a, fix16_t b)
{
    // This uses a hardware 32/32 bit division multiple times, until we have
    // computed all the bits in (a<<17)/b. Usually this takes 1-3 iterations.
    
    if (b == 0)
            return fix16_minimum;
    
    uint32_t remainder = fix_abs(a);
    uint32_t divider = fix_abs(b);
    uint64_t quotient = 0;
    int bit_pos = 17;

    // Kick-start the division a bit.
    // This improves speed in the worst-case scenarios where N and D are large
    // It gets a lower estimate for the result by N/(D >> 17 + 1).
    if (divider & 0xFFF00000)
    {
        uint32_t shifted_div = ((divider >> 17) + 1);
        quotient = remainder / shifted_div;
        uint64_t tmp = ((uint64_t)quotient * (uint64_t)divider) >> 17;
        remainder -= (uint32_t)(tmp);
    }
    
    // If the divider is divisible by 2^n, take advantage of it.
    while (!(divider & 0xF) && bit_pos >= 4)
    {
        divider >>= 4;
        bit_pos -= 4;
    }
    
    while (remainder && bit_pos >= 0)
    {
        // Shift remainder as much as we can without overflowing
        int shift = clz(remainder);
        if (shift > bit_pos) shift = bit_pos;
        remainder <<= shift;
        bit_pos -= shift;
        
        uint32_t div = remainder / divider;
        remainder = remainder % divider;
        quotient += (uint64_t)div << bit_pos;

        #ifndef FIXMATH_NO_OVERFLOW
        if (div & ~(0xFFFFFFFF >> bit_pos))
                return fix16_overflow;
        #endif
        
        remainder <<= 1;
        bit_pos--;
    }
    
    #ifndef FIXMATH_NO_ROUNDING
    // Quotient is always positive so rounding is easy
    quotient++;
    #endif
    
    fix16_t result = quotient >> 1;
    
    // Figure out the sign of the result
    if ((a ^ b) & 0x80000000)
    {
        #ifndef FIXMATH_NO_OVERFLOW
        if (result == fix16_minimum)
                return fix16_overflow;
        #endif
        
        result = -result;
    }
    
    return result;
}
#endif /* !defined(FIXMATH_NO_HARD_DIVISION) */

/* Alternative 32-bit implementation of fix16_div. Fastest on e.g. Atmel AVR.
 * This does the division manually, and is therefore good for processors that
 * do not have hardware division.
 */
#if defined(FIXMATH_NO_HARD_DIVISION)
FIXMATH_FUNC_ATTRS fix16_t fix16_div(fix16_t a, fix16_t b)
{
    // This uses the basic binary restoring division algorithm.
    // It appears to be faster to do the whole division manually than
    // trying to compose a 64-bit divide out of 32-bit divisions on
    // platforms without hardware divide.
    
    if (b == 0)
        return fix16_minimum;
    
    uint32_t remainder = fix_abs(a);
    uint32_t divider = fix_abs(b);

    uint32_t quotient = 0;
    uint32_t bit = 0x10000;
    
    /* The algorithm requires D >= R */
    while (divider < remainder)
    {
        divider <<= 1;
        bit <<= 1;
    }
    
    #ifndef FIXMATH_NO_OVERFLOW
    if (!bit)
        return fix16_overflow;
    #endif
    
    if (divider & 0x80000000)
    {
        // Perform one step manually to avoid overflows later.
        // We know that divider's bottom bit is 0 here.
        if (remainder >= divider)
        {
                quotient |= bit;
                remainder -= divider;
        }
        divider >>= 1;
        bit >>= 1;
    }
    
    /* Main division loop */
    while (bit && remainder)
    {
        if (remainder >= divider)
        {
                quotient |= bit;
                remainder -= divider;
        }
        
        remainder <<= 1;
        bit >>= 1;
    }     
            
    #ifndef FIXMATH_NO_ROUNDING
    if (remainder >= divider)
    {
        quotient++;
    }
    #endif
    
    fix16_t result = quotient;
    
    /* Figure out the sign of result */
    if ((a ^ b) & 0x80000000)
    {
        #ifndef FIXMATH_NO_OVERFLOW
        if (result == fix16_minimum)
                return fix16_overflow;
        #endif
        
        result = -result;
    }
    
    return result;
}
#endif /* defined(FIXMATH_NO_HARD_DIVISION) */

#ifndef FIXMATH_NO_OVERFLOW
/* Wrapper around fix16_div to add saturating arithmetic. */
FIXMATH_FUNC_ATTRS fix16_t fix16_sdiv(fix16_t inArg0, fix16_t inArg1)
{
    fix16_t result = fix16_div(inArg0, inArg1);
    
    if (result == fix16_overflow)
    {
        if ((inArg0 >= 0) == (inArg1 >= 0))
            return fix16_maximum;
        else
            return fix16_minimum;
    }
    
    return result;
}
#endif

FIXMATH_FUNC_ATTRS fix16_t fix16_mod(fix16_t x, fix16_t y)
{
    #ifdef FIXMATH_NO_HARD_DIVISION
        /* The reason we do this, rather than use a modulo operator
         * is that if you don't have a hardware divider, this will result
         * in faster operations when the angles are close to the bounds. 
         */
        while(x >=  y) x -= y;
        while(x <= -y) x += y;
    #else
        /* Note that in C90, the sign of result of the modulo operation is
         * undefined. in C99, it's the same as the dividend (aka numerator).
         */
        x %= y;
    #endif

    return x;
}

/*---------------------------------------------------------------------------*/
/* SQUARE ROOT FUNCTION                                                      */
/*---------------------------------------------------------------------------*/

FIXMATH_FUNC_ATTRS fix16_t fix16_sqrt(fix16_t inValue)
{
    uint8_t  neg = (inValue < 0);
    uint32_t num = (neg ? -inValue : inValue);
    uint32_t result = 0;
    uint32_t bit;
    uint8_t  n;
    
    // Many numbers will be less than 15, so
    // this gives a good balance between time spent
    // in if vs. time spent in the while loop
    // when searching for the starting value.
    if (num & 0xFFF00000)
        bit = (uint32_t)1 << 30;
    else
        bit = (uint32_t)1 << 18;
    
    while (bit > num) bit >>= 2;
    
    // The main part is executed twice, in order to avoid
    // using 64 bit values in computations.
    for (n = 0; n < 2; n++)
    {
        // First we get the top 24 bits of the answer.
        while (bit)
        {
            if (num >= result + bit)
            {
                num -= result + bit;
                result = (result >> 1) + bit;
            }
            else
            {
                result = result >> 1;
            }
            bit >>= 2;
        }
        
        if (n == 0)
        {
            // Then process it again to get the lowest 8 bits.
            if (num > 65535)
            {
                // The remainder 'num' is too large to be shifted left
                // by 16, so we have to add 1 to result manually and
                // adjust 'num' accordingly.
                // num = a - (result + 0.5)^2
                //     = num + result^2 - (result + 0.5)^2
                //     = num - result - 0.5
                num -= result;
                num = (num << 16) - 0x8000;
                result = (result << 16) + 0x8000;
            }
            else
            {
                num <<= 16;
                result <<= 16;
            }
            
            bit = 1 << 14;
        }
    }
    
    #ifndef FIXMATH_NO_ROUNDING
    // Finally, if next bit would have been 1, round the result upwards.
    if (num > result)
    {
        result++;
    }
    #endif
    
    return (neg ? -result : result);
}

/*---------------------------------------------------------------------------*/
/* INTERPOLATION FUNCTIONS                                                   */
/*---------------------------------------------------------------------------*/

FIXMATH_FUNC_ATTRS fix16_t fix16_lerp8(fix16_t inArg0, fix16_t inArg1, uint8_t inFract)
{
    int64_t tempOut = int64_mul_i32_i32(inArg0, (((int32_t)1 << 8) - inFract));
    tempOut = int64_add(tempOut, int64_mul_i32_i32(inArg1, inFract));
    tempOut = int64_shift(tempOut, -8);
    return (fix16_t)int64_lo(tempOut);
}

FIXMATH_FUNC_ATTRS fix16_t fix16_lerp16(fix16_t inArg0, fix16_t inArg1, uint16_t inFract)
{
    int64_t tempOut = int64_mul_i32_i32(inArg0, (((int32_t)1 << 16) - inFract));
    tempOut = int64_add(tempOut, int64_mul_i32_i32(inArg1, inFract));
    tempOut = int64_shift(tempOut, -16);
    return (fix16_t)int64_lo(tempOut);
}

FIXMATH_FUNC_ATTRS fix16_t fix16_lerp32(fix16_t inArg0, fix16_t inArg1, uint32_t inFract)
{
    if(inFract == 0)
        return inArg0;
    int64_t inFract64 = int64_const(0, inFract);
    int64_t subbed = int64_sub(int64_const(1,0), inFract64);
    int64_t tempOut  = int64_mul_i64_i32(subbed,  inArg0);
    tempOut    = int64_add(tempOut, int64_mul_i64_i32(inFract64, inArg1));
    return int64_hi(tempOut);
}

/*---------------------------------------------------------------------------*/
/* LOGARITHMIC/EXPONENTIAL FUNCTIONS                                         */
/*---------------------------------------------------------------------------*/

FIXMATH_FUNC_ATTRS fix16_t fix16_exp(fix16_t inValue) {
    if(inValue == 0        ) return fix16_one;
    if(inValue == fix16_one) return fix16_e;
    if(inValue >= 681391   ) return fix16_maximum;
    if(inValue <= -772243  ) return 0;

    #ifndef FIXMATH_NO_CACHE
    fix16_t tempIndex = (inValue ^ (inValue >> 4)) & 0x0FFF;
    if(_fix16_exp_cache_index[tempIndex] == inValue)
        return _fix16_exp_cache_value[tempIndex];
    #endif
                        
    /* The algorithm is based on the power series for exp(x):
     * http://en.wikipedia.org/wiki/Exponential_function#Formal_definition
     * 
     * From term n, we get term n+1 by multiplying with x/n.
     * When the sum term drops to zero, we can stop summing.
     */
            
    // The power-series converges much faster on positive values
    // and exp(-x) = 1/exp(x).
    bool neg = (inValue < 0);
    if (neg) inValue = -inValue;
            
    fix16_t result = inValue + fix16_one;
    fix16_t term = inValue;

    uint_fast8_t i;        
    for (i = 2; i < 30; i++)
    {
        term = fix16_mul(term, fix16_div(inValue, fix16_from_int(i)));
        result += term;
                
        if ((term < 500) && ((i > 15) || (term < 20)))
            break;
    }
            
    if (neg) result = fix16_div(fix16_one, result);
            
    #ifndef FIXMATH_NO_CACHE
    _fix16_exp_cache_index[tempIndex] = inValue;
    _fix16_exp_cache_value[tempIndex] = result;
    #endif

    return result;
}

FIXMATH_FUNC_ATTRS fix16_t fix16_log2(fix16_t x)
{
    // Note that a negative x gives a non-real result.
    // If x == 0, the limit of log2(x)  as x -> 0 = -infinity.
    // log2(-ve) gives a complex result.
    if (x <= 0) return fix16_overflow;

    // If the input is less than one, the result is -log2(1.0 / in)
    if (x < fix16_one)
    {
        // Note that the inverse of this would overflow.
        // This is the exact answer for log2(1.0 / 65536)
        if (x == 1) return fix16_from_int(-16);

        fix16_t inverse = fix16_div(fix16_one, x);
        return -fix16__log2_inner(inverse);
    }

    // If input >= 1, just proceed as normal.
    // Note that x == fix16_one is a special case, where the answer is 0.
    return fix16__log2_inner(x);
}

FIXMATH_FUNC_ATTRS fix16_t fix16_slog2(fix16_t x)
{
    fix16_t retval = fix16_log2(x);
    // The only overflow possible is when the input is negative.
    if(retval == fix16_overflow)
        return fix16_minimum;
    return retval;
}

FIXMATH_FUNC_ATTRS fix16_t fix16_log(fix16_t inValue)
{
    fix16_t guess = fix16_from_int(2);
    fix16_t delta;
    int scaling = 0;
    int count = 0;
    
    if (inValue <= 0)
        return fix16_minimum;
    
    // Bring the value to the most accurate range (1 < x < 100)
    const fix16_t e_to_fourth = 3578144;
    while (inValue > fix16_from_int(100))
    {
        inValue = fix16_div(inValue, e_to_fourth);
        scaling += 4;
    }
    
    while (inValue < fix16_one)
    {
        inValue = fix16_mul(inValue, e_to_fourth);
        scaling -= 4;
    }
    
    do
    {
        // Solving e(x) = y using Newton's method
        // f(x) = e(x) - y
        // f'(x) = e(x)
        fix16_t e = fix16_exp(guess);
        delta = fix16_div(inValue - e, e);
        
        // It's unlikely that logarithm is very large, so avoid overshooting.
        if (delta > fix16_from_int(3))
            delta = fix16_from_int(3);
        
        guess += delta;
    } while ((count++ < 10)
        && ((delta > 1) || (delta < -1)));
    
    return guess + fix16_from_int(scaling);
}

/*---------------------------------------------------------------------------*/
/* TRIGONOMETRIC FUNCTIONS                                                   */
/*---------------------------------------------------------------------------*/

FIXMATH_FUNC_ATTRS fix16_t fix16_sin_parabola(fix16_t inAngle)
{
    fix16_t abs_inAngle, retval;
    fix16_t mask;
    #ifndef FIXMATH_FAST_SIN
    fix16_t abs_retval;
    #endif

    /* Absolute function */
    mask = (inAngle >> (sizeof(fix16_t)*8-1));
    abs_inAngle = (inAngle + mask) ^ mask;
    
    /* On 0->PI, sin looks like x² that is :
       - centered on PI/2,
       - equals 1 on PI/2,
       - equals 0 on 0 and PI
      that means :  4/PI * x  - 4/PI² * x²
      Use abs(x) to handle (-PI) -> 0 zone.
     */
    retval = fix16_mul(FOUR_DIV_PI, inAngle) + fix16_mul( fix16_mul(_FOUR_DIV_PI2, inAngle), abs_inAngle );
    /* At this point, retval equals sin(inAngle) on important points ( -PI, -PI/2, 0, PI/2, PI),
       but is not very precise between these points
     */
    #ifndef FIXMATH_FAST_SIN
    /* Absolute value of retval */
    mask = (retval >> (sizeof(fix16_t)*8-1));
    abs_retval = (retval + mask) ^ mask;
    /* So improve its precision by adding some x^4 component to retval */
    retval += fix16_mul(X4_CORRECTION_COMPONENT, fix16_mul(retval, abs_retval) - retval );
    #endif
    return retval;
}

FIXMATH_FUNC_ATTRS fix16_t fix16_sin(fix16_t inAngle)
{
    fix16_t tempAngle = fix16_mod(inAngle, (fix16_pi << 1));

    #ifdef FIXMATH_SIN_LUT
    if(tempAngle < 0)
        tempAngle += (fix16_pi << 1);

    fix16_t tempOut;
    if(tempAngle >= fix16_pi) {
        tempAngle -= fix16_pi;
        if(tempAngle >= (fix16_pi >> 1))
            tempAngle = fix16_pi - tempAngle;
        tempOut = -(tempAngle >= _fix16_sin_lut_count ? fix16_one : _fix16_sin_lut[tempAngle]);
    } else {
        if(tempAngle >= (fix16_pi >> 1))
            tempAngle = fix16_pi - tempAngle;
        tempOut = (tempAngle >= _fix16_sin_lut_count ? fix16_one : _fix16_sin_lut[tempAngle]);
    }
    #else
    if(tempAngle > fix16_pi)
        tempAngle -= (fix16_pi << 1);
    else if(tempAngle < -fix16_pi)
        tempAngle += (fix16_pi << 1);

    #ifndef FIXMATH_NO_CACHE
    fix16_t tempIndex = ((inAngle >> 5) & 0x00000FFF);
    if(_fix16_sin_cache_index[tempIndex] == inAngle)
        return _fix16_sin_cache_value[tempIndex];
    #endif

    fix16_t tempAngleSq = fix16_mul(tempAngle, tempAngle);

    #ifndef FIXMATH_FAST_SIN // Most accurate version, accurate to ~2.1%
    fix16_t tempOut = tempAngle;
    tempAngle = fix16_mul(tempAngle, tempAngleSq);
    tempOut -= (tempAngle / 6);
    tempAngle = fix16_mul(tempAngle, tempAngleSq);
    tempOut += (tempAngle / 120);
    tempAngle = fix16_mul(tempAngle, tempAngleSq);
    tempOut -= (tempAngle / 5040);
    tempAngle = fix16_mul(tempAngle, tempAngleSq);
    tempOut += (tempAngle / 362880);
    tempAngle = fix16_mul(tempAngle, tempAngleSq);
    tempOut -= (tempAngle / 39916800);
    #else // Fast implementation, runs at 159% the speed of above 'accurate' version with an slightly lower accuracy of ~2.3%
    fix16_t tempOut;
    tempOut = fix16_mul(-13, tempAngleSq) + 546;
    tempOut = fix16_mul(tempOut, tempAngleSq) - 10923;
    tempOut = fix16_mul(tempOut, tempAngleSq) + 65536;
    tempOut = fix16_mul(tempOut, tempAngle);
    #endif

    #ifndef FIXMATH_NO_CACHE
    _fix16_sin_cache_index[tempIndex] = inAngle;
    _fix16_sin_cache_value[tempIndex] = tempOut;
    #endif
    #endif

    return tempOut;
}

FIXMATH_FUNC_ATTRS fix16_t fix16_cos(fix16_t inAngle)
{
    return fix16_sin(inAngle + (fix16_pi >> 1));
}

FIXMATH_FUNC_ATTRS fix16_t fix16_tan(fix16_t inAngle)
{
    #ifndef FIXMATH_NO_OVERFLOW
    return fix16_sdiv(fix16_sin(inAngle), fix16_cos(inAngle));
    #else
    return fix16_div(fix16_sin(inAngle), fix16_cos(inAngle));
    #endif
}

FIXMATH_FUNC_ATTRS fix16_t fix16_atan2(fix16_t inY , fix16_t inX)
{
    fix16_t abs_inY, mask, angle, r, r_3;

    #ifndef FIXMATH_NO_CACHE
    uintptr_t hash = (inX ^ inY);
    hash ^= hash >> 20;
    hash &= 0x0FFF;
    if((_fix16_atan_cache_index[0][hash] == inX) && (_fix16_atan_cache_index[1][hash] == inY))
        return _fix16_atan_cache_value[hash];
    #endif

    /* Absolute inY */
    mask = (inY >> (sizeof(fix16_t)*8-1));
    abs_inY = (inY + mask) ^ mask;

    if (inX >= 0)
    {
        r = fix16_div( (inX - abs_inY), (inX + abs_inY));
        r_3 = fix16_mul(fix16_mul(r, r),r);
        angle = fix16_mul(0x00003240 , r_3) - fix16_mul(0x0000FB50,r) + PI_DIV_4;
    } else {
        r = fix16_div( (inX + abs_inY), (abs_inY - inX));
        r_3 = fix16_mul(fix16_mul(r, r),r);
        angle = fix16_mul(0x00003240 , r_3)
            - fix16_mul(0x0000FB50,r)
            + THREE_PI_DIV_4;
    }
    if (inY < 0)
    {
        angle = -angle;
    }

    #ifndef FIXMATH_NO_CACHE
    _fix16_atan_cache_index[0][hash] = inX;
    _fix16_atan_cache_index[1][hash] = inY;
    _fix16_atan_cache_value[hash] = angle;
    #endif

    return angle;
}

FIXMATH_FUNC_ATTRS fix16_t fix16_atan(fix16_t x)
{
    return fix16_atan2(x, fix16_one);
}

FIXMATH_FUNC_ATTRS fix16_t fix16_asin(fix16_t x)
{
    if((x > fix16_one) || (x < -fix16_one))
        return 0;

    fix16_t out;
    out = (fix16_one - fix16_mul(x, x));
    out = fix16_div(x, fix16_sqrt(out));
    out = fix16_atan(out);
    return out;
}

FIXMATH_FUNC_ATTRS fix16_t fix16_acos(fix16_t x)
{
    return ((fix16_pi >> 1) - fix16_asin(x));
} 