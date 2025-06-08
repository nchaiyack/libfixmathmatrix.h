/* EXTRACTED IMPLEMENTATIONS */
static inline double  fix16_to_dbl(fix16_t a)   { return (double)a / fix16_one; }

static inline int fix16_to_int(fix16_t a)
{
#ifdef FIXMATH_NO_ROUNDING
    return (a >> 16);
#else
	if (a >= 0)
		return (a + (fix16_one >> 1)) / fix16_one;
	return (a - (fix16_one >> 1)) / fix16_one;
#endif
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE uint32_t uint32_log2(uint32_t inVal)
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
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE fract32_t fract32_create(uint32_t inNumerator, uint32_t inDenominator)
{
	if(inDenominator <= inNumerator)
		return 0xFFFFFFFF;
	uint32_t tempMod = (inNumerator % inDenominator);
	uint32_t tempDiv = (0xFFFFFFFF / (inDenominator - 1));
	return (tempMod * tempDiv);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE fract32_t fract32_invert(fract32_t inFract)
{
	return (0xFFFFFFFF - inFract);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE uint32_t fract32_usmul(uint32_t inVal, fract32_t inFract)
{
	return (uint32_t)(((uint64_t)inVal * (uint64_t)inFract) >> 32);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE int32_t fract32_smul(int32_t inVal, fract32_t inFract)
{
	if(inVal < 0)
        return -(int32_t)fract32_usmul(-inVal, inFract);
	return fract32_usmul(inVal, inFract);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_add(fix16_t a, fix16_t b)
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
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_sub(fix16_t a, fix16_t b)
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
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_sadd(fix16_t a, fix16_t b)
{
	fix16_t result = fix16_add(a, b);

	if (result == fix16_overflow)
		return (a >= 0) ? fix16_maximum : fix16_minimum;

	return result;
}	
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_ssub(fix16_t a, fix16_t b)
{
	fix16_t result = fix16_sub(a, b);

	if (result == fix16_overflow)
		return (a >= 0) ? fix16_maximum : fix16_minimum;

	return result;
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_mul(fix16_t inArg0, fix16_t inArg1)
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
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_mul(fix16_t inArg0, fix16_t inArg1)
{
	// Each argument is divided to 16-bit parts.
	//					AB
	//			*	 CD
	// -----------
	//					BD	16 * 16 -> 32 bit products
	//				 CB
	//				 AD
	//				AC
	//			 |----| 64 bit product
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
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_mul(fix16_t inArg0, fix16_t inArg1)
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
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_smul(fix16_t inArg0, fix16_t inArg1)
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
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_div(fix16_t a, fix16_t b)
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
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_div(fix16_t a, fix16_t b)
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
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_sdiv(fix16_t inArg0, fix16_t inArg1)
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
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_mod(fix16_t x, fix16_t y)
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
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_lerp8(fix16_t inArg0, fix16_t inArg1, uint8_t inFract)
{
	int64_t tempOut = int64_mul_i32_i32(inArg0, (((int32_t)1 << 8) - inFract));
	tempOut = int64_add(tempOut, int64_mul_i32_i32(inArg1, inFract));
	tempOut = int64_shift(tempOut, -8);
	return (fix16_t)int64_lo(tempOut);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_lerp16(fix16_t inArg0, fix16_t inArg1, uint16_t inFract)
{
	int64_t tempOut = int64_mul_i32_i32(inArg0, (((int32_t)1 << 16) - inFract));
	tempOut = int64_add(tempOut, int64_mul_i32_i32(inArg1, inFract));
	tempOut = int64_shift(tempOut, -16);
	return (fix16_t)int64_lo(tempOut);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_lerp32(fix16_t inArg0, fix16_t inArg1, uint32_t inFract)
{
	if(inFract == 0)
		return inArg0;
	int64_t inFract64 = int64_const(0, inFract);
	int64_t subbed = int64_sub(int64_const(1,0), inFract64);
	int64_t tempOut  = int64_mul_i64_i32(subbed,  inArg0);
	tempOut	= int64_add(tempOut, int64_mul_i64_i32(inFract64, inArg1));
	return int64_hi(tempOut);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_sqrt(fix16_t inValue)
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
				//	 = num + result^2 - (result + 0.5)^2
				//	 = num - result - 0.5
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
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_rs(fix16_t x)
{
	#ifdef FIXMATH_NO_ROUNDING
		return (x >> 1);
	#else
		fix16_t y = (x >> 1) + (x & 1);
		return y;
	#endif
}
FIXMATHMATRIX_DEF fix16_t fix16_exp(fix16_t inValue) {
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
FIXMATHMATRIX_DEF fix16_t fix16_log2(fix16_t x)
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
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_slog2(fix16_t x)
{
	fix16_t retval = fix16_log2(x);
	// The only overflow possible is when the input is negative.
	if(retval == fix16_overflow)
		return fix16_minimum;
	return retval;
}
FIXMATHMATRIX_DEF fix16_t fix16_log(fix16_t inValue)
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
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_sin_parabola(fix16_t inAngle)
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
FIXMATHMATRIX_DEF fix16_t fix16_sin(fix16_t inAngle)
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
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_cos(fix16_t inAngle)
{
	return fix16_sin(inAngle + (fix16_pi >> 1));
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_tan(fix16_t inAngle)
{
	#ifndef FIXMATH_NO_OVERFLOW
	return fix16_sdiv(fix16_sin(inAngle), fix16_cos(inAngle));
	#else
	return fix16_div(fix16_sin(inAngle), fix16_cos(inAngle));
	#endif
}
FIXMATHMATRIX_DEF fix16_t fix16_atan2(fix16_t inY , fix16_t inX)
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
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_atan(fix16_t x)
{
	return fix16_atan2(x, fix16_one);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_asin(fix16_t x)
{
	if((x > fix16_one) || (x < -fix16_one))
		return 0;

	fix16_t out;
	out = (fix16_one - fix16_mul(x, x));
	out = fix16_div(x, fix16_sqrt(out));
	out = fix16_atan(out);
	return out;
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE FIXMATH_FUNC_ATTRS fix16_t fix16_acos(fix16_t x)
{
	return ((fix16_pi >> 1) - fix16_asin(x));
}
FIXMATHMATRIX_DEF void fix16_to_str(fix16_t value, char *buf, int decimals)
{
    uint32_t uvalue = (value >= 0) ? value : -value;
    if (value < 0)
        *buf++ = '-';

    /* Separate the integer and decimal parts of the value */
    unsigned intpart = uvalue >> 16;
    uint32_t fracpart = uvalue & 0xFFFF;
    uint32_t scale = scales[decimals & 7];
    fracpart = fix16_mul(fracpart, scale);
    
    if (fracpart >= scale)
    {
        /* Handle carry from decimal part */
        intpart++;
        fracpart -= scale;    
    }
    
    /* Format integer part */
    buf = itoa_loop(buf, 10000, intpart, true);
    
    /* Format decimal part (if any) */
    if (scale != 1)
    {
        *buf++ = '.';
        buf = itoa_loop(buf, scale / 10, fracpart, false);
    }
    
    *buf = '\0';
}
FIXMATHMATRIX_DEF fix16_t fix16_from_str(const char *buf)
{
    while (isspace((unsigned char) *buf))
        buf++;
    
    /* Decode the sign */
    bool negative = (*buf == '-');
    if (*buf == '+' || *buf == '-')
        buf++;

    /* Decode the integer part */
    uint32_t intpart = 0;
    int count = 0;
    while (isdigit((unsigned char) *buf))
    {
        intpart *= 10;
        intpart += *buf++ - '0';
        count++;
    }
    
    #ifdef FIXMATH_NO_OVERFLOW
    if (count == 0)
        return fix16_overflow;
    #else
    if (count == 0 || count > 5
        || intpart > 32768 || (!negative && intpart > 32767))
        return fix16_overflow;
    #endif
    
    fix16_t value = intpart << 16;
    
    /* Decode the decimal part */
    if (*buf == '.' || *buf == ',')
    {
        buf++;
        
        uint32_t fracpart = 0;
        uint32_t scale = 1;
        while (isdigit((unsigned char) *buf) && scale < 100000)
        {
            scale *= 10;
            fracpart *= 10;
            fracpart += *buf++ - '0';
        }
        
        value += fix16_div(fracpart, scale);
    }
    
    /* Verify that there is no garbage left over */
    while (*buf != '\0')
    {
        if (!isdigit((unsigned char) *buf) && !isspace((unsigned char) *buf))
            return fix16_overflow;
        
        buf++;
    }
    
    return negative ? -value : value;
}
FIXMATHMATRIX_DEF void fix16_fft(INPUT_TYPE *input, fix16_t *real, fix16_t *imag, unsigned transform_length)
{
    int log_length = ilog2(transform_length);
    transform_length = 1 << log_length;

    unsigned i;
    for (i = 0; i < transform_length / 4; i++)
    {
        four_point_dft(input + INPUT_INDEX(rbit_n(i, log_length - 2)), transform_length / 4, real + 4*i, imag + 4*i);
    }

    for (i = 2; i < (unsigned) log_length; i++)
    {
        butterfly(real, imag, 1 << i, transform_length / (2 << i));
    }
    
#ifdef OUTPUT_SCALE
    fix16_t scale = OUTPUT_SCALE(transform_length);
    for (i = 0; i < transform_length; i++)
    {
        real[i] = fix16_mul(real[i], scale);
        imag[i] = fix16_mul(imag[i], scale);
    }
#endif
}
FIXMATHMATRIX_DEF void fa16_unalias(void *dest, void **a, void **b, void *tmp, unsigned size)
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
FIXMATHMATRIX_DEF fix16_t fa16_norm(const fix16_t *a, uint_fast8_t a_stride, uint_fast8_t n)
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
FIXMATHMATRIX_DEF fix16_t fa16_norm(const fix16_t *a, uint_fast8_t a_stride, uint_fast8_t n)
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
FIXMATHMATRIX_DEF void v2d_add(v2d *dest, const v2d *a, const v2d *b)
{
    dest->x = fix16_add(a->x, b->x);
    dest->y = fix16_add(a->y, b->y);
}
FIXMATHMATRIX_DEF void v2d_sub(v2d *dest, const v2d *a, const v2d *b)
{
    dest->x = fix16_sub(a->x, b->x);
    dest->y = fix16_sub(a->y, b->y);
}
FIXMATHMATRIX_DEF void v2d_mul_s(v2d *dest, const v2d *a, fix16_t b)
{
    dest->x = fix16_mul(a->x, b);
    dest->y = fix16_mul(a->y, b);
}
FIXMATHMATRIX_DEF void v2d_div_s(v2d *dest, const v2d *a, fix16_t b)
{
    dest->x = fix16_div(a->x, b);
    dest->y = fix16_div(a->y, b);
}
FIXMATHMATRIX_DEF fix16_t v2d_dot(const v2d *a, const v2d *b)
{
    return fix16_add(fix16_mul(a->x, b->x), fix16_mul(a->y, b->y));
}
FIXMATHMATRIX_DEF fix16_t v2d_norm(const v2d *a)
{
    return fa16_norm(&a->x, &a->y - &a->x, 2);
}
FIXMATHMATRIX_DEF void v2d_normalize(v2d *dest, const v2d *a)
{
    v2d_div_s(dest, a, v2d_norm(a));
}
FIXMATHMATRIX_DEF void v2d_rotate(v2d *dest, const v2d *a, fix16_t angle)
{
    fix16_t c = fix16_cos(angle);
    fix16_t s = fix16_sin(angle);
    
    dest->x = fix16_add(fix16_mul(c, a->x), fix16_mul(-s, a->y));
    dest->y = fix16_add(fix16_mul(s, a->x), fix16_mul(c, a->y));
}
FIXMATHMATRIX_DEF void v3d_add(v3d *dest, const v3d *a, const v3d *b)
{
    dest->x = fix16_add(a->x, b->x);
    dest->y = fix16_add(a->y, b->y);
    dest->z = fix16_add(a->z, b->z);
}
FIXMATHMATRIX_DEF void v3d_sub(v3d *dest, const v3d *a, const v3d *b)
{
    dest->x = fix16_sub(a->x, b->x);
    dest->y = fix16_sub(a->y, b->y);
    dest->z = fix16_sub(a->z, b->z);
}
FIXMATHMATRIX_DEF void v3d_mul_s(v3d *dest, const v3d *a, fix16_t b)
{
    dest->x = fix16_mul(a->x, b);
    dest->y = fix16_mul(a->y, b);
    dest->z = fix16_mul(a->z, b);
}
FIXMATHMATRIX_DEF void v3d_div_s(v3d *dest, const v3d *a, fix16_t b)
{
    dest->x = fix16_div(a->x, b);
    dest->y = fix16_div(a->y, b);
    dest->z = fix16_div(a->z, b);
}
FIXMATHMATRIX_DEF fix16_t v3d_dot(const v3d *a, const v3d *b)
{
    return fa16_dot(&a->x, &a->y - &a->x, &b->x, &b->y - &b->x, 3);
}
FIXMATHMATRIX_DEF fix16_t v3d_norm(const v3d *a)
{
    return fa16_norm(&a->x, &a->y - &a->x, 3);
}
FIXMATHMATRIX_DEF void v3d_normalize(v3d *dest, const v3d *a)
{
    v3d_div_s(dest, a, v3d_norm(a));
}
FIXMATHMATRIX_DEF void v3d_cross(v3d *dest, const v3d *a, const v3d *b)
{
    v3d tmp;
    fa16_unalias(dest, (void**)&a, (void**)&b, &tmp, sizeof(tmp));
    
    dest->x = fix16_sub(fix16_mul(a->y, b->z), fix16_mul(a->z, b->y));
    dest->y = fix16_sub(fix16_mul(a->z, b->x), fix16_mul(a->x, b->z));
    dest->z = fix16_sub(fix16_mul(a->x, b->y), fix16_mul(a->y, b->x));
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void mf16_fill(mf16 *dest, fix16_t value)
{
    int row, column;
    dest->errors = 0;
    
    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            dest->data[row][column] = value;
        }
    }
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void mf16_fill_diagonal(mf16 *dest, fix16_t value)
{
    int row;
    
    mf16_fill(dest, 0);
    
    for (row = 0; row < dest->rows; row++)
    {
        dest->data[row][row] = value;
    }
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void mf16_mul(mf16 *dest, const mf16 *a, const mf16 *b)
{
    int row, column;
    
    // If dest and input matrices alias, we have to use a temp matrix.
    mf16 tmp;
    fa16_unalias(dest, (void**)&a, (void**)&b, &tmp, sizeof(tmp));
    
    dest->errors = a->errors | b->errors;
    
    if (a->columns != b->rows)
        dest->errors |= FIXMATRIX_DIMERR;
    
    dest->rows = a->rows;
    dest->columns = b->columns;
    
    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            dest->data[row][column] = fa16_dot(
                &a->data[row][0], 1,
                &b->data[0][column], FIXMATRIX_MAX_SIZE,
                a->columns);
            
            if (dest->data[row][column] == fix16_overflow)
                dest->errors |= FIXMATRIX_OVERFLOW;
        }
    }
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void mf16_mul_at(mf16 *dest, const mf16 *at, const mf16 *b)
{
    int row, column;
    
    // If dest and input matrices alias, we have to use a temp matrix.
    mf16 tmp;
    fa16_unalias(dest, (void**)&at, (void**)&b, &tmp, sizeof(tmp));
    
    dest->errors = at->errors | b->errors;
    
    if (at->rows != b->rows)
        dest->errors |= FIXMATRIX_DIMERR;
    
    dest->rows = at->columns;
    dest->columns = b->columns;
    
    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            dest->data[row][column] = fa16_dot(
                &at->data[0][row], FIXMATRIX_MAX_SIZE,
                &b->data[0][column], FIXMATRIX_MAX_SIZE,
                at->rows);
            
            if (dest->data[row][column] == fix16_overflow)
                dest->errors |= FIXMATRIX_OVERFLOW;
        }
    }
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void mf16_mul_bt(mf16 *dest, const mf16 *a, const mf16 *bt)
{
    int row, column;
    
    // If dest and input matrices alias, we have to use a temp matrix.
    mf16 tmp;
    fa16_unalias(dest, (void**)&a, (void**)&bt, &tmp, sizeof(tmp));
    
    dest->errors = a->errors | bt->errors;
    
    if (a->columns != bt->columns)
        dest->errors |= FIXMATRIX_DIMERR;
    
    dest->rows = a->rows;
    dest->columns = bt->rows;
    
    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            dest->data[row][column] = fa16_dot(
                &a->data[row][0], 1,
                &bt->data[column][0], 1,
                a->columns);
            
            if (dest->data[row][column] == fix16_overflow)
                dest->errors |= FIXMATRIX_OVERFLOW;
        }
    }
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void mf16_add(mf16 *dest, const mf16 *a, const mf16 *b)
{
    mf16_addsub(dest, a, b, 1);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void mf16_sub(mf16 *dest, const mf16 *a, const mf16 *b)
{
    mf16_addsub(dest, a, b, 0);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void mf16_transpose(mf16 *dest, const mf16 *matrix)
{
    int row, column;
    
    // This code is a bit tricky in order to work
    // in the situation when dest = matrix.
    // Before writing a value in dest, we must copy
    // the corresponding value from matrix to a temporary
    // variable.
    
    // We actually transpose a n by n square matrix, because
    // that can be done in-place easily. Because mf16 always
    // allocates a square area even if actual matrix is smaller,
    // this is not a problem.
    int n = matrix->rows;
    if (matrix->columns > n) n = matrix->columns;
    
    uint8_t rows = matrix->rows;
    dest->rows = matrix->columns;
    dest->columns = rows;
    dest->errors = matrix->errors;
    
    for (row = 0; row < n; row++)
    {
        for (column = 0; column < row; column++)
        {
            fix16_t temp = matrix->data[row][column];
            dest->data[row][column] = matrix->data[column][row];
            dest->data[column][row] = temp;
        }
        
        dest->data[row][row] = matrix->data[row][row];
    }
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void mf16_mul_s(mf16 *dest, const mf16 *matrix, fix16_t scalar)
{
    mf16_divmul_s(dest, matrix, scalar, 1);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void mf16_div_s(mf16 *dest, const mf16 *matrix, fix16_t scalar)
{
    mf16_divmul_s(dest, matrix, scalar, 0);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void mf16_qr_decomposition(mf16 *q, mf16 *r, const mf16 *matrix, int reorthogonalize)
{
    int i, j, reorth;
    fix16_t dot, norm;
    
    uint8_t stride = FIXMATRIX_MAX_SIZE;
    uint8_t n = matrix->rows;
    
    // This uses the modified Gram-Schmidt algorithm.
    // subtract_projection takes advantage of the fact that
    // previous columns have already been normalized.
    
    // We start with q = matrix
    if (q != matrix)
    {
        *q = *matrix;
    }
    
    // R is initialized to have square size of cols(A) and zeroed.
    r->columns = matrix->columns;
    r->rows = matrix->columns;
    r->errors = 0;
    mf16_fill(r, 0);
    
    // Now do the actual Gram-Schmidt for the rows.
    for (j = 0; j < q->columns; j++)
    {
        for (reorth = 0; reorth <= reorthogonalize; reorth++)
        {
            for (i = 0; i < j; i++)
            {
                fix16_t *v = &q->data[0][j];
                fix16_t *u = &q->data[0][i];
                
                dot = fa16_dot(v, stride, u, stride, n);
                subtract_projection(v, u, dot, n, &q->errors);
                
                if (dot == fix16_overflow)
                    q->errors |= FIXMATRIX_OVERFLOW;
                
                r->data[i][j] += dot;
            }
        }
        
        // Normalize the row in q
        norm = fa16_norm(&q->data[0][j], stride, n);
        r->data[j][j] = norm;
        
        if (norm == fix16_overflow)
            q->errors |= FIXMATRIX_OVERFLOW;
        
        if (norm < 5 && norm > -5)
        {
            // Nearly zero norm, which means that the row
            // was linearly dependent.
            q->errors |= FIXMATRIX_SINGULAR;
            continue;
        }
        
        for (i = 0; i < n; i++)
        {
            // norm >= v[i] for all i, therefore this division
            // doesn't overflow unless norm approaches 0.
            q->data[i][j] = fix16_div(q->data[i][j], norm);
        }
    }
    
    r->errors = q->errors;
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void mf16_solve(mf16 *dest, const mf16 *q, const mf16 *r, const mf16 *matrix)
{
    int row, column, variable;
    
    if (r->columns != r->rows || r->columns != q->columns || r == dest)
    {
        dest->errors |= FIXMATRIX_USEERR;
        return;
    }
    
    // Ax=b <=> QRx=b <=> Q'QRx=Q'b <=> Rx=Q'b
    // Q'b is calculated directly and x is then solved row-by-row.
    mf16_mul_at(dest, q, matrix);
    
    for (column = 0; column < dest->columns; column++)
    {
        for (row = dest->rows - 1; row >= 0; row--)
        {
            fix16_t value = dest->data[row][column];
            
            // Subtract any already solved variables
            for (variable = row + 1; variable < r->columns; variable++)
            {
                fix16_t multiplier = r->data[row][variable];
                fix16_t known_value = dest->data[variable][column];
                fix16_t product = fix16_mul(multiplier, known_value);
                value = fix16_sub(value, product);
                
                if (product == fix16_overflow || value == fix16_overflow)
                {
                    dest->errors |= FIXMATRIX_OVERFLOW;
                }
            }
            
            // Now value = R_ij x_i <=> x_i = value / R_ij
            fix16_t divider = r->data[row][row];
            if (divider == 0)
            {
                dest->errors |= FIXMATRIX_SINGULAR;
                dest->data[row][column] = 0;
                continue;
            }
            
            fix16_t result = fix16_div(value, divider);
            dest->data[row][column] = result;
            
            if (result == fix16_overflow)
            {
                dest->errors |= FIXMATRIX_OVERFLOW;
            }
        }
    }
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void mf16_cholesky(mf16 *dest, const mf16 *matrix)
{
    // This is the Cholesky–Banachiewicz algorithm.
    // Refer to http://en.wikipedia.org/wiki/Cholesky_decomposition#The_Cholesky.E2.80.93Banachiewicz_and_Cholesky.E2.80.93Crout_algorithms
    
    int row, column, k;
    dest->errors = matrix->errors;
    
    if (matrix->rows != matrix->columns)
        dest->errors |= FIXMATRIX_DIMERR;
    
    dest->rows = dest->columns = matrix->rows;
    
    for (row = 0; row < dest->rows; row++)
    {
        for (column = 0; column < dest->columns; column++)
        {
            if (row == column)
            {
                // Value on the diagonal
                // Ljj = sqrt(Ajj - sum(Ljk^2, k = 1..(j-1))
                fix16_t value = matrix->data[row][column];
                for (k = 0; k < column; k++)
                {
                    fix16_t Ljk = dest->data[row][k];
                    Ljk = fix16_mul(Ljk, Ljk);
                    value = fix16_sub(value, Ljk);
                    
                    if (value == fix16_overflow || Ljk == fix16_overflow)
                        dest->errors |= FIXMATRIX_OVERFLOW;
                }
                
                if (value < 0)
                {
                    if (value < -65)
                        dest->errors |= FIXMATRIX_NEGATIVE;
                    value = 0;
                }
                
                dest->data[row][column] = fix16_sqrt(value);
            }
            else if (row < column)
            {
                // Value above diagonal
                dest->data[row][column] = 0;
            }
            else
            {
                // Value below diagonal
                // Lij = 1/Ljj (Aij - sum(Lik Ljk, k = 1..(j-1)))
                fix16_t value = matrix->data[row][column];
                for (k = 0; k < column; k++)
                {
                    fix16_t Lik = dest->data[row][k];
                    fix16_t Ljk = dest->data[column][k];
                    fix16_t product = fix16_mul(Lik, Ljk);
                    value = fix16_sub(value, product);
                    
                    if (value == fix16_overflow || product == fix16_overflow)
                        dest->errors |= FIXMATRIX_OVERFLOW;
                }
                fix16_t Ljj = dest->data[column][column];
                value = fix16_div(value, Ljj);
                dest->data[row][column] = value;
                
                if (value == fix16_overflow)
                    dest->errors |= FIXMATRIX_OVERFLOW;
            }
        }
    }
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void mf16_invert_lt(mf16 *dest, const mf16 *matrix)
{
    // This is port of the algorithm as found in the Efficient Java Matrix Library
    // https://code.google.com/p/efficient-java-matrix-library

    int_fast8_t i, j, k;
    const uint_fast8_t n = matrix->rows;

    // If dest and input matrices alias, we have to use a temp matrix.
    mf16 tmp;
    fa16_unalias(dest, (void**)&matrix, (void**)&matrix, &tmp, sizeof(tmp));

    dest->errors = dest->errors | matrix->errors;

    // TODO reorder these operations to avoid cache misses

    // inverts the lower triangular system and saves the result
    // in the upper triangle to minimize cache misses
    for (i = 0; i < n; ++i)
    {
        const fix16_t el_ii = matrix->data[i][i];
        for (j = 0; j <= i; ++j)
        {
            fix16_t sum = (i == j) ? fix16_one : 0;
            for (k = i - 1; k >= j; --k)
            {
                sum = fix16_sub(sum, fix16_mul(matrix->data[i][k], dest->data[j][k]));
            }
            dest->data[j][i] = fix16_div(sum, el_ii);
        }
    }
    // solve the system and handle the previous solution being in the upper triangle
    // takes advantage of symmetry
    for (i = n - 1; i >= 0; --i)
    {
        const fix16_t el_ii = matrix->data[i][i];
        for (j = 0; j <= i; ++j)
        {
            fix16_t sum = (i < j) ? 0 : dest->data[j][i];
            for (k = i + 1; k < n; ++k)
            {
                sum = fix16_sub(sum, fix16_mul(matrix->data[k][i], dest->data[j][k]));
            }
            dest->data[i][j] = dest->data[j][i] = fix16_div(sum, el_ii);
        }
    }
} 
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void qf16_conj(qf16 *dest, const qf16 *q)
{
    dest->a = q->a;
    dest->b = - q->b;
    dest->c = - q->c;
    dest->d = - q->d;
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void qf16_add(qf16 *dest, const qf16 *q, const qf16 *r)
{
    dest->a = q->a + r->a;
    dest->b = q->b + r->b;
    dest->c = q->c + r->c;
    dest->d = q->d + r->d;
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void qf16_mul_s(qf16 *dest, const qf16 *q, fix16_t s)
{
    dest->a = fix16_mul(q->a, s);
    dest->b = fix16_mul(q->b, s);
    dest->c = fix16_mul(q->c, s);
    dest->d = fix16_mul(q->d, s);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void qf16_div_s(qf16 *dest, const qf16 *q, fix16_t s)
{
    dest->a = fix16_div(q->a, s);
    dest->b = fix16_div(q->b, s);
    dest->c = fix16_div(q->c, s);
    dest->d = fix16_div(q->d, s);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE fix16_t qf16_dot(const qf16 *q, const qf16 *r)
{
    return fa16_dot(&q->a, &q->b - &q->a, &r->a, &r->b - &r->a, 4);    
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE fix16_t qf16_norm(const qf16 *q)
{
    return fa16_norm(&q->a, &q->b - &q->a, 4);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void qf16_normalize(qf16 *dest, const qf16 *q)
{
    qf16_div_s(dest, q, qf16_norm(q));
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void qf16_mul(qf16 *dest, const qf16 *q, const qf16 *r)
{
    qf16 tmp;
    fa16_unalias(dest, (void**)&q, (void**)&r, &tmp, sizeof(tmp));
    
    dest->a = fix16_mul(q->a, r->a) - fix16_mul(q->b, r->b) - fix16_mul(q->c, r->c) - fix16_mul(q->d, r->d);
    dest->b = fix16_mul(q->a, r->b) + fix16_mul(q->b, r->a) + fix16_mul(q->c, r->d) - fix16_mul(q->d, r->c);
    dest->c = fix16_mul(q->a, r->c) - fix16_mul(q->b, r->d) + fix16_mul(q->c, r->a) + fix16_mul(q->d, r->b);
    dest->d = fix16_mul(q->a, r->d) + fix16_mul(q->b, r->c) - fix16_mul(q->c, r->b) + fix16_mul(q->d, r->a);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void qf16_pow(qf16 *dest, const qf16 *q, fix16_t power)
{
    fix16_t old_half_angle = fix16_acos(q->a);
    fix16_t new_half_angle = fix16_mul(old_half_angle, power);
    fix16_t multiplier = 0;
    
    if (old_half_angle > 10) // Guard against almost-zero divider
    {
        multiplier = fix16_div(fix16_sin(new_half_angle),
                               fix16_sin(old_half_angle));
    }
    
    dest->a = fix16_cos(new_half_angle);
    dest->b = fix16_mul(q->b, multiplier);
    dest->c = fix16_mul(q->c, multiplier);
    dest->d = fix16_mul(q->d, multiplier);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void qf16_avg(qf16 *dest, const qf16 *q1, const qf16 *q2, fix16_t weight)
{
    // z = sqrt((w1 - w2)^2 + 4 w1 w2 (q1' q2)^2
    // <=>
    // z = sqrt((2 w1 - 1)^2 + 4 w1 (1 - w1) (q1' q2)^2)
    fix16_t dot = qf16_dot(q1, q2);
    fix16_t z = fix16_sq(2 * weight - fix16_one)
            + fix16_mul(4 * weight, fix16_mul((fix16_one - weight), fix16_sq(dot)));
    z = fix16_sqrt(z);
    
    // q = 2 * w1 * (q1' q2) q1 + (w2 - w1 + z) q2
    // <=>
    // q = 2 * w1 * (q1' q2) q1 + (1 - 2 * w1 + z) q2
    qf16 tmp1;
    qf16_mul_s(&tmp1, q1, fix16_mul(2 * weight, dot));
    
    qf16 tmp2;
    qf16_mul_s(&tmp2, q2, fix16_one - 2 * weight + z);
    
    qf16_add(dest, &tmp1, &tmp2);
    qf16_normalize(dest, dest);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void qf16_from_axis_angle(qf16 *dest, const v3d *axis, fix16_t angle)
{
    angle /= 2;
    fix16_t scale = fix16_sin(angle);
    
    dest->a = fix16_cos(angle);
    dest->b = fix16_mul(axis->x, scale);
    dest->c = fix16_mul(axis->y, scale);
    dest->d = fix16_mul(axis->z, scale);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void qf16_to_matrix(mf16 *dest, const qf16 *q)
{
    dest->rows = dest->columns = 3;
    dest->errors = 0;
    dest->data[0][0] = fix16_one - 2 * (fix16_sq(q->c) + fix16_sq(q->d));
    dest->data[1][1] = fix16_one - 2 * (fix16_sq(q->b) + fix16_sq(q->d));
    dest->data[2][2] = fix16_one - 2 * (fix16_sq(q->b) + fix16_sq(q->c));
    
    dest->data[1][0] = 2 * (fix16_mul(q->b, q->c) + fix16_mul(q->a, q->d));
    dest->data[0][1] = 2 * (fix16_mul(q->b, q->c) - fix16_mul(q->a, q->d));
    
    dest->data[2][0] = 2 * (fix16_mul(q->b, q->d) - fix16_mul(q->a, q->c));
    dest->data[0][2] = 2 * (fix16_mul(q->b, q->d) + fix16_mul(q->a, q->c));
    
    dest->data[2][1] = 2 * (fix16_mul(q->c, q->d) + fix16_mul(q->a, q->b));
    dest->data[1][2] = 2 * (fix16_mul(q->c, q->d) - fix16_mul(q->a, q->b));
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void print_fix16_t(FILE *stream, fix16_t value, uint_fast8_t width, uint_fast8_t decimals)
{
    char buf[13];
    fix16_to_str(value, buf, decimals);
    
    uint_fast8_t len = strlen(buf);
    if (len < width)
    {
        width -= len;
        while (width-- > 0)
            fputc(' ', stream);
    }

    fputs(buf, stream);
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void print_v2d(FILE *stream, const v2d *vector)
{
    fprintf(stream, "(");
    print_fix16_t(stream, vector->x, 9, 4);
    fprintf(stream, ", ");
    print_fix16_t(stream, vector->y, 9, 4);
    fprintf(stream, ")");
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void print_v3d(FILE *stream, const v3d *vector)
{
    fprintf(stream, "(");
    print_fix16_t(stream, vector->x, 9, 4);
    fprintf(stream, ", ");
    print_fix16_t(stream, vector->y, 9, 4);
    fprintf(stream, ", ");
    print_fix16_t(stream, vector->z, 9, 4);
    fprintf(stream, ")");
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void print_qf16(FILE *stream, const qf16 *quat)
{
    print_fix16_t(stream, quat->a, 9, 4);
    fprintf(stream, " ");
    print_fix16_t(stream, quat->b, 9, 4);
    fprintf(stream, "i ");
    print_fix16_t(stream, quat->c, 9, 4);
    fprintf(stream, "j ");
    print_fix16_t(stream, quat->d, 9, 4);
    fprintf(stream, "k");
}
FIXMATHMATRIX_DEF FIXMATHMATRIX_INLINE void print_mf16(FILE *stream, const mf16 *matrix)
{
    if (matrix->errors)
    {
        fprintf(stream, "MATRIX ERRORS: %d\n", matrix->errors);
    }
    
    int row, column;
    for (row = 0; row < matrix->rows; row++)
    {
        for (column = 0; column < matrix->columns; column++)
        {
            fix16_t value = matrix->data[row][column];
            print_fix16_t(stream, value, 9, 4);
            fprintf(stream, " ");
        }
        fprintf(stream, "\n");
    }
}
