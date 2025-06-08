/*===========================================================================*/
/* FIX16 FAST FOURIER TRANSFORM FUNCTIONS - GROUP H                         */
/*===========================================================================*/

/* Fast bit reversal for 32-bit numbers - optimized for ARM when available */
static uint32_t rbit_32(uint32_t x)
{
#if defined(__GNUC__) && defined(__ARM_ARCH_7M__)
    __asm__("rbit %0,%0" : "=r"(x) : "0"(x));
    return x;
#else
    x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
	x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
	x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
	x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
	return((x >> 16) | (x << 16));
#endif
}

/* Reverse bits in an n-bit number */
static uint32_t rbit_n(uint32_t x, unsigned n)
{
    return rbit_32(x << (32 - n));
}

/* Base-2 integer logarithm for FFT */
static int ilog2_fft(unsigned x)
{
    int result = -1;
    while (x)
    {
        x >>= 1;
        result++;
    }
    return result;
}

/* Fast calculation of DFT for a 4-point signal */
static void four_point_dft(INPUT_TYPE *input, unsigned input_stride,
                           fix16_t *real, fix16_t *imag)
{
    fix16_t x0 = INPUT_CONVERT(input[0 * input_stride]);
    fix16_t x1 = INPUT_CONVERT(input[1 * input_stride]);
    fix16_t x2 = INPUT_CONVERT(input[2 * input_stride]);
    fix16_t x3 = INPUT_CONVERT(input[3 * input_stride]);
    
    real[0] = x0 + x1 + x2 + x3;
    imag[0] = 0;
    real[1] = x0 - x2;
    imag[1] = -x1 + x3;
    real[2] = x0 - x1 + x2 - x3;
    imag[2] = 0;
    real[3] = x0 - x2;
    imag[3] = x1 - x3;
}

/* FFT butterfly operation - mix transforms pairwise */
static void butterfly(fix16_t *real, fix16_t *imag, unsigned blocksize, unsigned blockpairs)
{
    unsigned i, j;
    for (i = 0; i < blocksize; i++)
    {
        fix16_t angle = fix16_pi * i / blocksize;
        fix16_t c = fix16_cos(angle);
        fix16_t s = -fix16_sin(angle);
        
        fix16_t *rp = real + i;
        fix16_t *ip = imag + i;
        for (j = 0; j < blockpairs; j++)
        {
            /* Get the odd-indexed transform and multiply by sine */
            fix16_t re = fix16_mul(rp[blocksize], c) - fix16_mul(ip[blocksize], s);
            fix16_t im = fix16_mul(ip[blocksize], c) + fix16_mul(rp[blocksize], s);
            
            /* Update the transforms */
            rp[blocksize] = rp[0] - re;
            ip[blocksize] = ip[0] - im;
            rp[0] += re;
            ip[0] += im;
            
            rp += blocksize * 2;
            ip += blocksize * 2;
        }
    }
}

/* Real-input FFT implementation - main public API function */
void fix16_fft(INPUT_TYPE *input, fix16_t *real, fix16_t *imag, unsigned transform_length)
{
    int log_length = ilog2_fft(transform_length);
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