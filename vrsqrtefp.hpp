#pragma once
#include <cstdint> 
#include <cmath>
#include <cstddef>
#ifdef _MSC_VER
#include <intrin.h>
#endif

#define		VRSQRTE_SUBTABLE_LENGTH		16
static constexpr unsigned int g_vrsqrtefp_tables[] = {
    0x568B4FD,
    0x4F3AF97,
    0x48DAAA5,
    0x435A618,
    0x3E7A1E4,
    0x3A29DFE,
    0x3659A5C,
    0x32E96F8,
    0x2FC93CA,
    0x2D090CE,
    0x2A88DFE,
    0x2838B57,
    0x26188D4,
    0x2438673,
    0x2268431,
    0x20B820B,
    0x3D27FFA,
    0x3807C29,
    0x33878AA,
    0x2F97572,
    0x2C27279,
    0x2926FB7,
    0x2666D26,
    0x23F6AC0,
    0x21D6881,
    0x1FD6665,
    0x1E16468,
    0x1C76287,
    0x1AF60C1,
    0x1995F12,
    0x1855D79,
    0x1735BF4
};
#ifdef _MSC_VER
static uint32_t vrsqrtefp_lzcnt(uint32_t x) {
    return __lzcnt(x);
}
#else
static uint32_t vrsqrtefp_lzcnt(uint32_t x) {
    return __builtin_clz(x);
}
#endif
/*
Value	Result
-			Infinity QNaN
<0			QNaN
-0			-Infinity
+0			+Infinity
+Infinity	+0
NaN			QNaN
*/

float vrsqrtefp(float x, bool njm) {

restart_for_njm:

    unsigned int x_u = *reinterpret_cast<unsigned int*>(&x);

    unsigned int input_unbiased_exp = (unsigned char)(x_u >> 23);
    int biased_input_exp = input_unbiased_exp - 127;

    unsigned int input_mantissa = x_u & 0x7FFFFF;

    if ((x_u & 0x7F800000) != 0 || !input_mantissa) {

        if (input_unbiased_exp == 255) {
            if (x >= .0f) {
                if (!input_mantissa) {
                    return .0f;
                }
                x_u |= 0x400000;
                return *reinterpret_cast<float*>(&x_u);
            }
            if (x_u == 0xFF800000) {
                uint32_t tmpres = 0x7FC00000;
                return *reinterpret_cast<float*>(&tmpres);
            }
        }
    }
    else {
        //id recommend removing this and flushing to zero beforehand yourself
        if (njm) {
            *reinterpret_cast<unsigned*>(&x) &= 0x80000000;

            goto restart_for_njm;
        }
        uint32_t denormal_lz = vrsqrtefp_lzcnt(input_mantissa);
        input_unbiased_exp = 9 - denormal_lz;
        input_mantissa = (x_u << (denormal_lz - 8)) & 0x7FFFFE;
        biased_input_exp = 9 - denormal_lz - 127;

        if ((x_u & 0x7FFFFFFF) == 0x400000) {
            //0x80400000 is another bizarre special case
            if (x_u >> 31) {
                uint32_t tmpres = 0x7fc00000;
                return *reinterpret_cast<float*>(&tmpres);
            }
            else {
                uint32_t tmpres = 0x5F34FD00;
                return *reinterpret_cast<float*>(&tmpres);
            }
        }
    }
    if (x < .0f && biased_input_exp == 128) {
        if (!input_mantissa) {
            if ((int32_t)x_u < 0) {
                return -NAN;
            }
            else {
                return NAN;
            }
        }
        x_u |= 0x400000;
        return *reinterpret_cast<float*>(&x_u);
    }
    if (input_mantissa) {
        if (biased_input_exp == 128) {
            x_u |= 0x400000;
            return *reinterpret_cast<float*>(&x_u);
        }
    }
    else if (!input_unbiased_exp) {
        uint32_t infinity_unsigned = 0x7F800000U;
        infinity_unsigned |= (x_u & (1U << 31));

        return *reinterpret_cast<float*>(&infinity_unsigned);
    }
    if (x < 0.0f) {
        //qnan
        return -NAN;
    }
    /*
        The highorder bits Bh specify an interval within the range of values
        from 1.0 to 2.0, as illustrated in Figure 9.

        Because the
        exponent is processed separately, a limited range need only
        be considered.
        For each interval, the function is approximated by a straight line segment. The values along each
        segment are computed using the value y0 of that segment at
        its left-hand edge, the slope m of the segment, and the displacement along the horizontal axis. The displacement is
        equal to the low-order fraction bits Bl
         of the input B operand, just to the right of Bh. The Bh bits are used to access
        the table which contains the y0 values and the absolute values of m. Both of these estimate functions have negative
        slopes for all positive values of B; therefore, to obtain the
        value of the function along the line segment, the product of
        m and Bl
         must be subtracted from y0.
    */

    uint32_t Bh = input_mantissa >> 19;
    uint32_t Bl = ((input_mantissa >> 9) & 0x3FF);
    //select which table (odd or even) to use, Bh is index into selected table
    uint32_t table_index = (Bh | (VRSQRTE_SUBTABLE_LENGTH * (uint8_t)biased_input_exp) & VRSQRTE_SUBTABLE_LENGTH) ^ VRSQRTE_SUBTABLE_LENGTH;
    int32_t exponent_res = (int32_t)(127 - input_unbiased_exp) >> 1;
    //packed slope/offset
    uint32_t table_word = g_vrsqrtefp_tables[table_index];
    uint32_t slope = table_word >> 16;
    uint32_t offset = table_word & 0xFFFF;
    uint32_t gen_mantissa = (offset << 10) - Bl * slope;

    if ((gen_mantissa & 0x2000000) == 0) {
        //masking is probably unnecessary
        unsigned normalize_lz = vrsqrtefp_lzcnt(gen_mantissa & 0x1FFFFFF);
        gen_mantissa <<= normalize_lz - 6;
        exponent_res = exponent_res - normalize_lz + 6;
    }
    //round up
    if ((gen_mantissa & 5) != 0 && (gen_mantissa & 2) != 0) {
        gen_mantissa += 4;
    }
    float result;
    unsigned intermed_result = x_u & 0x80000000 | ((exponent_res << 23) + 0x3F800000) | (gen_mantissa >> 2) & 0x7FFFFF;
    
    if (njm) {
        if ((intermed_result & 0x7F800000) == 0 && (intermed_result & 0x7FFFFF)!=0) {
            intermed_result &= (1U << 31);
        }
    }
    *reinterpret_cast<unsigned int*>(&result) = intermed_result;
    
    return result;
}
