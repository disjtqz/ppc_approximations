#pragma once

/*

  This implementation is not verified to be bit-perfect. non-ieee mode has not been tested. Different exception settings has not been tested (they can produce different results?)
*/


#include <cstdint>
#if defined(_MSC_VER)
#include <intrin.h>
#endif
const unsigned char frsqrte_table2[] = {
    241u, 216u, 192u, 168u, 152u, 136u, 128u, 112u,
    96u,  76u,  60u,  48u,  32u,  24u,  16u,  8u };

static uint32_t lookup_frsqrte_table(unsigned v) {
    if (v < 16) {
        return frsqrte_table2[v];
    }
    else {
        return 0;
    }
}

static uint64_t frsqrte_clz64(uint64_t value) {
#if defined(_MSC_VER)
    return __lzcnt64(value);
#else
    return __builtin_clzll(value);
#endif
}
// checked against a real 360... for a single value (3.0)
double frsqrte(double x, bool non_ieee) {
    union {
        double d;
        uint64_t u;
    } result;
    uint64_t x_u;
reenter_for_denorm:
    x_u = *reinterpret_cast<uint64_t*>(&x);

    if (non_ieee) {

        if (x_u << 12 && !(x_u & 0x7FF0000000000000ULL)) {
            *reinterpret_cast<uint64_t*>(&x) &= (1ULL << 63);
            goto reenter_for_denorm;
        }
    }

    // is signed/unsigned zero?
    
    if (!(x_u << 1)) {
        result.u = (x_u & 0x8000000000000000ULL) | 0x7FF0000000000000ULL;
        return result.d;
    }

    if ((~x_u & 0x7FF0000000000000LL) == 0) {
        if (x_u == 0x7FF0000000000000ULL) {
            return .0;
        }
        if (x_u << 12 || x >= 0.0) {
            result.u = (~0x7FF0000000000000LL & x_u) | 0x7FF8000000000000ULL;
            return result.d;
        }

        result.u = 0x7FF8000000000000ULL;
        return result.d;
    }

    if (x < 0.0) {
        result.u = 0x7FF8000000000000;
        return result.d;
    }

    unsigned int biased_exponent = (x_u >> 52) & 0x7FFLL;
    uint64_t input_mantissa = x_u & 0xFFFFFFFFFFFFFLL;
    if (input_mantissa != 0 && !biased_exponent) {
        uint32_t mantissa_leading_zeros = frsqrte_clz64(input_mantissa);

        input_mantissa <<= mantissa_leading_zeros - 11;

        biased_exponent = ((x_u >> 52) & 0x700 | 0xC) - mantissa_leading_zeros;
    }
    uint32_t high_mantissa_index =
        ((input_mantissa >> 49) & 7 | (8 * (biased_exponent & 0xff)) & 8) ^ 8;

    uint32_t new_exponent = 1022 - ((biased_exponent - 1023) >> 1);

    uint64_t exponent52 = static_cast<uint64_t>(new_exponent) << 52;
    result.u =
        exponent52 |
        (static_cast<uint64_t>(lookup_frsqrte_table(high_mantissa_index)) << 44);

    if (non_ieee) {
        if (result.u << 12 && !(result.u & 0x7FF0000000000000LL)) {
            result.u &= (1ULL << 63);
        }
    }
    return result.d;
}
