#pragma once 
#include <cstdint> 

#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <x86intrin.h>
#endif
//only three possible increments that the approx advances by: 0, 0x20 and 0x40
//2 bits per element
//if bit 0 is set, advance by 0x20
//if bit 1 is set, advance by 0x40
//if neither set, advance by 0
//this allows us to calculate the value within a word using popcnt and multiply
static constexpr uint8_t BYTE_MASK_DIFF20 = 0b01010101;
static constexpr uint64_t PERMUTE8 = 0x0101010101010101ULL;
static constexpr uint64_t DIFF20_MASK = BYTE_MASK_DIFF20 * PERMUTE8;
static constexpr uint64_t DIFF40_MASK = ((~BYTE_MASK_DIFF20) & 0xff) * PERMUTE8;
static uint64_t get_diffset(uint32_t v) {
	if (v < 4) { return 0x1999999999999999ULL; }
	else if (v < 0xc) { return 0x1595959595959595ULL; }
	else if (v < 0x18) { return 0x1555955595559555ULL; }
	else if (v < 0x24) { return 0x1555555555555555ULL; }
	else if (v < 0x30) { return 0x1555155515551555ULL; }
	else { return 0x1515151515151515ULL; }
}
static uint32_t diffset_sum(uint64_t value) {
	return (__popcnt64(value & DIFF40_MASK) * 0x40) + (__popcnt64(value & DIFF20_MASK) * 0x20);
}
static constexpr uint16_t offset_table[] = {
	0 , 1536 , 3072 , 4608 , 6144 , 7424 , 8704 , 9984 , 11264 , 12544 , 13824 , 15104 , 16384 , 17536 , 18688 , 19840 , 20992 , 22144 , 23296 , 24448 , 25600 , 26752 , 27904 , 29056 , 30208 , 31232 , 32256 , 33280 , 34304 , 35328 , 36352 , 37376 , 38400 , 39424 , 40448 , 41472 , 42464 , 43360 , 44256 , 45152 , 46048 , 46944 , 47840 , 48736 , 49632 , 50528 , 51424 , 52320 , 53216 , 53984 , 54752 , 55520 , 56288 , 57056 , 57824 , 58592 , 59360 , 60128 , 60896 , 61664 , 62432 , 63200 , 63968 , 64736
};
static constexpr unsigned APPROX_TBL_BITS = 11;

float vlogefp(float f) {
  union {
    float _f;
    unsigned int _v;
  } bitu;
  bitu._f = f;
  uint32_t float_bits = bitu._v;

  if (float_bits >> 31) {
    // wronggg.

    // return .0f;
    bitu._v |= 0xFF800000U;
    return bitu._f;
  } else {
    uint32_t exp = float_bits >> 23;

    int32_t exp_unbiased = static_cast<int32_t>(exp) - 127;

    // todo: check out nan payloads

    if (!exp) {
      bitu._v = 0xFF800000U;
      return bitu._f;

    } else if (exp == 255) {
      bitu._v = 0x7FC00000;
      return bitu._f;
    } else {
      uint32_t mant = float_bits & ((1U << 23) - 1U);

      uint32_t approxbits = mant >> (23 - APPROX_TBL_BITS);
      uint32_t entry_word = (approxbits * 2) / 64;
      uint16_t offset = offset_table[entry_word];

      if (approxbits & 31) {
        uint32_t bit_begin = (approxbits * 2) & 63;

        uint64_t diffset = get_diffset(entry_word);

        offset += diffset_sum(
            diffset &
            ((1ULL << (bit_begin)) -
             1ULL));  // sum all preceeding entries in this 32 entry block
      }
      // i don't think double is actually necessary here.

      double fractional_part = static_cast<float>(static_cast<int>(
                                   static_cast<unsigned int>(offset))) /
                               65536.0;

      fractional_part += static_cast<double>(exp_unbiased);

      return static_cast<float>(fractional_part);
    }
  }
}
