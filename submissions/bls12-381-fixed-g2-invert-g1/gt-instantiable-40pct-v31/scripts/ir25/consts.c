/* consts.c -- the extern constants kernels.c expects (v27 ships this file: v25/v26
   cite scripts/ir25/consts.c from verify25.py but did not include it, so the kernel
   tests and the price listing were not rebuildable from the package).  Values are the
   real BLS12-381 base field: P6 = p little-endian, NP0 = -p^{-1} mod 2^64,
   ONE6 = 2^384 mod p (Montgomery 1).  ZETA6 is a fixed nonzero residue used only by
   the f6 kernels; it is not on the priced path. */
#include <stdint.h>
const uint64_t P6[6] = {0xb9feffffffffaaabULL, 0x1eabfffeb153ffffULL, 0x6730d2a0f6b0f624ULL,
                        0x64774b84f38512bfULL, 0x4b1ba7b6434bacd7ULL, 0x1a0111ea397fe69aULL};
const uint64_t NP0 = 0x89f3fffcfffcfffdULL;
const uint64_t ZETA6[6] = {0x7609f4a70f9a1b2eULL, 0x1a3f5e2c8d4b6c11ULL, 0x4c9e0d33a7b8c2f5ULL,
                           0x0d2b7e61f4a3c985ULL, 0x62c8b40e5d1f7a93ULL, 0x0000000000000001ULL};
const uint64_t ONE6[6] = {0x760900000002fffdULL, 0xebf4000bc40c0002ULL, 0x5f48985753c758baULL,
                          0x77ce585370525745ULL, 0x5c071a97a256ec6dULL, 0x15f65ec3fa80e493ULL};
