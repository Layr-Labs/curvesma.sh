/* consts.c -- extern tables for the counted kernels (defeats constant folding). */
#include <stdint.h>
const uint64_t P6[6] = {0xffffffffffffaaabULL, 0xb9feffffffffffffULL, 0x1eabfffeb153ffffULL,
                        0x6730d2a0f6b0f624ULL, 0x64774b84f38512bfULL, 0x1a0111ea397fe69aULL};
/* -p^{-1} mod 2^64 */
const uint64_t NP0 = 0x89f3fffcfffcfffdULL;
/* zeta = gamma, a fixed base-field constant (canonical, non-Montgomery for counting only) */
const uint64_t ZETA6[6] = {0x3d3d4b3b1b1b1b1bULL, 0x9a0d8f6f5e4c3b2aULL, 0x7f6e5d4c3b2a1908ULL,
                           0x1f2e3d4c5b6a7988ULL, 0x1122334455667788ULL, 0x0a0b0c0d0e0f1011ULL};
