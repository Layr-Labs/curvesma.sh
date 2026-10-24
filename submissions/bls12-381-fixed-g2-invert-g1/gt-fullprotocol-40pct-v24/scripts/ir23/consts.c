/* consts.c -- extern tables for the counted kernels (defeats constant folding).
 *
 * v23 FIX 3 (harness): the v22 table transcribed the BLS12-381 modulus shifted by
 * one limb (limb0 0xffffffffffffaaab, limb1 0xb9feffffffffffff, limb4 0x64774b84...),
 * so no kernel could match a bigint reference through this harness.  The limbs below
 * are the canonical little-endian limbs of p (BLS12-381), and NP0 = -p^{-1} mod 2^64,
 * both re-derived and asserted in verify23.py.
 */
#include <stdint.h>
const uint64_t P6[6] = {0xb9feffffffffaaabULL, 0x1eabfffeb153ffffULL, 0x6730d2a0f6b0f624ULL,
                        0x64774b84f38512bfULL, 0x4b1ba7b6434bacd7ULL, 0x1a0111ea397fe69aULL};
const uint64_t NP0 = 0x89f3fffcfffcfffdULL;
const uint64_t ZETA6[6] = {0x3d3d4b3b1b1b1b1bULL, 0x9a0d8f6f5e4c3b2aULL, 0x7f6e5d4c3b2a1908ULL,
                           0x1f2e3d4c5b6a7988ULL, 0x1122334455667788ULL, 0x0a0b0c0d0e0f1011ULL};
/* Montgomery representation of 1: R mod p, R = 2^384 (the identity element of the
 * Montgomery domain -- v22 wrote the canonical constant 1 here, which is NOT the
 * identity in that domain). */
const uint64_t ONE6[6] = {0x760900000002fffdULL, 0xebf4000bc40c0002ULL,
                          0x5f48985753c758baULL, 0x77ce585370525745ULL,
                          0x5c071a97a256ec6dULL, 0x15f65ec3fa80e493ULL};
