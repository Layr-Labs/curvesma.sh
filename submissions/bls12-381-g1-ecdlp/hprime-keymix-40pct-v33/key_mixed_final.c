/* v32 KEY KERNEL: genuine 168-bit key over the FULL 381-bit residue (6 x 64-bit limbs).
 *
 * v31 fault (verdict 3e66901 flaw 1): mix64() masked its input to 64 bits, so the
 * "168-bit mixed digest" WAS a 64-bit key: residues differing by 2^64 collided and the
 * collision ledger was priced on a key space that did not exist.
 *
 * v32: three 64-bit lanes, each a Horner fold of ALL SIX residue limbs with its own
 * odd multiplier (no input mask anywhere), then a 2-round xorshift-multiply finaliser
 * per lane.  The stored 168-bit key is 56 bits from each finalised lane at disjoint
 * windows, so every stored bit is a function of all six limbs (avalanche MEASURED in
 * key_collision_test.py, not asserted).
 *
 * ONE call per record computes all three lanes and shares the limb loads: the
 * instruction count of this function's compiled body is the charged per-record cost
 * (otool -tV; the fleet's counted-kernel method).
 */
#include <stdint.h>

#define RES_LIMBS 6

static const uint64_t KL_MUL[3] = {0x9E3779B97F4A7C15ULL,
                                   0xC2B2AE3D27D4EB4FULL,
                                   0xD6E8FEB86659FD93ULL};
static const uint64_t KL_FIN[3] = {0xFF51AFD7ED558CCDULL,
                                   0xC4CEB9FE1A85EC53ULL,
                                   0x94D049BB133111EBULL};

/* lane = ((((l5*M + l4)*M + l3)*M + l2)*M + l1)*M + l0  -- all six limbs, one lane */
static inline uint64_t kl_lane(const uint64_t *limb, int lane) {
    uint64_t m = KL_MUL[lane], a = limb[5];
    for (int i = 4; i >= 0; i--) a = a * m + limb[i];
    a ^= a >> 33; a *= KL_FIN[lane];
    a ^= a >> 29; a *= KL_FIN[(lane + 1) % 3];
    a ^= a >> 31;
    return a;
}

/* 168 stored bits: bits [8,64) of lane0, [8,64) of lane1, [0,56) of lane2 -> out[3] */
void key_mixed_final(const uint64_t *limb, uint64_t out[3]) {
    out[0] = (kl_lane(limb, 0) >> 8) & 0xFFFFFFFFFFFFFFULL;
    out[1] = (kl_lane(limb, 1) >> 8) & 0xFFFFFFFFFFFFFFULL;
    out[2] = kl_lane(limb, 2) & 0xFFFFFFFFFFFFFFULL;
}
