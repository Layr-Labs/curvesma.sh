#!/usr/bin/env python3
"""v32 SINGLE SOURCE OF TRUTH for the key rule.

verdict 3e66901 flaw 1 (FATAL, again): v31's key_mixed() called mix64(), whose FIRST
statement masks the input to 64 bits, so the "168-bit mixed digest" was a 64-bit key:
residues differing by 2^64 got identical keys.  It also left record_of() calling
key_scalar(), i.e. the charged key was never the key in the hot path.

v32: key_mixed_final() folds ALL SIX 64-bit limbs of the 381-bit residue into each of
three 64-bit lanes with distinct odd multipliers (no input mask anywhere), finalises
each lane with a 2-round xorshift-multiply, and returns 168 bits taken from disjoint
windows of the three finalised lanes.  record_of() calls THIS function.

mix64()/key_mixed() are kept ONLY as the planted negative control: the collision test
asserts that they DO collide on the adversarial pair class and that v32 does not.
"""
import math

# ---- production group (BLS12-381 G1, y^2 = x^3 + 4)
P = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001

KEY_BITS = 168
KEY_IMAGE_BITS = 168          # v32 image is FULL width (v31's even-rep low bit forced 2^167)
ROUTE_BITS = 84
ROUTE_HOPS = 4
HOP_BITS = 21
SORT_BITS = 63
SORT_PASSES = 3
SORT_DIGIT = 21
SCAN_BITS = 21
KEY_LIMBS = 3
LIMB_BITS = 56
RES_LIMBS = 6                 # ceil(381/64): the FULL production residue
MASK64 = (1 << 64) - 1

RECORD_WORDS = 6
RECORD_BYTES = 48
INDEX_BITS = 128
INDEX_BOUND_LOG2 = 116
WALK_LOG2 = 50
HOME_BITS = 12

# constants: identical to key_mixed_final.c (preflight cross-checks the C literals)
KL_MUL = (0x9E3779B97F4A7C15, 0xC2B2AE3D27D4EB4F, 0xD6E8FEB86659FD93)
KL_FIN = (0xFF51AFD7ED558CCD, 0xC4CEB9FE1A85EC53, 0x94D049BB133111EB)


def residue_limbs(residue, limbs=RES_LIMBS):
    """Split the FULL residue into `limbs` 64-bit limbs.  No truncation: the caller may
    pass any width up to 384 bits and every bit lands in limb[0..5]."""
    assert residue >= 0
    return [(residue >> (64 * i)) & MASK64 for i in range(limbs)]


def _lane(limb, lane):
    m, a = KL_MUL[lane], limb[RES_LIMBS - 1]
    for i in range(RES_LIMBS - 2, -1, -1):
        a = (a * m + limb[i]) & MASK64          # Horner: ALL limbs, this lane's constant
    a ^= a >> 33
    a = (a * KL_FIN[lane]) & MASK64
    a ^= a >> 29
    a = (a * KL_FIN[(lane + 1) % 3]) & MASK64
    a ^= a >> 31
    return a


def key_mixed_final(residue):
    """168-bit key = 56 bits from each finalised lane at disjoint windows.
    Every output bit is a function of all six residue limbs (avalanche measured)."""
    l = residue_limbs(residue)
    k0 = (_lane(l, 0) >> 8) & ((1 << 56) - 1)
    k1 = (_lane(l, 1) >> 8) & ((1 << 56) - 1)
    k2 = _lane(l, 2) & ((1 << 56) - 1)
    return (k2 << 112) | (k1 << 56) | k0


def key_limbs(k, limb_bits=LIMB_BITS, limbs=KEY_LIMBS):
    return [(k >> (limb_bits * i)) & ((1 << limb_bits) - 1) for i in range(limbs)]


def even_member(y, p=P):
    return y if (y & 1) == 0 else p - y


def key_scalar(y, p=P, bits=KEY_BITS):
    """v30 key.  NOT used by record_of in v32; kept for the paired diagnostic only."""
    return even_member(y, p) & ((1 << bits) - 1)


# ---- planted negative control (v31's defect, reproduced on purpose) ----
MIX_M1 = 0xFF51AFD7ED558CCD
MIX_M2 = 0xC4CEB9FE1A85EC53


def mix64(x):
    x &= MASK64                       # <-- THE BUG: input masked to 64 bits
    x ^= x >> 33
    x = (x * MIX_M1) & MASK64
    x ^= x >> 33
    x = (x * MIX_M2) & MASK64
    x ^= x >> 33
    return x


def key_mixed_v31(residue, bits=KEY_BITS):
    return mix64(residue) & ((1 << bits) - 1)


# ---- production record: the hot path calls key_mixed_final, not key_scalar ----
HOT_PATH_CALLS = {"key_mixed_final": 0, "key_scalar": 0, "key_mixed_v31": 0}


def record_of(y, side, index, p=P, count=True):
    assert 0 <= index < (1 << INDEX_BITS)
    assert side in (0, 1)
    k = key_mixed_final(even_member(y, p))
    if count:
        HOT_PATH_CALLS["key_mixed_final"] += 1
    limbs = key_limbs(k)
    return [limbs[0], limbs[1], limbs[2], side,
            (index >> 64) & MASK64, index & MASK64]


def unpack(rec):
    k = rec[0] | (rec[1] << LIMB_BITS) | (rec[2] << (2 * LIMB_BITS))
    return k, rec[3], (rec[4] << 64) | rec[5]


def pack_index(a_pos, w_pos, walk_log2=WALK_LOG2):
    return (a_pos << walk_log2) | w_pos


LAYOUT_BITS = 3 * LIMB_BITS + 1 + INDEX_BITS     # 297 payload bits in 6 words


def layout_ok():
    assert KEY_BITS == ROUTE_BITS + SORT_BITS + SCAN_BITS == KEY_LIMBS * LIMB_BITS
    assert ROUTE_BITS == ROUTE_HOPS * HOP_BITS and SORT_BITS == SORT_PASSES * SORT_DIGIT
    assert LAYOUT_BITS == 297 and RECORD_WORDS * 64 == 384 >= LAYOUT_BITS
    assert KEY_IMAGE_BITS == KEY_BITS
    return True


# ---- charged per-record key cost: measured from the compiled C body ----
KEY_KERNEL_INSTRUCTIONS = 70      # otool -tV key_mixed_final: 70 instructions (all 3 lanes)
V31_KEY_ITEMS = 7                 # v31 charged mixed_digest_finalize 3 + key_pack 4
KEY_NET_DELTA = KEY_KERNEL_INSTRUCTIONS - V31_KEY_ITEMS     # +63 instructions/record

# self-test vectors: the C kernel must reproduce these (preflight checks a build)
SELFTEST = [
    (0, 0),
    (1, key_mixed_final(1)),
    (P - 1, key_mixed_final(P - 1)),
    (1 << 64, key_mixed_final(1 << 64)),
]

if __name__ == "__main__":
    print("layout_ok", layout_ok(), "KEY_BITS", KEY_BITS,
          "kernel", KEY_KERNEL_INSTRUCTIONS, "net_delta", KEY_NET_DELTA)
    for r in SELFTEST:
        print(r[0], hex(r[1]))
