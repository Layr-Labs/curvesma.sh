#!/usr/bin/env python3
"""SINGLE SOURCE OF TRUTH for the key rule, the production record layout and the
position-index reconstruction.

Imported by toy.py (executes it), build_certificate.py (prices it) and
verify_certificate.py (re-derives it).  There is exactly ONE definition of the
key: the EVEN member of {y, p-y}, truncated to KEY_BITS bits.  p is odd, so
exactly one of y, p-y is even and the representative is unique.  Every other
file in the package calls these functions; no file re-implements the rule.
"""
import math

# ---- production group (BLS12-381 G1, y^2 = x^3 + 4)
P = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
GX = 0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb
GY = 0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1

# ---- key geometry.  KEY_BITS = ROUTE + SORT + SCAN = 84 + 63 + 21.
KEY_BITS = 168
KEY_IMAGE_BITS = 167          # the even representative's low bit is ALWAYS 0
ROUTE_BITS = 84
ROUTE_HOPS = 4
HOP_BITS = 21                 # 4 x 21 = 84 routed bits
SORT_BITS = 63
SORT_PASSES = 3
SORT_DIGIT = 21               # 3 x 21 = 63 sorted bits
SCAN_BITS = 21                # the remaining low 21 bits, explicit compare
KEY_LIMBS = 3
LIMB_BITS = 56                # 3 x 56 = 168

# ---- production record: 6 words = 48 bytes
RECORD_WORDS = 6
RECORD_BYTES = 48
WORD_BYTES = 8
INDEX_BITS = 128              # words 4..5
INDEX_BOUND_LOG2 = 116        # declared bound on a list position

# ---- deterministic list geometry
WALK_LOG2 = 50                # |W| = 2^50 shared walk-list elements
A_POOL_LOG2 = 66              # |A| per side (2^66 >= 2^65.65 used)
LIST_LOG2 = 116               # N < 2^116, so a position fits INDEX_BITS
HOME_BITS = 12                # gap-checkpoint spacing: 2^12 gaps per checkpoint


def even_member(y, p=P):
    """The unique EVEN member of {y, p-y}.  p odd: exactly one is even."""
    return y if (y & 1) == 0 else p - y


def key_scalar(y, p=P, bits=KEY_BITS):
    """Declared key: even representative of y, truncated to `bits` bits."""
    return even_member(y, p) & ((1 << bits) - 1)


def key_limbs(k, limb_bits=LIMB_BITS, limbs=KEY_LIMBS):
    return [(k >> (limb_bits * i)) & ((1 << limb_bits) - 1) for i in range(limbs)]


def route_of(k):
    return k >> (KEY_BITS - ROUTE_BITS)


def local_of(k):
    return k & ((1 << (KEY_BITS - ROUTE_BITS)) - 1)


def record_of(y, side, index, p=P):
    """Pack the 6 production words: 3 key limbs, side, 128-bit position index."""
    assert 0 <= index < (1 << INDEX_BITS)
    assert side in (0, 1)
    k = key_scalar(y, p)
    limbs = key_limbs(k)
    return [limbs[0], limbs[1], limbs[2], side,
            (index >> 64) & ((1 << 64) - 1), index & ((1 << 64) - 1)]


def unpack(rec):
    k = rec[0] | (rec[1] << LIMB_BITS) | (rec[2] << (2 * LIMB_BITS))
    side = rec[3]
    index = (rec[4] << 64) | rec[5]
    return k, side, index


def pack_index(a_pos, w_pos, walk_log2=WALK_LOG2):
    """Position = (pool slot, walk offset); the list element is a deterministic
    function of the position, so the 174-bit scalar is NOT stored."""
    assert 0 <= a_pos < (1 << (INDEX_BITS - walk_log2))
    assert 0 <= w_pos < (1 << walk_log2)
    return (a_pos << walk_log2) | w_pos


def split_index(index, walk_log2=WALK_LOG2):
    return index >> walk_log2, index & ((1 << walk_log2) - 1)


# ---- reconstruction cost, itemised (AArch64 instructions per candidate)
RECON_ITEMS = [
    ("index_load_and_mask", 3, "2 LDP of the two index words, AND/shift"),
    ("pool_load_2words", 4, "2 LDP from the stored 97-bit pool entry, addressing 2"),
    ("gap_index_shift", 2, "shift to checkpoint byte offset"),
    ("checkpoint_load_2words", 4, "2 LDP of the prefix sum at the checkpoint"),
    ("local_gap_walk", 6 * (1 << HOME_BITS),
     "2^12 gap adds: 2 loads, add, carry, store, loop = 6 each"),
    ("scalar_product_3x3", 9 * 383 + 24,
     "9 fp_mul (383 Ir each, COST_MODEL.md) + 24 accumulate/addressing for a*w"),
    ("ratio_and_six_test", 6 * (2 * 383 + 30),
     "6 candidate multipliers: one fp_mul + compare/select each, 796 per u"),
]
RECON_ITEMISED = sum(v for _, v, _ in RECON_ITEMS)
RECON_STRESSED = -(-RECON_ITEMISED * 6 // 5)      # family 6/5 stress rule


def recon_per_candidate():
    return RECON_STRESSED


def layout_rows():
    """Machine-generated record layout: word, byte offset, field, bits, range."""
    rows = []
    for i in range(KEY_LIMBS):
        rows.append(dict(word=i, byte=i * WORD_BYTES, field="key_limb%d" % i,
                         bits=LIMB_BITS, span="[%d,%d)" % (LIMB_BITS * i, LIMB_BITS * (i + 1)),
                         note="even-rep y bits; limb 0 carries the always-zero low bit"))
    rows.append(dict(word=3, byte=24, field="side", bits=1, span="bit 0",
                     note="0 = baby/X side, 1 = giant/Y side; 63 bits reserved"))
    rows.append(dict(word=4, byte=32, field="index_hi", bits=64, span="[64,128)",
                     note="position index high word"))
    rows.append(dict(word=5, byte=40, field="index_lo", bits=64, span="[0,64)",
                     note="position index low word; position < 2^116"))
    return rows


LAYOUT_BITS = 3 * LIMB_BITS + 1 + INDEX_BITS   # 297 payload bits; 6 words = 384 bits = 48 B


def layout_ok():
    rows = layout_rows()
    assert len(rows) == RECORD_WORDS
    assert max(r["byte"] for r in rows) + WORD_BYTES == RECORD_BYTES
    assert KEY_BITS == ROUTE_BITS + SORT_BITS + SCAN_BITS
    assert ROUTE_BITS == ROUTE_HOPS * HOP_BITS
    assert SORT_BITS == SORT_PASSES * SORT_DIGIT
    assert KEY_BITS == KEY_LIMBS * LIMB_BITS
    assert INDEX_BOUND_LOG2 < INDEX_BITS
    assert LAYOUT_BITS == 297 and RECORD_WORDS * 64 == 384 >= LAYOUT_BITS
    assert WALK_LOG2 + (INDEX_BITS - WALK_LOG2) == INDEX_BITS
    return True


if __name__ == "__main__":
    print("layout_ok", layout_ok(), "recon itemised", RECON_ITEMISED,
          "stressed", RECON_STRESSED, "log2 %.4f" % math.log2(RECON_STRESSED))
    for r in layout_rows():
        print(r)
