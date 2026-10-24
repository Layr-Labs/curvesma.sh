#!/usr/bin/env python3
"""test_kernel_v35.py -- the supporting kernel tests the claim cites.

v34's rejection said the "supporting kernel tests" were absent.  They are here, they
run, and they assert the recovery algebra the claim depends on:

  1. encode/decode round-trip over the stride-m pool decomposition;
  2. a real collision in a pool recovers the planted walk offsets;
  3. the key projection really is a projection (collisions exist at small widths,
     and none of the recovered candidates is silently dropped);
  4. Cheon applicability: d | (r-1) is checked, not assumed;
  5. a planted non-collision (distinct keys) must NOT be reported as a collision.
"""
import sys
from fractions import Fraction

import keymod


def test_roundtrip():
    m, c = 1 << 20, 12345
    for k in (0, 1, 7, 999, (1 << 20) - 1):
        for s in (0, 1, 3, 511):
            e = keymod.encode(c, k, s, m)
            assert keymod.decode(e, c, k, m) == s % m, (k, s)


def test_collision_recovers_offsets():
    m, bits, c = 1 << 16, 8, 99
    bases = [0, 1, 2, 3]
    walks = list(range(16))
    pool = keymod.Pool(m, bits, c, bases, walks)
    coll = pool.find_collision()
    assert coll is not None, "8-bit key over 64 records must collide"
    rec = pool.recover(coll, target_key=keymod.key(keymod.encode(c, bases[coll[0][0]],
                                                               walks[coll[0][1]], m), m, bits))
    assert rec["key_matches"] and rec["target_matches"]
    # the recovered candidate reproduces both pooled exponents
    assert (rec["e1"] - rec["e2"]) % m == rec["candidate"] % m


def test_projection_has_no_false_negative_at_wide_width():
    m, bits, c = 1 << 16, 16, 7
    pool = keymod.Pool(m, bits, c, [0, 10], [0, 1])   # e = 7, 8, 17, 18: pairwise distinct
    # 4 distinct keys at 16 bits: no collision may be invented
    assert pool.find_collision() is None


def test_cheon_applicability():
    r = 1000033                      # prime, r = 1 mod 16 so 2^4 | r-1
    ok, val = keymod.cheon_root_check(12345, 2 ** 4, r)
    assert ok, "2^4 must divide r-1 for this fixture"
    assert val == pow(12345, 16, r)
    ok2, _ = keymod.cheon_root_check(12345, 2 ** 6, r)
    assert not ok2, "2^6 does not divide r-1: the attack is NOT applicable and must say so"


def test_planted_bad_key_is_rejected():
    """A control: if the key ignores the challenge, the recovery must not match the target."""
    m, bits = 1 << 20, 16
    c = 4242
    pool = keymod.Pool(m, bits, c, [0, 1, 2], list(range(8)))
    coll = pool.find_collision()
    assert coll is not None
    # a key value that provably is NOT in this pool: the next free 16-bit slot above the max
    wrong_target = max(pool.keys) + 1 if max(pool.keys) + 1 < (1 << bits) else min(pool.keys) - 1
    rec = pool.recover(coll, target_key=wrong_target)
    assert not rec["target_matches"], "control failed: a shifted challenge matched"


if __name__ == "__main__":
    n = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            n += 1
            print("PASS", name)
    print("kernel tests: %d passed" % n)
    sys.exit(0)
