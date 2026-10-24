#!/usr/bin/env python3
"""test_kernel_v36.py -- kernel tests that assert the PRODUCTION predicates, not offset
identities (the v35 defect the verdict named: "the kernel tests use low exponent bits
instead of EC coordinates and test offset identities, not recovery verified against T
and A_d").

Every test below runs on a real curve (ec.Curve), uses the real key x([e]P0) mod 2^b,
and the recovery tests check [s]P0 == T and [s^d]P0 == A_d.
"""
import random, sys

import ec, keymod

cv, G, r = ec.Curve(ec.P37), ec.G37, ec.R37


def test_lift_recovers_full_secret():
    rng = random.Random(11)
    for D in (8, 10):
        m = (r - 1) >> D
        for _ in range(20):
            s = rng.randrange(1, r)
            T = cv.mul(s, G)
            ad = cv.mul(pow(s, 1 << D, r), G)
            o = keymod.lift(cv, G, T, s % m, m, D, ad=ad, d=1 << D, r=r)
            assert o["found"] and o["s"] == s, (D, s)
            assert o["t_check"] and o["ad_check"]
            assert len(o["passes"]) == 1, "the lift must be unique"


def test_cartesian_product_is_the_whole_grid():
    m, bits, c = (r - 1) >> 8, 12, 12345
    nb, ell = 6, 9
    bases = [1, 2, 3, 5, 8, 13]
    walks = list(range(ell))
    pool = keymod.Pool(m, bits, c, bases, walks, keyfn=keymod.KeyFn(cv, G, bits, r))
    assert len(pool.records) == nb * ell, (len(pool.records), nb * ell)
    # every base really walks the SAME list: record (i,l) == [c+k_i+s_l]P0 for all i,l
    for i in range(nb):
        for l in range(ell):
            e = c + bases[i] + walks[l]
            assert pool.records[(i, l)] == cv.mul(e, G), (i, l)


def test_target_binding_relation():
    """a key match binds c, the offsets and s: s = +-(c+k_i+s_l) mod r, sign decided by
    the full point check.  This is the relation v35 did not have."""
    rng = random.Random(12)
    m, bits = (r - 1) >> 10, 12
    for _ in range(10):
        c = rng.randrange(1, m // 4)
        bases = rng.sample(range(0, m // 2), 20)
        walks = list(range(30))
        i, l = rng.randrange(20), rng.randrange(30)
        s = c + bases[i] + walks[l]
        T = cv.mul(s, G)
        pool = keymod.Pool(m, bits, c, bases, walks, keyfn=keymod.KeyFn(cv, G, bits, r))
        hits, _ = pool.hits_for_target(T)
        assert hits, "the planted target must be found"
        ok = any(cv.mul(cand, G) == T for (ii, ll) in hits
                 for cand in ((c + bases[ii] + walks[ll]) % r, (r - c - bases[ii] - walks[ll]) % r))
        assert ok


def test_truncated_key_alone_never_verifies():
    """a truncated coincidence must never be accepted as a secret"""
    rng = random.Random(13)
    m, bits = (r - 1) >> 10, 12
    c = rng.randrange(1, m // 4)
    bases = rng.sample(range(0, m // 2), 30)
    walks = list(range(40))
    pool = keymod.Pool(m, bits, c, bases, walks, keyfn=keymod.KeyFn(cv, G, bits, r))
    su = rng.randrange(1, r)
    Tu = cv.mul(su, G)
    hits, _ = pool.hits_for_target(Tu)
    for (ii, ll) in hits:
        e = c + bases[ii] + walks[ll]
        assert cv.mul(e % r, G) != Tu, "a truncated hit must not equal the target"
        o = keymod.lift(cv, G, Tu, e % m, m, 10)
        assert not (o["found"] and o.get("s") == su)


def test_corrupt_seed_fails():
    rng = random.Random(14)
    m = (r - 1) >> 10
    for _ in range(5):
        s = rng.randrange(1, r)
        T = cv.mul(s, G)
        bad = (s % m) ^ (1 << rng.randrange(0, m.bit_length()))
        assert not keymod.lift(cv, G, T, bad, m, 10)["found"]


if __name__ == "__main__":
    n = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            n += 1
            print("PASS", name)
    print("kernel tests v36: %d passed" % n)
    sys.exit(0)
