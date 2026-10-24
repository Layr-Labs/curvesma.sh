#!/usr/bin/env python3
"""Routing: one mixed 256-bit word per CLASS KEY, with bit-disjoint fields.

Layout (least significant bit first), fixed at production width:
    bits  0..13    14  table index          (2^14 multipliers)
    bits 16..95    80  distinguished point  (density 2^-80)
    bits 112..114   3  step word            (squaring when the word is < 3)

v12c has ONE multiplier per table index.  v12b carried a 12-column table and a
32-bit column field, because its multiplier had to be re-derived through the
representative; here the multiplier is a function of the class key alone, so no
column is derived, no held offset exists and the field is DELETED rather than
declared-but-unused.  Storage drops from 2^14 x 12 to 2^14 elements.

The key is kappa(y) = c_0(y)^(1+p) in Fp, a 381-bit integer (constants.KAPPA_BITS).
Nothing is truncated before mixing.  Every branch input of the walk is a function
of the class key alone: NO window, NO checkpoint field, NO escape counter, NO
per-trail state beyond the labels (a,b).  A step is therefore a stationary map on
the classes, and two trails that meet at a class take the identical next step for
ever -- the merging lemma the solve's success model rests on.
"""
import constants as K

MASK = (1 << 256) - 1
TABLE_BITS, TABLE_OFF = 14, 0
DP_BITS, DP_OFF = 80, 16
STEP_BITS, STEP_OFF = 3, 112

_M1 = 0x9E3779B97F4A7C15F39CC0605CEDC8341082276BF3A27251F86C6A11D0C18E7D
_M2 = 0xC2B2AE3D27D4EB4F165667B19E3779F9B4A7C15F39CC0605CEDC8341082276BF


def mix(key, bits=None):
    """Independent 256-bit mixing of the class key.

    Two multiplicative rounds over disjoint halves of the key, then a
    splitmix-style avalanche.  The key is an integer in [0, 2^381).
    """
    if bits is None:
        bits = K.KAPPA_BITS
    key %= (1 << bits)
    lo = key & ((1 << 192) - 1)
    hi = key >> 192
    x = (lo * _M1 + 0x2545F4914F6CDD1D) & MASK
    y = (hi * _M2 + 0x9E3779B97F4A7C15) & MASK
    h = (x ^ ((y << 97) | (y >> 159))) & MASK
    h ^= h >> 33
    h = (h * 0xFF51AFD7ED558CCD) & MASK
    h ^= h >> 29
    h = (h * 0xC4CEB9FE1A85EC53) & MASK
    h ^= h >> 32
    return h


def fields(key, table_bits=None):
    """All three routing fields from ONE mix.  Production path: one mix per step.

    table_bits is a declared WIDTH parameter, not a different function: the toy
    runs the identical code at 6 bits because its whole table is 64 entries.
    """
    if table_bits is None:
        table_bits = K.TABLE_BITS
    h = mix(key)
    return ((h >> TABLE_OFF) & ((1 << table_bits) - 1),
            (h >> DP_OFF) & ((1 << DP_BITS) - 1),
            (h >> STEP_OFF) & ((1 << STEP_BITS) - 1))


def dp_hit(key, dp_bits=None):
    """The distinguished-point test, density 2^-dp_bits (2^-80 at production).

    A property of the class alone, so two trails at a common class agree on it
    and store the same record: the two reports describe one class reached by two
    label paths, which is exactly the solvable event.
    """
    _t, d, _st = fields(key)
    if dp_bits is None:
        return d == 0
    return (d & ((1 << dp_bits) - 1)) == 0


def disjoint():
    spans = [(TABLE_OFF, TABLE_BITS), (DP_OFF, DP_BITS), (STEP_OFF, STEP_BITS)]
    for i in range(len(spans)):
        for j in range(i + 1, len(spans)):
            a, la = spans[i]
            b, lb = spans[j]
            if a < b + lb and b < a + la:
                return False
    return True


def selftest(n=120000, seed=20261024):
    """Histogram of the routing fields over n synthetic production-width keys."""
    import random
    rng = random.Random(seed)
    out = dict(keys=n, disjoint=disjoint(), table_distinct=0, table_chi2=None,
               dp_distinct=0, dp_zero=0, step_hist=[0] * 8,
               table_min=0, table_max=0, dp80_fires=0, key_bits=K.KAPPA_BITS,
               has_column_field=False, has_ckpt_field=False, window=K.WINDOW)
    seen_t = set()
    seen_d = set()
    hist = [0] * (1 << TABLE_BITS)
    for _ in range(n):
        key = rng.randrange(1, 1 << K.KAPPA_BITS)
        t, dp, st = fields(key)
        out['dp80_fires'] += int(dp_hit(key))
        seen_t.add(t)
        hist[t] += 1
        seen_d.add(dp)
        if dp == 0:
            out['dp_zero'] += 1
        out['step_hist'][st] += 1
    out['table_distinct'] = len(seen_t)
    out['dp_distinct'] = len(seen_d)
    out['table_min'] = min(hist)
    out['table_max'] = max(hist)
    exp = n / (1 << TABLE_BITS)
    out['table_chi2'] = sum((c - exp) ** 2 / exp for c in hist)
    ok = (out['disjoint'] and out['table_distinct'] > 0.95 * (1 << TABLE_BITS)
          and out['dp_distinct'] == n and out['dp_zero'] == 0
          and min(out['step_hist']) > 0.9 * n / 8 and max(out['step_hist']) < 1.1 * n / 8)
    out['collapse'] = 0 if ok else 1
    return out


if __name__ == "__main__":
    import json
    o = selftest()
    json.dump(o, open("routing-out.json", "w"), indent=1)
    print("keys %d  disjoint %s  table_distinct %d / %d  dp_distinct %d  dp_zero %d"
          % (o['keys'], o['disjoint'], o['table_distinct'], 1 << TABLE_BITS,
             o['dp_distinct'], o['dp_zero']))
    print("table chi2 %.1f (df %d)  step min/max %d/%d"
          % (o['table_chi2'], (1 << TABLE_BITS) - 1, min(o['step_hist']),
             max(o['step_hist'])))
    print("ROUTING", "PASS" if o['collapse'] == 0 else "FAIL")
