#!/usr/bin/env python3
"""Routing: one mixed 256-bit word per class key, with bit-disjoint fields.

Layout (least significant bit first), fixed at production width:
    bits  0..13   14  table index          (2^14 multipliers)
    bits 14..45   32  Frobenius twist      (reduced mod 12; 2^-32 selection bias)
    bits 48..127  80  distinguished point  (density 2^-80)
    bits 128..130  3  step type            (squaring if the low three bits are < 3)

The table index and the distinguished test therefore come from disjoint bits of
one independently mixed word; no field is a truncation of another and no field
shares a limb with another.  A 64-bit mixer would force the 80-bit distinguished
field and the table index into the same 64 bits, which collapses one of them; the
mixer below keeps 256 bits of state and is checked for collapse by selftest().
"""
import constants as K

MASK = (1 << 256) - 1
TABLE_BITS, TABLE_OFF = 14, 0
TWIST_BITS, TWIST_OFF = 32, 14
DP_BITS, DP_OFF = 80, 48
STEP_BITS, STEP_OFF = 3, 128

_M1 = 0x9E3779B97F4A7C15F39CC0605CEDC8341082276BF3A27251F86C6A11D0C18E7D
_M2 = 0xC2B2AE3D27D4EB4F165667B19E3779F9B4A7C15F39CC0605CEDC8341082276BF


def mix(key, bits=381):
    """Independent 256-bit mixing of the class key.

    Two multiplicative rounds over disjoint halves of the key, then a splitmix-style
    avalanche.  The key is taken as an integer in [0, 2^bits); nothing is truncated
    before mixing, so no two distinct keys collide by construction of the fields.
    """
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


def fields(key):
    h = mix(key)
    table = (h >> TABLE_OFF) & ((1 << TABLE_BITS) - 1)
    twist = ((h >> TWIST_OFF) & ((1 << TWIST_BITS) - 1)) % 12
    dp = (h >> DP_OFF) & ((1 << DP_BITS) - 1)
    step = (h >> STEP_OFF) & ((1 << STEP_BITS) - 1)
    return table, twist, dp, step


def disjoint():
    spans = [(TABLE_OFF, TABLE_BITS), (TWIST_OFF, TWIST_BITS),
             (DP_OFF, DP_BITS), (STEP_OFF, STEP_BITS)]
    for i in range(len(spans)):
        for j in range(i+1, len(spans)):
            a, la = spans[i]; b, lb = spans[j]
            if a < b + lb and b < a + la:
                return False
    return True


def selftest(n=120000, seed=20261023):
    """Histogram of the routing fields over n synthetic production-width keys."""
    import random
    rng = random.Random(seed)
    out = dict(keys=n, disjoint=disjoint(), table_distinct=0, table_chi2=None,
               dp_distinct=0, dp_zero=0, step_hist=[0]*8, twist_hist=[0]*12,
               table_min=0, table_max=0)
    seen_t = set(); seen_d = set(); hist = [0]*(1 << TABLE_BITS)
    for _ in range(n):
        key = rng.randrange(1, K.R_ORDER)
        t, tw, dp, st = fields(key)
        seen_t.add(t); hist[t] += 1
        seen_d.add(dp)
        if dp == 0:
            out['dp_zero'] += 1
        out['step_hist'][st] += 1
        out['twist_hist'][tw] += 1
    out['table_distinct'] = len(seen_t)
    out['dp_distinct'] = len(seen_d)
    out['table_min'] = min(hist)
    out['table_max'] = max(hist)
    # chi-square of the table index against uniform
    exp = n/(1 << TABLE_BITS)
    out['table_chi2'] = sum((c-exp)**2/exp for c in hist)
    ok = (out['disjoint'] and out['table_distinct'] > 0.95*(1 << TABLE_BITS)
          and out['dp_distinct'] == n and out['dp_zero'] == 0
          and min(out['step_hist']) > 0.9*n/8 and max(out['step_hist']) < 1.1*n/8
          and min(out['twist_hist']) > 0.85*n/12 and max(out['twist_hist']) < 1.15*n/12)
    out['collapse'] = 0 if ok else 1
    return out


if __name__ == "__main__":
    import json
    o = selftest()
    json.dump(o, open("routing-out.json", "w"), indent=1)
    print("keys %d  disjoint %s  table_distinct %d / %d  dp_distinct %d  dp_zero %d"
          % (o['keys'], o['disjoint'], o['table_distinct'], 1 << TABLE_BITS,
             o['dp_distinct'], o['dp_zero']))
    print("table chi2 %.1f (df %d)  step min/max %d/%d  twist min/max %d/%d"
          % (o['table_chi2'], (1 << TABLE_BITS)-1, min(o['step_hist']), max(o['step_hist']),
             min(o['twist_hist']), max(o['twist_hist'])))
    print("ROUTING", "PASS" if o['collapse'] == 0 else "FAIL", "collapse=%d" % o['collapse'])
