"""Recompute W and recover scalars with the production adding walk."""
import json
import sys
from decimal import Decimal, localcontext

TRACK = 'bls12-381-g1-ecdlp-plain'
ORDER = 52435875175126190479447740508185965837690552500527637822603658699938581184513
C = 2000
M_INS = 512
A_INS = 128
GM = 16
GA = 24
EXPECTED_W = 346538841911651323066977771571551249526064
EXPECTED_NUM = 137992069
EXPECTED_DEN = 1000000
EXPECTED_HITS = {'10007': 80, '100003': 79, '1000003': 80}
RUN_CURVE = True

"""Partition, jump coefficients, and Floyd recovery for the adding walk.

The same statements are what the proof specifies. The integer walk uses the
group element as its canonical word. The curve walk uses canonical_word on
an affine point.
"""

import math

MASK64 = (1 << 64) - 1
GOLDEN = 0x9E3779B97F4A7C15
MIX_A = 0xBF58476D1CE4E5B9
MIX_B = 0x94D049BB133111EB


def mix64(z):
    z = (z + GOLDEN) & MASK64
    z = ((z ^ (z >> 30)) * MIX_A) & MASK64
    z = ((z ^ (z >> 27)) * MIX_B) & MASK64
    return (z ^ (z >> 31)) & MASK64


def partition(word, seed):
    """Index in 0..31. Production, the proof, and this file use this body."""
    mixed = mix64(word & MASK64)
    k = (2 * (seed & MASK64) + 1) & MASK64
    return ((mixed * k) >> 32) & 31


def coefficient_block(seed, index):
    acc = 0
    seed &= MASK64
    for j in range(8):
        w = mix64((seed + ((index * 8 + j + 1) * GOLDEN)) & MASK64)
        acc = (acc << 64) | w
    return acc


def jump_coefficients(seed, order):
    table = []
    for i in range(32):
        u = coefficient_block(seed, 2 * i)
        v = coefficient_block(seed, 2 * i + 1)
        c = 1 + (u % (order - 1))
        d = v % order
        table.append((c, d))
    return table


def canonical_word(point):
    """Affine (x, y), or None for the identity. Low 64 bits only."""
    if point is None:
        return 0
    x, y = point
    return ((x & MASK64) ^ (((y & MASK64) << 1) & MASK64)) & MASK64


def trial_stream(master, samples, order):
    for i in range(samples):
        seed = mix64((master + (i + 1) * GOLDEN) & MASK64)
        secret = coefficient_block(master ^ 0xA5A5A5A5A5A5A5A5, i) % order
        yield seed, secret


def floyd_exponent(order, secret, seed):
    """Floyd on Z/nZ. Representative of k is k itself.

    Recovery is (alpha_h - alpha_t) * (beta_t - beta_h)^(-1).
    The negated numerator is computed only to show that it is -secret.
    """
    secret %= order
    if secret == 0:
        return {"ok": True, "recovered": 0, "negated": 0, "iters": 0, "short": True}
    table = jump_coefficients(seed, order)
    tp = ta = tb = 0
    hp = ha = hb = 0
    budget = 3 * _ceil_sqrt(order)
    for i in range(1, budget + 1):
        tp, ta, tb = _step_exp(tp, ta, tb, table, seed, secret, order)
        hp, ha, hb = _step_exp(hp, ha, hb, table, seed, secret, order)
        hp, ha, hb = _step_exp(hp, ha, hb, table, seed, secret, order)
        if tp == hp and (ta, tb) != (ha, hb):
            db = (tb - hb) % order
            if db == 0:
                continue
            da = (ha - ta) % order
            recovered = da * pow(db, -1, order) % order
            negated = (-da) * pow(db, -1, order) % order
            return {
                "ok": recovered == secret,
                "recovered": recovered,
                "negated": negated,
                "iters": i,
                "short": False,
                "neg_differs": negated != secret,
            }
    return {"ok": False, "recovered": None, "negated": None, "iters": budget, "short": False}


def _step_exp(point, a, b, table, seed, secret, order):
    c, d = table[partition(point, seed)]
    return (point + c + d * secret) % order, (a + c) % order, (b + d) % order


def _ceil_sqrt(n):
    s = math.isqrt(n)
    return s if s * s == n else s + 1


class Affine:
    """Short Weierstrass group in canonical affine coordinates."""

    def __init__(self, p, a):
        self.p = p
        self.a = a % p

    def add(self, p1, p2):
        if p1 is None:
            return p2
        if p2 is None:
            return p1
        p = self.p
        x1, y1 = p1
        x2, y2 = p2
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None
            lam = (3 * x1 * x1 + self.a) % p * pow((2 * y1) % p, p - 2, p) % p
        else:
            lam = (y2 - y1) % p * pow((x2 - x1) % p, p - 2, p) % p
        x3 = (lam * lam - x1 - x2) % p
        y3 = (lam * (x1 - x3) - y1) % p
        return (x3, y3)

    def mul(self, k, point):
        result = None
        base = point
        while k:
            if k & 1:
                result = self.add(result, base)
            base = self.add(base, base)
            k >>= 1
        return result


# Prime-order curve y^2 = x^3 + 2 over the field of 10099 elements.
# Group order 9967. Generator (3, 7278). Both facts are rechecked at runtime.
CURVE_P = 10099
CURVE_A = 0
CURVE_B = 2
CURVE_N = 9967
CURVE_G = (3, 7278)


def floyd_curve(secret, seed, curve=None):
    curve = curve or Affine(CURVE_P, CURVE_A)
    order = CURVE_N
    secret %= order
    base = CURVE_G
    target = curve.mul(secret, base)
    if secret == 0:
        return {
            "ok": target is None,
            "recovered": 0,
            "negated": 0,
            "iters": 0,
            "short": True,
        }
    table = jump_coefficients(seed, order)
    jumps = []
    for c, d in table:
        jumps.append(curve.add(curve.mul(c, base), curve.mul(d, target)))
    tortoise = hare = None
    ta = tb = ha = hb = 0
    budget = 3 * _ceil_sqrt(order)
    for i in range(1, budget + 1):
        tortoise, ta, tb = _step_curve(tortoise, ta, tb, jumps, table, seed, order, curve)
        hare, ha, hb = _step_curve(hare, ha, hb, jumps, table, seed, order, curve)
        hare, ha, hb = _step_curve(hare, ha, hb, jumps, table, seed, order, curve)
        if tortoise == hare and (ta, tb) != (ha, hb):
            db = (tb - hb) % order
            if db == 0:
                continue
            da = (ha - ta) % order
            recovered = da * pow(db, -1, order) % order
            negated = (-da) * pow(db, -1, order) % order
            check = curve.mul(recovered, base) == target
            return {
                "ok": recovered == secret and check,
                "recovered": recovered,
                "negated": negated,
                "iters": i,
                "short": False,
                "neg_differs": negated != secret,
            }
    return {"ok": False, "recovered": None, "negated": None, "iters": budget, "short": False}


def _step_curve(point, a, b, jumps, table, seed, order, curve):
    h = partition(canonical_word(point), seed)
    c, d = table[h]
    return curve.add(point, jumps[h]), (a + c) % order, (b + d) % order


def summarize_exponent(order, master, samples):
    hits = 0
    neg_bad = 0
    ratios = []
    examples = []
    root = math.sqrt(order)
    for seed, secret in trial_stream(master, samples, order):
        result = floyd_exponent(order, secret, seed)
        if not result["ok"]:
            continue
        hits += 1
        if result["short"]:
            continue
        if not result["neg_differs"]:
            neg_bad += 1
        ratios.append(result["iters"] / root)
        if len(examples) < 2:
            examples.append(
                {
                    "secret": result["recovered"],
                    "negated": result["negated"],
                    "iters": result["iters"],
                }
            )
    ratios.sort()
    if not ratios:
        median = p90 = None
    else:
        mid = len(ratios) // 2
        if len(ratios) % 2:
            median = ratios[mid]
        else:
            median = (ratios[mid - 1] + ratios[mid]) / 2
        p90 = ratios[min(len(ratios) - 1, int(0.90 * (len(ratios) - 1)))]
    return {
        "order": order,
        "samples": samples,
        "hits": hits,
        "neg_bad": neg_bad,
        "mean": sum(ratios) / len(ratios) if ratios else None,
        "median": median,
        "p90": p90,
        "max": ratios[-1] if ratios else None,
        "lengths": len(ratios),
    }


def summarize_curve(master, samples):
    hits = 0
    neg_bad = 0
    ratios = []
    examples = []
    root = math.sqrt(CURVE_N)
    curve = Affine(CURVE_P, CURVE_A)
    for seed, secret in trial_stream(master, samples, CURVE_N):
        result = floyd_curve(secret, seed, curve)
        if not result["ok"]:
            continue
        hits += 1
        if result["short"]:
            continue
        if not result["neg_differs"]:
            neg_bad += 1
        ratios.append(result["iters"] / root)
        if len(examples) < 2:
            examples.append(
                {
                    "secret": secret % CURVE_N,
                    "recovered": result["recovered"],
                    "negated": result["negated"],
                    "iters": result["iters"],
                }
            )
    ratios.sort()
    if not ratios:
        median = p90 = None
    else:
        mid = len(ratios) // 2
        median = ratios[mid] if len(ratios) % 2 else (ratios[mid - 1] + ratios[mid]) / 2
        p90 = ratios[min(len(ratios) - 1, int(0.90 * (len(ratios) - 1)))]
    return {
        "order": CURVE_N,
        "samples": samples,
        "hits": hits,
        "neg_bad": neg_bad,
        "mean": sum(ratios) / len(ratios) if ratios else None,
        "median": median,
        "p90": p90,
        "max": ratios[-1] if ratios else None,
        "examples": examples,
        "lengths": len(ratios),
    }


def ceil_sqrt(n):
    s = math.isqrt(n)
    return s if s * s == n else s + 1


def recompute():
    prime = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
    exponent = prime - 2
    fermat = (exponent.bit_length() - 1) + (exponent.bit_count() - 1)
    if fermat != 608:
        raise SystemExit("Fermat multiplication count changed")
    iterations = 3 * ceil_sqrt(ORDER)
    steps = 9 * ceil_sqrt(ORDER)
    inv = fermat * M_INS + fermat * 32
    step = inv + 4 * M_INS + 8 * A_INS
    walk_ins = steps * step
    group_ins = GM * M_INS + GA * A_INS
    jac_ins = (64 * 520 + 520) * group_ins
    norm_ins = 32 * inv + 32 * 6 * M_INS
    eq = iterations * (8 * M_INS + 8 * A_INS)
    coeff = iterations * 3 * 2 * 256
    hash_ins = iterations * 3 * 256
    derived = 2 * 64 * (8 * 40 + 16384)
    total = walk_ins + jac_ins + norm_ins + eq + coeff + hash_ins + derived + (1 << 40)
    return {"W": (total + C - 1) // C, "iterations": iterations, "derived": derived}


def log_ok(value, num, den):
    with localcontext() as ctx:
        ctx.prec = 80
        log2 = Decimal(value).ln() / Decimal(2).ln()
        claimed = Decimal(num) / Decimal(den)
    if claimed < log2:
        raise SystemExit("declared score is below log2(W)")
    if claimed - log2 > Decimal("0.000002"):
        raise SystemExit("declared score is more than two ticks above log2(W)")


def main():
    payload = json.load(sys.stdin)
    samples = int(payload["samples"])
    master = int(payload["master"])
    if not 40 <= samples <= 80:
        raise SystemExit("samples out of declared range")
    if not 0 <= master <= MASK64:
        raise SystemExit("master out of declared range")
    if partition(1, 1) != 5:
        raise SystemExit("partition vector changed")
    if partition(canonical_word((1, 1)), 1) != 14:
        raise SystemExit("canonical partition vector changed")
    if canonical_word(None) != 0:
        raise SystemExit("identity word changed")
    got = recompute()
    if got["W"] != EXPECTED_W or got["derived"] != 2138112:
        raise SystemExit("recomputed bound does not match the claim")
    log_ok(got["W"], EXPECTED_NUM, EXPECTED_DEN)
    zero = floyd_exponent(10007, 0, 1)
    if not (zero["ok"] and zero["short"] and zero["recovered"] == 0):
        raise SystemExit("identity secret was not recovered")
    reports = []
    for order in (10007, 100003, 1000003):
        summary = summarize_exponent(order, master, samples)
        if summary["neg_bad"]:
            raise SystemExit("opposite numerator matched a nonzero secret")
        if summary["hits"] * 4 < samples * 3:
            raise SystemExit("success rate is below 0.75")
        if master == 12648430 and samples == 80 and summary["hits"] != int(EXPECTED_HITS[str(order)]):
            raise SystemExit("fixed sample hit count changed")
        reports.append({
            "order": summary["order"],
            "hits": summary["hits"],
            "samples": summary["samples"],
            "mean": summary["mean"],
            "median": summary["median"],
            "p90": summary["p90"],
            "maximum": summary["max"],
            "neg_bad": summary["neg_bad"],
        })

    rhs = (CURVE_G[0] ** 3 + CURVE_B) % CURVE_P
    if pow(CURVE_G[1], 2, CURVE_P) != rhs:
        raise SystemExit("generator is not on the curve")
    group = Affine(CURVE_P, CURVE_A)
    if group.mul(CURVE_N, CURVE_G) is not None or group.mul(1, CURVE_G) != CURVE_G:
        raise SystemExit("curve order check failed")
    curve = summarize_curve(master, 40)
    if curve["neg_bad"] or curve["hits"] * 4 < 40 * 3:
        raise SystemExit("affine curve success rate is below 0.75")
    if master == 12648430 and curve["hits"] != 40:
        raise SystemExit("affine curve hit count changed")
    example = curve["examples"][0]
    if master == 12648430 and (example["secret"], example["recovered"], example["negated"]) != (2444, 2444, 7523):
        raise SystemExit("sign example failed")
    zero_c = floyd_curve(0, 1, group)
    if not (zero_c["ok"] and zero_c["short"] and zero_c["recovered"] == 0):
        raise SystemExit("curve identity was not recovered")
    curve_out = {
        "order": curve["order"],
        "hits": curve["hits"],
        "samples": curve["samples"],
        "mean": curve["mean"],
        "median": curve["median"],
        "p90": curve["p90"],
        "maximum": curve["max"],
        "example": example,
    }

    print(json.dumps({
        "track": TRACK,
        "W": str(got["W"]),
        "score": [EXPECTED_NUM, EXPECTED_DEN],
        "exponent_walk": reports,
        "affine_curve": curve_out,
    }))


if __name__ == "__main__":
    main()
