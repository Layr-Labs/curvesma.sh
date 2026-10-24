#!/usr/bin/env python3
"""v26 transfer arm (ratio layer): measure the nw=1 hitting probability at the
SAME design density N^2/m = 0.085 as the small-m curve arm, on the H subgroup of
the 268M-bit toy field (d=8, m=2^23.0).  24.4x the small arm's m.

Why the ratio layer is the right object: a solve requires v*b_j/a_i in mu6 (the
walk entries cancel at nw=1).  The curve arm's measured q (0.39907) equals that
ratio-hitting probability (law 0.39848); key truncation and affine verification
only filter candidates, they cannot create a solve.  So p_solve = p_ratio, and
the m-transfer question is a question about the ratio set.

Byproduct (verdict flaw 2): the script records the EXACT set cardinality |S| of
distinct ratios per trial, so the diagonal/ordered-pair counting question is
answered by a set count, not by a formula: p_pair uses N^2, p_set uses |S|.
"""
import json, math, os, random, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
R = 67122049
D = 8
M = (R - 1) // D
N = 844
assert 6 * N * N / M <= 0.51
TRIALS = int(sys.argv[1]) if len(sys.argv) > 1 else 4096
CHUNK = 256


def mu6_gen(r):
    # zeta has order 6 in Z_r^* (the automorphism scalars); it need not lie in H.
    g = 2
    while True:
        z = pow(g, (r - 1) // 6, r)
        if pow(z, 6, r) == 1 and pow(z, 3, r) != 1 and pow(z, 2, r) != 1:
            return z
        g += 1


def main():
    rng = random.SystemRandom()
    z = mu6_gen(R)
    zug = [pow(z, j, R) for j in range(6)]
    out = dict(kind="v26 TRANSFER ARM (ratio layer), fixed design density",
               p_field=R, d=D, m=M, m_log2=math.log2(M), nb=N, nw=1,
               K=6 * N * N / M, trials_planned=TRIALS, partial=0, hits=0,
               distinct_sum=0, distinct_min=None, distinct_max=None,
               p_pair_law=1 - math.exp(-6 * N * N / M),
               p_set_law_expected=None, seconds=0.0, done=False,
               scope="ratio-hitting probability; equals p_solve because every "
                     "true ratio hit is the same point and passes affine verification")
    t0 = time.time()
    for t in range(TRIALS):
        a = [pow(rng.randrange(2, R), D, R) for _ in range(N)]
        b = [pow(rng.randrange(2, R), D, R) for _ in range(N)]
        ib = [pow(x, R - 2, R) for x in b]
        S = set()
        add = S.add
        for x in a:
            for y in ib:
                add(x * y % R)
        v = pow(rng.randrange(2, R), D, R)
        hit = any(v * u % R in S for u in zug)
        out["hits"] += hit
        out["partial"] = t + 1
        out["distinct_sum"] += len(S)
        out["distinct_min"] = len(S) if out["distinct_min"] is None else min(out["distinct_min"], len(S))
        out["distinct_max"] = len(S) if out["distinct_max"] is None else max(out["distinct_max"], len(S))
        if (t + 1) % CHUNK == 0 or t + 1 == TRIALS:
            n = t + 1
            out["seconds"] = time.time() - t0
            out["q"] = out["hits"] / n
            out["distinct_mean"] = out["distinct_sum"] / n
            out["p_set_law_expected"] = 1 - math.exp(-6 * out["distinct_mean"] / M)
            import direct_v25 as dv
            lo = dv.cp95_lower(out["hits"], n)
            hi = dv.cp95_upper(out["hits"], n)
            out["cp95_lower"] = lo
            out["cp95_upper"] = hi
            out["prediction_small_arm"] = 0.39907073974609375
            out["prediction_in_interval"] = bool(lo <= 0.39907073974609375 <= hi)
            out["pair_law_in_interval"] = bool(lo <= out["p_pair_law"] <= hi)
            out["set_law_in_interval"] = bool(lo <= out["p_set_law_expected"] <= hi)
            out["done"] = (t + 1 == TRIALS)
            json.dump(out, open(HERE + "/transfer-v26.json", "w"), indent=1)
            print(json.dumps({k: out[k] for k in ("partial", "q", "cp95_lower",
                             "cp95_upper", "distinct_mean", "p_pair_law",
                             "p_set_law_expected", "prediction_in_interval",
                             "pair_law_in_interval", "set_law_in_interval",
                             "seconds")}, indent=1), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
