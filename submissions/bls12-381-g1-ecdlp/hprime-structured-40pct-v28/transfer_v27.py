#!/usr/bin/env python3
"""v27 TRANSFER ARM, big-n replication of the v26 ratio-layer arm.

Exact same design as transfer_v26.py (R=67122049, D=8, m=(R-1)/8=2^23.0003,
N=844, nw=1, design density N^2/m = 0.08507), executed over ALL physical cores
with independent SystemRandom streams per worker, so the pooled hit count is a
genuine n-trial binomial sample of the same event.

Why bigger n: v26 ran 2048 trials -> q=0.403809, CP-95 [0.385841, 0.421980].
The lower endpoint is below the 0.39 floor, so the BAND, not the point estimate,
failed the floor test. Required n for q~0.404 to clear 0.39:
  n >= (1.645^2 * 0.404 * 0.596) / 0.014^2 ~= 4600.
We run 20000 (4.3x the minimum) so the lower endpoint clears with margin even
if the true q sits at 0.399.
"""
import json, math, multiprocessing as mp, os, random, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
R = 67122049
D = 8
M = (R - 1) // D
N = 844
WORKERS = 24
TOTAL = int(sys.argv[1]) if len(sys.argv) > 1 else 20000


def mu6_gen(r):
    g = 2
    while True:
        z = pow(g, (r - 1) // 6, r)
        if pow(z, 6, r) == 1 and pow(z, 3, r) != 1 and pow(z, 2, r) != 1:
            return z
        g += 1


def chunk(args):
    wid, trials = args
    rng = random.SystemRandom()
    z = mu6_gen(R)
    zug = [pow(z, j, R) for j in range(6)]
    hits = 0
    dsum = 0
    dmin = None
    dmax = None
    t0 = time.time()
    for _ in range(trials):
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
        hits += hit
        c = len(S)
        dsum += c
        dmin = c if dmin is None else min(dmin, c)
        dmax = c if dmax is None else max(dmax, c)
    return dict(worker=wid, trials=trials, hits=hits, distinct_sum=dsum,
                distinct_min=dmin, distinct_max=dmax, seconds=time.time() - t0)


def main():
    base = TOTAL // WORKERS
    jobs = [(i, base + (1 if i < TOTAL - base * WORKERS else 0))
            for i in range(WORKERS)]
    t0 = time.time()
    with mp.Pool(WORKERS) as pool:
        parts = pool.map(chunk, jobs)
    n = sum(p["trials"] for p in parts)
    hits = sum(p["hits"] for p in parts)
    import direct_v25 as dv
    q = hits / n
    out = dict(kind="v27 TRANSFER ARM big-n (ratio layer, fixed design density)",
               p_field=R, d=D, m=M, m_log2=math.log2(M), nb=N, nw=1,
               K=6 * N * N / M, trials_planned=TOTAL, partial=n, hits=hits,
               q=q, cp95_lower=dv.cp95_lower(hits, n),
               cp95_upper=dv.cp95_upper(hits, n),
               distinct_mean=sum(p["distinct_sum"] for p in parts) / n,
               distinct_min=min(p["distinct_min"] for p in parts),
               distinct_max=max(p["distinct_max"] for p in parts),
               p_pair_law=1 - math.exp(-6 * N * N / M),
               workers=WORKERS, core_seconds=sum(p["seconds"] for p in parts),
               wall_seconds=time.time() - t0, per_worker=parts, done=True,
               floor=0.39,
               floor_met=bool(dv.cp95_lower(hits, n) >= 0.39),
               prediction_small_arm=0.39907073974609375,
               scope="ratio-hitting probability; equals p_solve at nw=1 because "
                     "every true ratio hit is the same point and passes affine "
                     "verification; identical generative design to transfer_v26.py")
    out["p_set_law_expected"] = 1 - math.exp(-6 * out["distinct_mean"] / M)
    out["prediction_in_interval"] = bool(
        out["cp95_lower"] <= 0.39907073974609375 <= out["cp95_upper"])
    json.dump(out, open(HERE + "/transfer-v27.json", "w"), indent=1)
    print(json.dumps({k: out[k] for k in ("partial", "hits", "q", "cp95_lower",
          "cp95_upper", "distinct_mean", "p_pair_law", "p_set_law_expected",
          "floor_met", "wall_seconds")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
