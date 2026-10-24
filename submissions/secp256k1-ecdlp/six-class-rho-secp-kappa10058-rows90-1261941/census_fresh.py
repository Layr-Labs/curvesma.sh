"""Fresh-table hierarchical census of the six-class rho walk (offline evidence).

Each table instance is one cold-attack table (fresh random table points and a
fresh secret), and B independent collision searches are run on it with fresh
start coins.  This measures both the trial-to-trial and the table-to-table
variation, which the pooled one-table census does not.  Stdlib only.
"""
import argparse, json, math, random, time
import census as C

PS = [0.390, 0.392, 0.394, 0.396, 0.398]

def quantile(sorted_vals, p):
    n = len(sorted_vals)
    x = p * (n - 1)
    i = int(math.floor(x)); f = x - i
    if i + 1 >= n: return sorted_vals[-1]
    return sorted_vals[i] * (1 - f) + sorted_vals[i + 1] * f

def one_table(E, rng, R=2048, dp=5, trials=1000):
    r = E.r
    xs = rng.randrange(1, r)
    W = C.Walker(E, R, dp, E.fmul(xs), rng, lambda c, d: E.fmul((c + d * xs) % r))
    tcs = []
    for _ in range(trials):
        first_seen = {}
        st = {"visits": 0, "tc": None, "trail": 0}
        def on_class(K):
            st["visits"] += 1
            prev = first_seen.get(K)
            if prev is None: first_seen[K] = st["trail"]
            elif prev != st["trail"] and st["tc"] is None: st["tc"] = st["visits"]
        while st["tc"] is None:
            a, b = rng.randrange(r), rng.randrange(r)
            P = E.fmul((a + b * xs) % r)
            if P is None: continue
            st["trail"] += 1
            W.run_trail(P, on_class)
        tcs.append(st["tc"])
    return tcs

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--curve", required=True)
    ap.add_argument("--tables", type=int, default=200)
    ap.add_argument("--trials", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    t0 = time.time()
    rng = random.Random(a.seed)
    E = C.make_toy(a.curve, rng)
    N = (E.r - 1) // 6
    sq = math.sqrt(N)
    rows = []
    for t in range(a.tables):
        tcs = one_table(E, rng, trials=a.trials)
        s = sorted(v / sq for v in tcs)
        rows.append({"mean": sum(s) / len(s),
                     "q": {("%.3f" % p): quantile(s, p) for p in PS},
                     "nv": len(s)})
        if t % 20 == 0:
            print("table %d/%d  %.1fs" % (t, a.tables, time.time() - t0), flush=True)
    json.dump({"curve": a.curve, "q": E.q, "r": E.r, "N": N, "R": 2048,
               "dp_bits": 5, "tables": a.tables, "trials": a.trials,
               "seed": a.seed, "ps": PS, "rows": rows,
               "seconds": round(time.time() - t0, 1)}, open(a.out, "w"))

if __name__ == "__main__":
    main()
