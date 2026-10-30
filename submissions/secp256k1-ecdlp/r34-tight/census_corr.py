"""Borrow census of the three branchy fe_sub sites (dx = x_T - x_W, dy = y_T - y_W, u = x3 - x_W) on real secp256k1
under the production walk rule: six classes, canonical index j = argmin(x, bx, b^2x), key K = the even one of y and p - y,
h = bits 1..rbits of K, Tab[h][3s + j] = (-1)^s phi^-j(M_h), the stored member -(W + T), checkpoint anchors every 8192
operations and compares every 256 (phase phi per trail), traversal of <= 32 steps, doubling escape from the
minimum-fp member, history cleared. A doubling round still runs its slot's addition: its borrows are counted.

Table rows M_h are uniform random points derived from (seed, h) by BLAKE2b (try-and-increment on x, sign from the
hash), so a 2^26-row table needs no memory. Walk starts are uniform random points from the worker's own RNG.

usage: census_corr.py rbits steps_per_walk walks_per_worker workers seed out.json
"""
import hashlib, json, os, random, sys, time
from multiprocessing import Pool

try:
    from gmpy2 import mpz, powmod, invert
except ImportError:  # pure-Python fallback (slower)
    mpz = int
    powmod = pow
    def invert(a, m):
        return pow(a, -1, m)

P = 2 ** 256 - 2 ** 32 - 977
BETA = int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee", 16)
BETA2 = BETA * BETA % P
M64 = (1 << 64) - 1
C0 = 2 ** 32 + 977
STRIDE, PERIOD, TRAV_CAP = 256, 8192, 32
Pm, BETAm, BETA2m, SQE = mpz(P), mpz(BETA), mpz(BETA2), mpz((P + 1) // 4)


def row_point(seed, h):
    """Uniform random point M_h (x and y), redrawn while an orbit x-value is >= p - C0 (production rule)."""
    ctr = 0
    while True:
        d = hashlib.blake2b(b"%d/%d/%d" % (seed, h, ctr), digest_size=33).digest()
        ctr += 1
        x = mpz(int.from_bytes(d[:32], "little") % P)
        y2 = (x * x * x + 7) % Pm
        y = powmod(y2, SQE, Pm)
        if y * y % Pm != y2 or y == 0:
            continue
        if max(x, BETAm * x % Pm, BETA2m * x % Pm) >= P - C0:
            continue
        return x, (y if d[32] & 1 else Pm - y)


def fp(K):
    return (K & M64) | ((K >> 128) << 64)  # limbs 0, 2, 3 of K


def walk(args):
    """One trail sequence in one slot: steps executed slot-steps. Returns per-site borrow counts and diagnostics."""
    seed, rbits, steps, wid, cache = args
    R = 1 << rbits
    rr = random.Random(seed * 1000003 + wid)
    tab = cache.get(seed) if cache is not None else None

    def M(h):
        return tab[h] if tab is not None else row_point(seed, h)

    def rand_point():
        while True:
            x = mpz(rr.randrange(P))
            y2 = (x * x * x + 7) % Pm
            y = powmod(y2, SQE, Pm)
            if y * y % Pm == y2 and y != 0:
                return x, (y if rr.random() < 0.5 else Pm - y)

    def arrive(W):
        x, y = W
        v1, v2 = BETAm * x % Pm, BETA2m * x % Pm
        j = 0 if (x < v1 and x < v2) else (1 if v1 < v2 else 2)
        s = int(y & 1)
        K = Pm - y if s else y
        return int(K), j, s

    def target(W, K, j, s):
        xm, ym = M((K >> 1) & (R - 1))
        xt = xm if j == 0 else (BETA2m * xm % Pm if j == 1 else BETAm * xm % Pm)  # beta^-j x_M
        return xt, (Pm - ym if s else ym)

    def add(W, T):
        (x, y), (xt, yt) = W, T
        lam = (yt - y) * invert(xt - x, Pm) % Pm
        x3 = (lam * lam - x - xt) % Pm
        y3 = (lam * (x - x3) - y) % Pm
        return x3, y3

    cnt = dict(steps=0, dx=0, dy=0, u=0, doublings=0, traversal_steps=0, trails=1, exceptional=0,
               j=[0, 0, 0], dx_by_j=[0, 0, 0], xw_sum_by_j=[0.0, 0.0, 0.0])
    prefix = None

    def slot_step(W, K, j, s):
        """The slot's addition (always run): count borrows, return -(W + T)."""
        T = target(W, K, j, s)
        if T[0] == W[0]:
            return None
        x3, y3 = add(W, T)
        bdx, bdy, bu = T[0] < W[0], T[1] < W[1], x3 < W[0]
        cnt["dx"] += bdx; cnt["dy"] += bdy; cnt["u"] += bu
        cnt["j"][j] += 1; cnt["dx_by_j"][j] += bdx; cnt["xw_sum_by_j"][j] += float(W[0]) / P
        cnt["steps"] += 1
        return x3, Pm - y3

    W = rand_point()
    length, phase, cp = 0, rr.randrange(PERIOD), None
    while cnt["steps"] < steps:
        if prefix is None and cnt["steps"] >= 4096:
            prefix = [cnt["steps"], cnt["dx"], cnt["dy"], cnt["u"]]
        K, j, s = arrive(W)
        if (K >> 64) & M64 == 0:  # distinguished (2^-64): new trail
            W, length, phase, cp = rand_point(), 0, rr.randrange(PERIOD), None
            cnt["trails"] += 1
            continue
        if cp is not None and (length + phase) % STRIDE == 0 and fp(K) == cp:
            best, Q, trav, ok = (fp(K), W), W, 0, True
            while True:
                trav += 1
                if trav > TRAV_CAP:
                    ok = False
                    break
                kq, jq, sq = arrive(Q)
                Q = slot_step(Q, kq, jq, sq)
                cnt["traversal_steps"] += 1
                length += 1
                if Q is None:
                    ok = False
                    break
                fq = fp(arrive(Q)[0])
                if fq == fp(K):
                    break
                if fq < best[0]:
                    best = (fq, Q)
            if not ok:
                W, length, phase, cp = rand_point(), 0, rr.randrange(PERIOD), None
                cnt["trails"] += 1
                continue
            Wm = best[1]
            km, jm, sm = arrive(Wm)
            slot_step(Wm, km, jm, sm)  # the doubling round's discarded addition
            x, y = Wm
            lam = 3 * x * x * invert(2 * y, Pm) % Pm
            x3 = (lam * lam - 2 * x) % Pm
            W = (x3, (lam * (x - x3) - y) % Pm)
            cnt["doublings"] += 1
            cp = None
            length += 1
            continue
        if (length + phase) % PERIOD == 0:
            cp = fp(K)
        Wn = slot_step(W, K, j, s)
        if Wn is None:
            cnt["exceptional"] += 1
            W, length, phase, cp = rand_point(), 0, rr.randrange(PERIOD), None
            cnt["trails"] += 1
            continue
        W = Wn
        length += 1
    cnt["prefix4096"] = prefix
    return cnt


_CACHE = {}


def init(seed, rbits):
    if rbits <= 16:
        _CACHE[seed] = [row_point(seed, h) for h in range(1 << rbits)]


def job(a):
    return walk(a[:4] + (_CACHE,))


def main():
    rbits, steps, wpw, workers, seed, outp = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
    t0 = time.time()
    rows = []
    # one table per seed (seed + k, k < tables); each table's walks are independent trails on it
    tasks = [(seed + k, rbits, steps, w) for k in range(workers) for w in range(wpw)]
    res = []
    if rbits <= 16:  # small tables are built once per process
        for k in range(workers):
            with Pool(workers, initializer=init, initargs=(seed + k, rbits)) as p2:
                res += p2.map(job, [t for t in tasks if t[0] == seed + k])
    else:
        with Pool(workers) as pool:
            res = pool.map(job, tasks, chunksize=1)
    for t, c in zip(tasks, res):
        rows.append(dict(table=t[0], walk=t[3], **c))
    tot = {k: sum(r[k] for r in rows) for k in ("steps", "dx", "dy", "u", "doublings", "traversal_steps", "trails", "exceptional")}
    n = tot["steps"]
    summ = {}
    for site in ("dx", "dy", "u"):
        q = tot[site] / n
        ss = sum((r[site] - q * r["steps"]) ** 2 for r in rows)
        se = (len(rows) / (len(rows) - 1) * ss) ** 0.5 / n
        summ[site] = dict(rate=q, se=se)
    out = dict(config=dict(rbits=rbits, steps_per_walk=steps, walks_per_table=wpw, tables=workers, seed=seed),
               totals=tot, summary=summ, seconds=time.time() - t0, rows=rows)
    with open(outp, "w") as fh:
        json.dump(out, fh)
    print(json.dumps(dict(config=out["config"], totals=tot, summary=summ, seconds=round(out["seconds"], 1))))


if __name__ == "__main__":
    main()
