"""Filtering simulation for the stnfs-archd relation system (stdlib only).

Synthetic relation hypergraph with the production per-ideal occurrence profile:
an ideal orbit of norm q on side s is contained in a random relation orbit with
probability min(1, 2*kappa_s(q)/q), kappa_s(q) = rho(u_s - log_B q)/rho(u_s)
(Dickman conditioning of B-smooth integers), where u_s = log N_s / log B.
Then CADO-style purge (singleton removal, clique removal to a target excess) and
Markowitz merge with |pivot| = 1 only, reporting rows N(W) at row-weight targets W
and the maximum coefficient magnitude.
"""
import math, random, sys, bisect, heapq
LN2 = math.log(2)
from decimal import Decimal, getcontext

_RHO_KMAX = 40; _RHO_NT = 160; _RHO_STEP = 0.001
def _build_rho_coefs():
    getcontext().prec = 150
    zero, one = Decimal(0), Decimal(1)
    coefs = [None, [one] + [zero] * (_RHO_NT - 1)]
    for k in range(2, _RHO_KMAX + 1):
        prev = coefs[k - 1]
        c = [zero] * _RHO_NT
        for m in range(_RHO_NT - 1):
            c[m + 1] = (prev[m] + m * c[m]) / (k * (m + 1))
        c[0] = prev[0] - sum(c[1:])
        coefs.append(c)
    return [None] + [[float(x) for x in c] for c in coefs[1:]]
_RHO_COEFS = _build_rho_coefs()
def rho_series(u):
    if u <= 1.0: return 1.0
    k = math.ceil(u)
    if k > _RHO_KMAX or k < 1: return 0.0
    xi = k - u; s = 0.0
    for ci in reversed(_RHO_COEFS[k]): s = s * xi + ci
    return s if s > 0.0 else 0.0
_RHO_NMAX = int((_RHO_KMAX - 0.001) / _RHO_STEP)
_RHO_TAB = [1.0 if i * _RHO_STEP <= 1.0 else rho_series(i * _RHO_STEP) for i in range(_RHO_NMAX + 2)]
def rho(u):
    if u <= 1.0: return 1.0
    x = u / _RHO_STEP
    if x >= _RHO_NMAX: return 0.0
    i = int(x); f = x - i
    return _RHO_TAB[i] * (1.0 - f) + _RHO_TAB[i + 1] * f

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

def kappa(u, t):
    return rho(u - t) / rho(u)

def build_system(L, us, excess, seed, log2U, t_dense=0.8, nt=4000, lsmall=None):
    """Production-profile hypergraph at scale 2^log2U sampled columns, with the dense part made explicit.
    Columns of a production relation orbit, by norm q = B^t:
      (a) t >= t_dense: sampled columns (per-side count 2^log2U/2, density B^t/t, rate lambda_s(t) = R 2kappa_s(t)/q);
      (b) q < 2^lsmall ("small primes", lsmall = log2(R'/25) so that these columns keep weight >= ~50 at the simulated scale): EXPLICIT columns, one per prime q and side (so their row combinations and
          coefficient growth are simulated exactly like the sampled ones); inclusion probability min(1, 2 kappa_s/q),
          coefficient = valuation ~ 1 + Geometric(1/q);
      (c) 2^lsmall <= q < B^t_dense ("mid band"): too heavy to be eliminated and too rare to be sampled at small scale;
          their expected count per row w_mid is carried as a per-row number that ADDS under row combination (their
          supports are disjoint with probability 1 - O(sum (2kappa/q)^2) ~ 1 - 2^-10), which is exact for the row weight
          and leaves their coefficients at the parents' valuations.
    Rows: dict col -> coefficient for (a) and (b); row[-1] = w_mid (float)."""
    rng = random.Random(seed)
    LiB = 2.0 ** L / (L * LN2)
    Rprod = excess * LiB
    Up = 2 ** log2U
    Rp = int(round(excess * Up))
    if lsmall is None:
        # explicit columns only where the small-scale weight R' 2kappa/q stays >= ~50 (never eliminable, never singleton):
        lsmall = max(8.0, math.log2(Rp / 25.0))
    ts = [t_dense + (1 - t_dense) * (i + 0.5) / nt for i in range(nt)]
    dens = [2.0 ** (L * (t - 1)) / t for t in ts]
    cum = []; tot = 0.0
    for d in dens: tot += d; cum.append(tot)
    cum = [c / tot for c in cum]
    cols = []; w = []
    for s_, u in enumerate(us):
        for _ in range(Up // 2):
            i = bisect.bisect_left(cum, rng.random())
            t = ts[min(i, nt - 1)] + (rng.random() - 0.5) * (1 - t_dense) / nt
            lam = Rprod * 2 * kappa(u, t) / 2.0 ** (L * t)
            cols.append((s_, t)); w.append(min(1.0, lam / Rp))
    # (b) explicit small-prime columns
    small = []   # (col index, prob, q)
    for s_, u in enumerate(us):
        for q in primes_upto(int(2 ** lsmall)):
            t = math.log(q) / (L * LN2)
            pr = min(1.0, 2 * kappa(u, t) / q)
            cols.append((s_, t)); w.append(pr); small.append((len(cols) - 1, pr, q))
    # (c) mid band expected count per row
    tmin = lsmall / L
    w_mid = 0.0
    m = 2000
    for s_, u in enumerate(us):
        for i in range(m):
            t = tmin + (t_dense - tmin) * (i + 0.5) / m
            w_mid += kappa(u, t) * ((t_dense - tmin) / m) / t
    U = len(cols)
    nsamp = (Up // 2) * len(us)
    dense = [i for i in range(nsamp) if w[i] >= 0.2]
    sparse = [i for i in range(nsamp) if w[i] < 0.2]
    scum = []; stot = 0.0
    for i in sparse:
        stot += w[i]; scum.append(stot)
    rows = []
    for _ in range(Rp):
        row = {}
        for i in dense:
            if rng.random() < w[i]: row[i] = 1
        k = rng_poisson(rng, stot)
        for _ in range(k):
            x = rng.random() * stot
            j = sparse[bisect.bisect_left(scum, x)]
            row[j] = row.get(j, 0) + 1
        for (j, pr, q) in small:
            if rng.random() < pr:
                v = 1
                while rng.random() < 1.0 / q: v += 1       # valuation
                row[j] = v
        row[-1] = w_mid
        rows.append(row)
    return cols, w, rows, w_mid

def rng_poisson(rng, lam):
    # Knuth for small lam, normal approx otherwise
    if lam < 30:
        L = math.exp(-lam); k = 0; p = 1.0
        while True:
            p *= rng.random()
            if p <= L: return k
            k += 1
    return max(0, int(round(rng.gauss(lam, math.sqrt(lam)))))

def _singletons(rows, alive, colw, nalive):
    changed = True
    while changed:
        changed = False
        for i, r in enumerate(rows):
            if not alive[i]: continue
            for c in r:
                if c >= 0 and colw[c] == 1:
                    alive[i] = False; nalive -= 1; changed = True
                    for cc in r:
                        if cc >= 0: colw[cc] -= 1
                    break
    return nalive

def rowweight(r):
    return (len(r) - 1) + r.get(-1, 0.0)      # explicit entries + mid-band count


def purge(rows, U, target_excess, verbose=True):
    """CADO-style purge.  (1) singleton removal to a fixpoint.  (2) clique removal: the
    weight-2 columns define a graph on rows; its connected components ("cliques") are
    removed heaviest-first (by total row weight), each followed by singleton cleanup,
    until rows - cols <= target_excess (an absolute row count)."""
    colw = [0] * U
    for r in rows:
        for c in r:
            if c >= 0: colw[c] += 1
    alive = [True] * len(rows)
    nalive = _singletons(rows, alive, colw, len(rows))
    ncols = sum(1 for x in colw if x > 0)
    if verbose: print("  purge/singletons: rows %d cols %d excess %d (%.3f)" % (nalive, ncols, nalive - ncols, (nalive - ncols) / max(ncols, 1)))
    rounds = 0
    while nalive - ncols > target_excess:
        rounds += 1
        # union-find over weight-2 columns
        parent = list(range(len(rows)))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        col_rows = {}
        for i, r in enumerate(rows):
            if not alive[i]: continue
            for c in r:
                if c >= 0 and colw[c] == 2:
                    col_rows.setdefault(c, []).append(i)
        for c, rs in col_rows.items():
            if len(rs) == 2:
                a, b = find(rs[0]), find(rs[1])
                if a != b: parent[a] = b
        comp = {}
        for i in range(len(rows)):
            if alive[i]:
                comp.setdefault(find(i), []).append(i)
        comps = sorted(comp.values(), key=lambda rs: -sum(rowweight(rows[i]) for i in rs))
        need = nalive - ncols - target_excess
        # each removed component reduces the excess by ~1 (k rows, k-1 connecting columns)
        # plus cascades; remove a batch of min(need/2, 1%) components, then recompute exactly
        nb = max(1, min(need // 2 + 1, max(1, len(comps) // 100)))
        removed = 0
        for rs in comps:
            if removed >= nb: break
            if len(rs) < 2 and need > 0 and nb > 1: pass
            for i in rs:
                if alive[i]:
                    alive[i] = False; nalive -= 1
                    for c in rows[i]:
                        if c >= 0: colw[c] -= 1
            removed += 1
        nalive = _singletons(rows, alive, colw, nalive)
        ncols = sum(1 for x in colw if x > 0)
    if verbose: print("  purge/cliques:    rows %d cols %d excess %d (%d rounds)" % (nalive, ncols, nalive - ncols, rounds))
    return [rows[i] for i in range(len(rows)) if alive[i]], ncols

def merge(rows, U, targets, kmax=32, cmax=255, verbose=True):
    """Markowitz merge: eliminate lightest columns first using a pivot row with coefficient +-1.
    Reports rows and max|coef| when the average row weight crosses each target W."""
    rows = [dict(r) for r in rows]
    alive = [True] * len(rows)
    colrows = {}
    for i, r in enumerate(rows):
        for c in r:
            if c >= 0: colrows.setdefault(c, set()).add(i)
    nrows = len(rows); ncols = len(colrows)
    totw = sum(rowweight(r) for r in rows)
    heap = [(len(s), c) for c, s in colrows.items()]
    heapq.heapify(heap)
    results = {}; tlist = sorted(targets); ti = 0
    maxc = 1
    eliminated_cols = 0
    while heap and ti < len(tlist):
        wgt, c = heapq.heappop(heap)
        if c not in colrows: continue
        s = colrows[c]
        if len(s) != wgt:
            heapq.heappush(heap, (len(s), c)); continue
        if wgt == 0:
            del colrows[c]; ncols -= 1; continue
        if wgt > kmax: break
        if wgt == 1:
            i = next(iter(s))
            alive[i] = False; nrows -= 1; totw -= rowweight(rows[i])
            for cc in rows[i]:
                if cc < 0: continue
                colrows[cc].discard(i)
                if not colrows[cc]: del colrows[cc]; ncols -= 1
                else: heapq.heappush(heap, (len(colrows[cc]), cc))
            rows[i] = {}
            continue
        piv = [i for i in s if abs(rows[i][c]) == 1]
        if not piv:
            # no unit pivot: skip this column (keep it)
            continue
        p = min(piv, key=lambda i: rowweight(rows[i]))
        pr = rows[p]; pc = pr[c]
        # coefficient bound: skip the column if any resulting entry would exceed cmax
        ok = True
        for i in s:
            if i == p: continue
            r = rows[i]; m = r[c] * pc
            for cc, v in pr.items():
                if cc < 0: continue
                if abs(r.get(cc, 0) - m * v) > cmax: ok = False; break
            if not ok: break
        if not ok:
            continue
        for i in list(s):
            if i == p: continue
            r = rows[i]; m = r[c] * pc   # r - m*pr kills column c (pc = +-1)
            r[-1] = r.get(-1, 0.0) + pr.get(-1, 0.0); totw += pr.get(-1, 0.0)   # mid band adds (disjoint supports)
            for cc, v in pr.items():
                if cc < 0: continue
                nv = r.get(cc, 0) - m * v
                if nv == 0:
                    if cc in r:
                        del r[cc]; colrows[cc].discard(i); totw -= 1
                        if not colrows[cc]: del colrows[cc]; ncols -= 1
                        else: heapq.heappush(heap, (len(colrows[cc]), cc))
                else:
                    if cc not in r:
                        colrows[cc].add(i); totw += 1
                        heapq.heappush(heap, (len(colrows[cc]), cc))
                    r[cc] = nv
                    if abs(nv) > maxc: maxc = abs(nv)
        # remove pivot row and column c
        alive[p] = False; nrows -= 1; totw -= rowweight(pr)
        for cc in pr:
            if cc == c or cc < 0: continue
            colrows[cc].discard(p)
            if not colrows[cc]: del colrows[cc]; ncols -= 1
            else: heapq.heappush(heap, (len(colrows[cc]), cc))
        rows[p] = {}
        del colrows[c]; ncols -= 1; eliminated_cols += 1
        avg = totw / max(nrows, 1)
        while ti < len(tlist) and avg >= tlist[ti]:
            results[tlist[ti]] = (nrows, ncols, avg, maxc, wgt)
            if verbose: print("  merge W>=%4d: rows %d cols %d avg %.1f maxcoef %d (col weight %d)" % (tlist[ti], nrows, ncols, avg, maxc, wgt))
            ti += 1
    if ti < len(tlist):
        avg = totw / max(nrows, 1)
        if verbose: print("  merge stopped: rows %d cols %d avg %.1f maxcoef %d" % (nrows, ncols, avg, maxc))
        for t in tlist[ti:]:
            results[t] = (nrows, ncols, avg, maxc, -1)
    return results

def run(L, uf, ug, excess, seed, log2U, pexc=0.002, kmax=32, targets=(100, 150, 200, 300, 400, 600), verbose=True):
    cols, w, rows, w_mid = build_system(L, (uf, ug), excess, seed, log2U)
    Usamp = sum(1 for (s_, t) in cols if t >= 0.8)          # sampled (top) columns: the scale reference U
    U = len(cols)
    lam = [0.0] * U
    for r in rows:
        for c in r:
            if c >= 0: lam[c] += 1
    zero = sum(1 for x in lam[:Usamp] if x == 0)
    w0 = sum(rowweight(r) for r in rows) / len(rows)
    if verbose:
        print("L=%g uf=%g ug=%g excess=%g seed=%d: U=%d (+%d explicit small-prime columns) R=%d w0(full)=%.1f w_mid=%.1f zero-cols %.4f" % (L, uf, ug, excess, seed, Usamp, U - Usamp, len(rows), w0, w_mid, zero / Usamp))
    rows2, ncols = purge(rows, U, max(2, int(pexc * Usamp)), verbose)
    res = merge(rows2, U, targets, kmax=kmax, verbose=verbose)
    out = {}
    for W, (n, nc, avg, maxc, cw) in sorted(res.items()):
        out[W] = (Usamp / n, n, maxc, avg)
        if verbose: print("  W=%d: U/N = %.2f  (N=%d, maxcoef=%d, full row weight %.1f)" % (W, Usamp / n, n, maxc, avg))
    return Usamp, out

if __name__ == "__main__":
    L = float(sys.argv[1]); uf = float(sys.argv[2]); ug = float(sys.argv[3]); exc = float(sys.argv[4]); seed = int(sys.argv[5]); log2U = int(sys.argv[6])
    kmax = int(sys.argv[7]) if len(sys.argv) > 7 else 32
    run(L, uf, ug, exc, seed, log2U, kmax=kmax)


# ---- calibration against the De Micheli-Gaudry-Pierrot 521-bit F_{p^6} record (ASIACRYPT 2021) ----
# Reported: factor base 15.21M ideals, 13.63M unique relations (0.896 per ideal), 87% of the ideals appear,
# singleton removal leaves 5.21M rows (38%), merge to density 100 / 150 leaves 1.73M / 1.51M rows (33% / 29% of
# the purged matrix).  Their norms: ~80 bits (f1 = x^4+1 side) and ~170 bits (f2 side, special-q of ~23 bits
# removed), B = 2^27, no automorphism quotient.  u-values are therefore ~3.0 and ~5.4 (our best estimate).
def calibrate_dgp21(log2U=16, seed=1, verbose=True):
    L = 27.0; us = (3.0, 5.4); excess = 0.896
    LiB = 2.0 ** L / (L * LN2)
    rng = random.Random(seed)
    Up = 2 ** log2U; Rp = int(round(excess * Up))
    t_dense = 0.8; nt = 4000
    ts = [t_dense + (1 - t_dense) * (i + 0.5) / nt for i in range(nt)]
    dens = [2.0 ** (L * (t - 1)) / t for t in ts]
    cum = []; tot = 0.0
    for d in dens: tot += d; cum.append(tot)
    cum = [c / tot for c in cum]
    cols = []; w = []
    for s, u in enumerate(us):
        for _ in range(Up // 2):
            i = bisect.bisect_left(cum, rng.random())
            t = ts[min(i, nt - 1)] + (rng.random() - 0.5) * (1 - t_dense) / nt
            lam = excess * (2 * LiB) * kappa(u, t) / 2.0 ** (L * t)     # no quotient: one column per ideal, rate kappa/q
            cols.append((s, t)); w.append(min(1.0, lam / Rp))
    U = len(cols)
    dense = [i for i in range(U) if w[i] >= 0.2]
    sparse = [i for i in range(U) if w[i] < 0.2]
    scum = []; stot = 0.0
    for i in sparse:
        stot += w[i]; scum.append(stot)
    rows = []
    for _ in range(Rp):
        row = {}
        for i in dense:
            if rng.random() < w[i]: row[i] = 1
        k = rng_poisson(rng, stot)
        for _ in range(k):
            x = rng.random() * stot
            j = sparse[bisect.bisect_left(scum, x)]
            row[j] = row.get(j, 0) + 1
        rows.append(row)
    lam = [0] * U
    for r in rows:
        for c in r: lam[c] += 1
    appear = 1 - sum(1 for x in lam if x == 0) / U
    rows2, ncols = purge(rows, U, 10 ** 9, verbose=False)      # singleton removal only (no excess to spend)
    res = merge(rows2, U, (100, 150), kmax=64, verbose=False)
    out = dict(appearing=appear, purged_frac=len(rows2) / len(rows), merge100=res[100][0] / len(rows2), merge150=res[150][0] / len(rows2))
    if verbose:
        print("DGP21 calibration (scale 2^%d): ideals appearing %.3f (reported 0.87); rows after singletons %.3f (reported 0.38); "
              "merge to density 100: %.3f, 150: %.3f of purged rows (reported 0.33, 0.29)" % (log2U, appear, out['purged_frac'], out['merge100'], out['merge150']))
    return out
