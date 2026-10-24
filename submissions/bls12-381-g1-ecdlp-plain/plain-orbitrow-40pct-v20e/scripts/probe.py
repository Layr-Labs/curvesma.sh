import sys, math, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, qwalk
from collections import Counter

SQRTS = 23700.0
LAMF = math.log(1 / 0.6097)
KAPPA = 1.015
d, bt, w = 5, 100, 64
U = qwalk.Universe(int(round(6 * SQRTS * SQRTS)), 1 << 16, seed=20250)
q = U.q
N = U.classes
s_ref = math.sqrt(2 * N * LAMF)
budget = int(math.ceil(KAPPA * s_ref))
per = budget // w
lanes = bt * w
rng = np.random.default_rng(99)
xt = rng.integers(1, q - 1, size=bt).astype(np.uint64)
xv = np.repeat(xt, w).astype(np.uint64)
a = rng.integers(0, q, size=lanes).astype(np.uint64)
cc = rng.integers(0, q, size=lanes).astype(np.uint64)
t0 = np.array([(int(ai) + int(ci) * int(xi)) % q for ai, ci, xi in
               zip(a.tolist(), cc.tolist(), xv.tolist())], dtype=np.uint64)
canon, e, _ = U.canon(t0)
a = (e * a) % q
cc = (e * cc) % q
trace = np.empty((per, lanes), dtype=np.uint64)
rl, rcls, rw, ra, rcc_l = [], [], [], [], []
for r in range(per):
    trace[r] = canon
    canon, a, cc, hit = qwalk.step(U, canon, a, cc, xv, d, want_dp=True)
    if hit.any():
        idx = np.flatnonzero(hit)
        rl.append(np.full(idx.size, r))
        rcls.append(canon[idx])
        ra.append(a[idx])
        rcc_l.append(cc[idx])
        rw.append(idx)
        a0 = rng.integers(0, q, size=idx.size).astype(np.uint64)
        c0 = rng.integers(0, q, size=idx.size).astype(np.uint64)
        s0 = np.array([(int(x) + int(y) * int(z)) % q for x, y, z in
                       zip(a0.tolist(), c0.tolist(), xv[idx].tolist())], dtype=np.uint64)
        cn2, e2, _ = U.canon(s0)
        canon[idx] = cn2
        a[idx] = (e2 * a0) % q
        cc[idx] = (e2 * c0) % q
rlv = np.concatenate(rl); rc = np.concatenate(rcls); rwv = np.concatenate(rw)
rav = np.concatenate(ra); rccv = np.concatenate(rcc_l)
tri = rwv // w
c = Counter(zip(tri.tolist(), rc.tolist()))
dups = {k: v for k, v in c.items() if v >= 2}
print("records", rlv.size, "| (trial,class) with >=2 recs", len(dups),
      "| trials with a dup", len({k[0] for k in dups}))

ctr = np.repeat(np.arange(bt), w)
wid = np.tile(np.arange(lanes), per)
key = (ctr[None, :] * np.int64(q + 1) + trace.astype(np.int64)).reshape(-1)
lvl = np.repeat(np.arange(per)[:, None], lanes, axis=1).reshape(-1)
o = np.argsort(key, kind="stable")
sk = key[o]
dup = np.flatnonzero(sk[1:] == sk[:-1])
merge = set()
sw = wid[o]
for kk in np.unique(sk[dup]).tolist():
    g0 = int(np.searchsorted(sk, kk, "left"))
    g1 = int(np.searchsorted(sk, kk, "right"))
    if len(set(sw[g0:g1].tolist())) > 1:
        merge.add(int(kk // (q + 1)))
duptrials = {k[0] for k in dups}
print("merge trials", len(merge), "| dup trials", len(duptrials),
      "| merge & ~dup", len(merge - duptrials), "| dup & ~merge", len(duptrials - merge))

# among (trial,class) groups with >=2 records: how many are from 2 walkers, and how
# many solve correctly?
ok = bad = dead = diffwalkers = 0
for (tr, cls), nrec in dups.items():
    sel = np.flatnonzero((tri == tr) & (rc == cls))
    ws = rwv[sel]
    if len(set(ws.tolist())) < 2:
        continue
    diffwalkers += 1
    k1 = int(sel[0])
    pos = 1
    while pos < len(sel) and rwv[int(sel[pos])] == rwv[k1]:
        pos += 1
    if pos >= len(sel):
        continue
    k2 = int(sel[pos])
    den = (int(rccv[k2]) - int(rccv[k1])) % q
    if den == 0:
        dead += 1
        continue
    xr = (int(rav[k1]) - int(rav[k2])) * pow(den, -1, q) % q
    if xr == int(xt[tr]):
        ok += 1
    else:
        bad += 1
print("dup groups with 2 walkers", diffwalkers, "| solve ok", ok, "| solve wrong", bad,
      "| denom dead", dead)
# check invariant a + c*x == canon at a record
k = 0
print("invariant check", (int(rav[k]) + int(rccv[k]) * int(xt[tri[k]])) % q == int(rc[k]),
      "| class at record", int(rc[k]), "| q", q)
