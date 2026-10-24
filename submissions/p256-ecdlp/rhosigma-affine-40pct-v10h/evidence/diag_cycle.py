#!/usr/bin/env python3
"""diag_cycle.py — exact-arithmetic reproduction of the p256 v10f 2^40 walk map.

Builds the same 2048-row canonical table, walks the same map, and (a) reports the
cycle structure exactly and (b) verifies the cycle relation in bigint arithmetic,
so the attractor diagnosis cannot be an artifact of the C field code.
"""
import sys

P = 1099511629627
Q = 1099513640503
NB = 2048
GX, GY = 2, 185782196093
A = 0

M64 = (1 << 64) - 1

def splitmix(state):
    state = (state + 0x9E3779B97F4A7C15) & M64
    z = state
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M64
    return state, z ^ (z >> 31)

def mix64(v):
    v = (v * 0x9E3779B97F4A7C15) & M64
    v ^= v >> 29
    v = (v * 0xBF58476D1CE4E5B9) & M64
    v ^= v >> 32
    return v

def add(Ap, Bp):
    if Ap is None: return Bp
    if Bp is None: return Ap
    (x1, y1), (x2, y2) = Ap, Bp
    if x1 == x2:
        if (y1 + y2) % P == 0: return None
        lam = (3 * x1 * x1) * pow(2 * y1, -1, P) % P
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, P) % P
    x3 = (lam * lam - x1 - x2) % P
    y3 = (lam * (x1 - x3) - y1) % P
    return (x3, y3)

def mul(k, Pnt=(GX, GY)):
    R = None
    while k:
        if k & 1: R = add(R, Pnt)
        Pnt = add(Pnt, Pnt)
        k >>= 1
    return R

def canon(Pnt):
    if Pnt is None: return None
    x, y = Pnt
    return (x, y if y % 2 == 0 else P - y)

def build(seed=0xABCDEF0123456789):
    s = seed; T = []; L = []
    for j in range(NB):
        while True:
            s, t = splitmix(s); t = t % (Q - 1) + 1
            if t not in L and (Q - t) % Q not in L: break
        Pt = mul(t)
        if Pt[1] % 2:
            Pt = (Pt[0], P - Pt[1]); t = Q - t
        T.append(Pt); L.append(t)
    return T

def idx(x, mode):
    if mode == 0: return x & (NB - 1)
    if mode == 1: return mix64(x) & (NB - 1)
    if mode == 4: return (x ^ (x >> 20)) & (NB - 1)
    raise SystemExit("mode")

def step(S, T, mode):
    j = idx(S[0], mode)
    Tj = T[j]
    R = add(S, S) if S[0] == Tj[0] else add(S, Tj)
    return canon(R), j

def main():
    mode = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    nst = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    T = build()
    print("table distinct rows:", len(set(T)), "min x size 2^%.2f" % (max(t[0] for t in T).bit_length()))
    s = 777
    hist = {}
    for t in range(nst):
        s, k = splitmix(s); k = k % (Q - 1) + 1
        P0 = canon(mul(k))
        seen = {P0: 0}
        seq = [P0]
        cur = P0
        cyc = None
        for n in range(1, 6000):
            cur, j = step(cur, T, mode)
            if cur in seen:
                cyc = (seen[cur], n - seen[cur])
                break
            seen[cur] = n
            seq.append(cur)
        hist[cyc[1] if cyc else None] = hist.get(cyc[1] if cyc else None, 0) + 1
        if t < 4:
            print(f"start {t}: cycle entry mu={cyc[0] if cyc else '-'} lam={cyc[1] if cyc else '-'}")
            if cyc:
                mu, lam = cyc
                c0 = seq[mu]
                c1, j0 = step(c0, T, mode)
                print(f"   cycle p.x={c0[0]} q.x={c1[0]}  lam={lam}")
                # exact relation: q + T_j' = +-p
                j1 = idx(c1[0], mode)
                lhs = add(c1, T[j1])
                print(f"   q+T[j'] == p ? {lhs == c0}   == -p ? {lhs == (c0[0], P - c0[1])}")
                print(f"   T[j'] == -(2p+T[j]) ? {T[j1] == add(canon(mul((2*(0)))), None) if False else (T[j1] == (lambda r: (r[0], P - r[1]))(add(add(c0, c0), T[j0])))}")
    print("lam histogram:", hist)

main()
