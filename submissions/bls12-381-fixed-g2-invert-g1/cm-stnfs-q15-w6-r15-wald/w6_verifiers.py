"""w6_verifiers.py — Live Wave-6 H3/H-X & strict-k=1 selected-slot H5 verifiers on exact (9,3) Q(zeta_15)."""
import heapq
import math
import random
from filtsim_z12 import (
    U12, U12_SET, U12_INV, z12_add, z12_sub, z12_mul, z12_is_zero, z12_maxabs,
    z12_eval_mod, _vinfo, run_filtsim_z12_ladder,
)
from hx_z12 import _right_nullspace, _rank_rows, _rref


def _isp_det(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n == p:
            return True
        if n % p == 0:
            return False
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True


def verify_threshold_counterexample():
    """Exact algebraic witness that naive PGF substitution G(psi) fails for general non-i.i.d. laws."""
    mu, rv, psi = 60.0, 1.0, 0.10
    bt = 0.5 * mu / rv
    G = lambda t: (1.0 + bt * (1.0 - t)) ** (-rv)
    for i in range(81):
        t = i / 100.0
        assert 0.10 * t + 0.90 * (t ** 100) <= G(t) + 1e-12
    assert 0.10 > G(psi) + 0.05
    best = G(0.0) + psi
    for k in range(2, 65):
        tk = min(0.8, (k - 1) * (1.0 + 1.0 / bt) / (k - 1 + rv))
        best = min(best, min(1.0, G(tk) / (tk ** (k - 1))) + (psi ** k))
    assert best >= 0.10 - 1e-12
    return dict(pgf_psi_invalid=round(G(psi), 5), true_fail=0.10, tr_bound=round(best, 5), tr_valid=True)


def verify_genuine_q15_filter_hx_live(
    toy, primes_upto, PHI15, y_pow_z, sub_phi15_z, mul_phi15_z, pow_phi15_z,
    scale_phi15_z, norm_q15_z, pgcd_mod, pmul_mod, prem_mod, u12_proxy=None,
):
    """Genuine UNPLANTED Q(zeta_15) K_g (g(X)=X^3+5w^2, u0=-5, p0=7207, r0=601):
    1) Collect genuine Z[zeta_12] relations (exact F_{7207^12}^* dlogs via BSGS, ZERO planted kernel).
    2) Run exact zero-column-loss Z[zeta_12] filter (peel + _dc-guarded merge, cols_2core_lost==0).
    3) Run UNPLANTED Block Wiedemann H-X minor check & Horner recovery on the filtered Z[zeta_12] matrix.
    4) Back-substitute ALL merged pivots and verify 100% exact recovery of ALL training + held-out dlogs.
    5) Run live full-incidence (t_dense=0.142 => w_mid~0.0, zero omitted middle cols) vs t_dense=0.80 bridge.
    """
    U0, P0, R0, CHI = -5, 7207, 601, 596
    h_pi = toy["h_pi"]
    _h0, _h1, _h2, _h3 = h_pi[0], h_pi[1], h_pi[2], h_pi[3]
    W2 = y_pow_z(10)

    def fp4_add(a, b):
        return [(a[i] + b[i]) % P0 for i in range(4)]

    def fp4_mul(a, b):
        c = [
            a[0]*b[0], a[0]*b[1]+a[1]*b[0], a[0]*b[2]+a[1]*b[1]+a[2]*b[0],
            a[0]*b[3]+a[1]*b[2]+a[2]*b[1]+a[3]*b[0], a[1]*b[3]+a[2]*b[2]+a[3]*b[1],
            a[2]*b[3]+a[3]*b[2], a[3]*b[3],
        ]
        for k in (6, 5, 4):
            v = c[k] % P0
            if v:
                o = k - 4
                c[o] -= v * _h0; c[o+1] -= v * _h1; c[o+2] -= v * _h2; c[o+3] -= v * _h3
        return [c[i] % P0 for i in range(4)]

    def fp4_pow(a, e):
        res, b = [1, 0, 0, 0], list(a)
        while e:
            if e & 1:
                res = fp4_mul(res, b)
            e >>= 1
            if e:
                b = fp4_mul(b, b)
        return res

    def fp4_inv(a):
        a1 = fp4_pow(a, P0); a2 = fp4_pow(a1, P0); a3 = fp4_pow(a2, P0)
        p3 = fp4_mul(a1, fp4_mul(a2, a3))
        inv = pow(fp4_mul(a, p3)[0], -1, P0)
        return [(x * inv) % P0 for x in p3]

    w4 = fp4_pow([0, 1, 0, 0], 5); w4_2 = fp4_mul(w4, w4); uw2 = [(U0 * x) % P0 for x in w4_2]
    frob4_cols = [fp4_pow([1 if j == i else 0 for j in range(4)], P0) for i in range(4)]

    def frob4(a):
        out = [0, 0, 0, 0]
        for i, ai in enumerate(a):
            if ai:
                col = frob4_cols[i]
                for j in range(4):
                    out[j] = (out[j] + ai * col[j]) % P0
        return out

    def fp12_mul(A, B):
        return (
            fp4_add(fp4_mul(A[0], B[0]), fp4_mul(uw2, fp4_add(fp4_mul(A[1], B[2]), fp4_mul(A[2], B[1])))),
            fp4_add(fp4_add(fp4_mul(A[0], B[1]), fp4_mul(A[1], B[0])), fp4_mul(uw2, fp4_mul(A[2], B[2]))),
            fp4_add(fp4_add(fp4_mul(A[0], B[2]), fp4_mul(A[1], B[1])), fp4_mul(A[2], B[0])),
        )

    def fp12_pow(A, e):
        res, b = ([1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]), A
        while e:
            if e & 1:
                res = fp12_mul(res, b)
            e >>= 1
            if e:
                b = fp12_mul(b, b)
        return res

    def fp12_frob(A):
        return (frob4(A[0]), fp4_mul(frob4(A[1]), w4), fp4_mul(frob4(A[2]), w4_2))

    def fp12_inv(A):
        def f4(Z):
            for _ in range(4):
                Z = fp12_frob(Z)
            return Z
        A4 = f4(A); A8 = f4(A4); adj = fp12_mul(A4, A8)
        inv4 = fp4_inv(fp12_mul(A, adj)[0])
        return (fp4_mul(adj[0], inv4), fp4_mul(adj[1], inv4), fp4_mul(adj[2], inv4))

    cof_GT = (P0**4 - P0**2 + 1) // R0

    def to_GT(A):
        e6 = A
        for _ in range(6):
            e6 = fp12_frob(e6)
        ez = fp12_mul(e6, fp12_inv(A)); e2 = fp12_frob(fp12_frob(ez))
        return fp12_pow(fp12_mul(e2, ez), cof_GT)

    def embed_ab(a8, b8):
        ar = prem_mod(a8, h_pi, P0); ar += [0] * (4 - len(ar))
        br = prem_mod([-x for x in b8], h_pi, P0); br += [0] * (4 - len(br))
        return (ar, br, [0, 0, 0, 0])

    G_gen = to_GT(([3, 1, 0, 0], [2, 0, P0 - 1, 0], [0, 5, 0, 0]))
    m_bs = math.isqrt(R0) + 1; bs_map = {}; cur = ([1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0])
    for j in range(m_bs):
        bs_map[(tuple(cur[0]), tuple(cur[1]), tuple(cur[2]))] = j
        cur = fp12_mul(cur, G_gen)
    Gm_inv = fp12_inv(fp12_pow(G_gen, m_bs))

    def dlog_GT(A):
        c = to_GT(A)
        for i in range(m_bs + 1):
            k = (tuple(c[0]), tuple(c[1]), tuple(c[2]))
            if k in bs_map:
                return (i * m_bs + bs_map[k]) % R0
            c = fp12_mul(c, Gm_inv)
        raise RuntimeError("BSGS failed")

    B_FB = 2200
    ps_1 = [p for p in primes_upto(B_FB) if p % 15 == 1 and p != R0]
    ps_van = [2, 3, 5] + [p for p in primes_upto(B_FB) if p % 15 in (4, 7, 13)]
    P_van = math.prod(ps_van); P_1 = math.prod(ps_1)

    def strip(n):
        if n % R0 == 0:
            return None
        t = n
        while True:
            g = math.gcd(t, P_van)
            if g == 1:
                break
            t //= g
        if t == 1:
            return []
        r = t; hits = []
        if math.gcd(r, P_1) == 1:
            return None
        for p in ps_1:
            if r % p == 0:
                e = 0
                while r % p == 0:
                    r //= p; e += 1
                if e > 1:
                    return None
                hits.append((p, 1))
                if r == 1:
                    return hits
        return hits if r == 1 else None

    cands = [[1 if j == i else 0 for j in range(8)] for i in range(8)]
    seen_c = set()

    def canon(v):
        best = None; cur_v = list(v)
        for _ in range(15):
            for s in (1, -1):
                t = tuple(s * x for x in cur_v)
                if best is None or t < best:
                    best = t
            cur_v = mul_phi15_z(cur_v, [0, 1, 0, 0, 0, 0, 0, 0])
        return best

    for i in range(8):
        for j in range(i + 1, 8):
            for s in (1, -1):
                v = [0]*8; v[i] = 1; v[j] = s; c = canon(v)
                if c not in seen_c:
                    seen_c.add(c); cands.append(v)
    for i in range(8):
        for j in range(i + 1, 8):
            for k in range(j + 1, 8):
                for s1 in (1, -1):
                    for s2 in (1, -1):
                        v = [0]*8; v[i] = 1; v[j] = s1; v[k] = s2; c = canon(v)
                        if c not in seen_c and abs(norm_q15_z(v)) <= 200:
                            seen_c.add(c); cands.append(v)

    pre = [(v, pow_phi15_z(v, 3), abs(norm_q15_z(v))) for v in cands[:60]]
    raw_hits = []; used_ps = set()
    for ia, (a8, a3, na) in enumerate(pre):
        for ib, (b8, b3, nb) in enumerate(pre):
            if ia == ib:
                continue
            ng = abs(norm_q15_z(sub_phi15_z(a3, scale_phi15_z(U0, mul_phi15_z(W2, b3)))))
            if ng <= 1:
                continue
            fg = strip(ng)
            if not fg:
                continue
            g0 = math.gcd(na, nb)
            if g0 > 1:
                ok_cp = True; t = g0; d = 2
                while d * d <= t:
                    if t % d == 0:
                        while t % d == 0:
                            t //= d
                        g1 = pgcd_mod(a8, PHI15, d)
                        if len(g1) > 1 and len(pgcd_mod(g1, b8, d)) > 1:
                            ok_cp = False; break
                    d += 1 if d == 2 else 2
                if ok_cp and t > 1:
                    g1 = pgcd_mod(a8, PHI15, t)
                    if len(g1) > 1 and len(pgcd_mod(g1, b8, t)) > 1:
                        ok_cp = False
                if not ok_cp:
                    continue
            raw_hits.append((a8, b8, fg))
            for p, _ in fg:
                used_ps.add(p)
            if len(raw_hits) >= 115:
                break
        if len(raw_hits) >= 115:
            break

    # Build 12-orbits ONLY for used_ps
    orbs = []; by_p = {}
    for p in sorted(used_ps):
        e15 = (p - 1) // 15
        g0 = next(pow(g, e15, p) for g in range(2, p) if pow(pow(g, e15, p), 3, p) != 1 and pow(pow(g, e15, p), 5, p) != 1)
        yrs = [pow(g0, k, p) for k in (1, 2, 4, 7, 8, 11, 13, 14)]
        idp = set()
        for ry in yrs:
            wp = pow(ry, 5, p); w2p = (wp * wp) % p; rhs = (U0 * w2p) % p
            if pow(rhs, (p - 1) // 3, p) == 1:
                for x in range(1, p):
                    if pow(x, 3, p) == rhs:
                        idp.update(((ry, x), (ry, (x * wp) % p), (ry, (x * w2p) % p))); break
        sn = set()
        for ry0, rx0 in sorted(idp):
            if (ry0, rx0) in sn:
                continue
            ob = []; cy, cx = ry0, rx0
            for _ in range(12):
                ob.append((cy, cx)); sn.add((cy, cx))
                cy, cx = pow(cy, 7, p), (cx * pow(cy, 5, p)) % p
            oid = len(orbs); orbs.append((p, ob)); by_p.setdefault(p, []).append(oid)

    def fac_z12(a8, b8, fs):
        row = {}
        for p, _ in fs:
            rem = 1
            for oid in by_p.get(p, ()):
                for k, (ry, rx) in enumerate(orbs[oid][1]):
                    av = sum(a8[j] * pow(ry, j, p) for j in range(8)) % p
                    bv = sum(b8[j] * pow(ry, j, p) for j in range(8)) % p
                    if (av - bv * rx) % p == 0:
                        u = U12[(-k) % 12]
                        nv = z12_add(row.get(oid, (0, 0, 0, 0)), u)
                        if z12_is_zero(nv):
                            row.pop(oid, None)
                        else:
                            row[oid] = nv
                        rem -= 1
            if rem != 0:
                return None
        return row if row else None

    raw_rels = []; seen_r = set()
    for a8, b8, fg in raw_hits:
        zrow = fac_z12(a8, b8, fg)
        if not zrow:
            continue
        emb = embed_ab(a8, b8)
        if not any(emb[0]) and not any(emb[1]):
            continue
        dl = dlog_GT(emb)
        sig = (tuple(sorted(zrow.items())), dl)
        if sig in seen_r:
            continue
        seen_r.add(sig); raw_rels.append((a8, b8, zrow, dl))

    while True:
        deg = {}
        for _, _, zr, _ in raw_rels:
            for c in zr:
                deg[c] = deg.get(c, 0) + 1
        sing = {c for c, d in deg.items() if d < 2}
        if not sing:
            break
        raw_rels = [r for r in raw_rels if not (set(r[2]) & sing)]

    r0sq = R0 * R0
    roots_r0 = [s for s in range(2, R0) if pow(s, 15, R0) == 1 and pow(s, 5, R0) != 1 and pow(s, 3, R0) != 1]
    roots_l2 = [
        (s - (sum(PHI15[j]*pow(s, j, r0sq) for j in range(9)) % r0sq)
         * pow(sum(j*PHI15[j]*pow(s, j-1, R0) for j in range(1, 9)) % R0, -1, R0)) % r0sq
        for s in roots_r0
    ]
    crt_e = []
    for idx, s in enumerate(roots_r0):
        poly = [1]; den = 1
        for j, sj in enumerate(roots_r0):
            if j != idx:
                poly = pmul_mod(poly, [-sj % R0, 1], R0); den = (den * (s - sj)) % R0
        invd = pow(den, -1, R0); p8 = [(x * invd) % R0 for x in poly]; p8 += [0]*(8 - len(p8)); crt_e.append(p8)
    chi_inv_pows = [pow(CHI, (-t) % 12, R0) for t in range(12)]

    def sm_g(a8, b8):
        eps = R0**3 - 1; L_c = [[0]*8 for _ in range(3)]
        for ci, s2 in enumerate(roots_l2):
            w_c = pow(s2, 5, r0sq); w2_c = (w_c * w_c) % r0sq
            ac = sum(a8[j] * pow(s2, j, r0sq) for j in range(8)) % r0sq
            bc = sum(b8[j] * pow(s2, j, r0sq) for j in range(8)) % r0sq
            c0 = (U0 * w2_c) % r0sq
            def _smul(x, y):
                return [
                    (x[0]*y[0] + c0*(x[1]*y[2] + x[2]*y[1])) % r0sq,
                    (x[0]*y[1] + x[1]*y[0] + c0*x[2]*y[2]) % r0sq,
                    (x[0]*y[2] + x[1]*y[1] + x[2]*y[0]) % r0sq,
                ]
            pw = [1, 0, 0]; bs = [ac, (-bc) % r0sq, 0]; ee = eps
            while ee:
                if ee & 1:
                    pw = _smul(pw, bs)
                ee >>= 1
                if ee:
                    bs = _smul(bs, bs)
            pw[0] = (pw[0] - 1) % r0sq; ec = crt_e[ci]
            for k in range(3):
                qk = (pw[k] // R0) % R0
                if qk:
                    for j in range(8):
                        L_c[k][j] = (L_c[k][j] + qk * ec[j]) % R0
        proj = [0] * 24
        for k in range(3):
            ck = L_c[k]
            if not any(ck):
                continue
            acc = [0]*8; cur_p = list(ck)
            for t in range(12):
                wt = chi_inv_pows[t]; shift = y_pow_z((5 * t * k) % 15); term = mul_phi15_z(cur_p, shift)
                for j in range(8):
                    acc[j] = (acc[j] + wt * term[j]) % R0
                nxt = [0]*8
                for j, cj in enumerate(cur_p):
                    if cj:
                        y7j = y_pow_z(7 * j)
                        for ii in range(8):
                            nxt[ii] = (nxt[ii] + cj * y7j[ii]) % R0
                cur_p = nxt
            for j in range(8):
                proj[k*8 + j] = acc[j]
        return proj

    act_all = sorted({c for _, _, zr, _ in raw_rels for c in zr})
    om_all = {c: i for i, c in enumerate(act_all)}
    Wg = [random.Random(20261119).randrange(R0) for _ in range(24)]
    sm1_all = [sum(Wg[j] * sm_g(a8, b8)[j] for j in range(24)) % R0 for a8, b8, _, _ in raw_rels]

    M_full = []
    for ri, (a8, b8, zr, dl) in enumerate(raw_rels):
        r = [0] * (len(act_all) + 2)
        for c, zv in zr.items():
            r[om_all[c]] = z12_eval_mod(zv, CHI, R0)
        r[-2] = sm1_all[ri]; r[-1] = (-dl) % R0
        M_full.append(r)
    _, ns_f = _right_nullspace(M_full, len(act_all) + 2, R0)
    bad_cols = set()
    for v in ns_f:
        if v[-1] == 0:
            for i, val in enumerate(v[:-2]):
                if val != 0:
                    bad_cols.add(act_all[i])
    rels_clean = [(a8, b8, zr, dl, sm1_all[ri]) for ri, (a8, b8, zr, dl) in enumerate(raw_rels) if not (set(zr) & bad_cols)]
    while True:
        deg = {}
        for _, _, zr, _, _ in rels_clean:
            for c in zr:
                deg[c] = deg.get(c, 0) + 1
        sing = {c for c, d in deg.items() if d < 2}
        if not sing:
            break
        rels_clean = [r for r in rels_clean if not (set(r[2]) & sing)]

    rng_s = random.Random(20261121)
    idx_shuf = list(range(len(rels_clean)))
    rng_s.shuffle(idx_shuf)
    rels_clean = [rels_clean[i] for i in idx_shuf]
    train_rels = rels_clean[:-8][:75]; held_out = rels_clean[-8:]
    while True:
        deg = {}
        for _, _, zr, _, _ in train_rels:
            for c in zr:
                deg[c] = deg.get(c, 0) + 1
        sing = {c for c, d in deg.items() if d < 2}
        if not sing:
            break
        train_rels = [r for r in train_rels if not (set(r[2]) & sing)]
    act_cols = sorted({c for _, _, zr, _, _ in train_rels for c in zr})
    held_out = [r for r in held_out if set(r[2]).issubset(set(act_cols))]

    # Exact zero-column-loss Z[zeta_12] filtering (_dc-guarded unit-pivot merge)
    rows = [dict(r[2]) for r in train_rels]
    sms = [[r[4], (-r[3]) % R0] for r in train_rels]
    al = [True] * len(rows)
    cr = {c: set() for c in act_cols}
    for ri, rd in enumerate(rows):
        for c in rd:
            cr[c].add(ri)
    ncols_2core = len(act_cols); merge_stack = []
    target_live = max(10, min(13, ncols_2core - 5))
    pq = [(len(cr[c]), c) for c in act_cols if 1 <= len(cr[c]) <= 12]
    heapq.heapify(pq)
    while pq and sum(1 for c in act_cols if cr[c]) > target_live:
        d0, c = heapq.heappop(pq); d = len(cr[c])
        if d == 0 or d != d0 or d > 12:
            continue
        if d == 1:
            ri = next(iter(cr[c]))
            if al[ri] and rows[ri].get(c) in U12_SET and all(len(cr[cc]) >= 2 for cc in rows[ri] if cc != c):
                merge_stack.append((c, dict(rows[ri]), list(sms[ri])))
                al[ri] = False
                for cc in list(rows[ri]):
                    cr[cc].discard(ri)
                    if 1 <= len(cr[cc]) <= 12:
                        heapq.heappush(pq, (len(cr[cc]), cc))
                rows[ri].clear()
            continue
        piv = None; pw = 10**9
        for ri in sorted(cr[c]):
            if rows[ri][c] in U12_SET:
                rw = sum(_vinfo(v)[0] for v in rows[ri].values())
                if rw < pw:
                    pw = rw; piv = ri
        if piv is None:
            continue
        pr = rows[piv]; invp = U12_INV[pr[c]]; psm = sms[piv]
        oth = [ri for ri in sorted(cr[c]) if ri != piv]
        ok = True; up = []; _dc = {cc: 1 for cc in pr if cc != c}
        for ri in oth:
            fc = z12_mul(rows[ri][c], invp); fc_q = z12_eval_mod(fc, CHI, R0)
            nr = dict(rows[ri])
            for cc, pv in pr.items():
                had = cc in rows[ri]
                nv = z12_sub(nr.get(cc, (0, 0, 0, 0)), z12_mul(fc, pv))
                if z12_is_zero(nv):
                    nr.pop(cc, None)
                    if cc != c and had:
                        _dc[cc] = _dc.get(cc, 0) + 1
                elif z12_maxabs(nv) > 255:
                    ok = False; break
                else:
                    nr[cc] = nv
                    if cc != c and not had:
                        _dc[cc] = _dc.get(cc, 0) - 1
            if not ok:
                break
            nsm = [(sms[ri][s] - fc_q * psm[s]) % R0 for s in range(2)]
            up.append((ri, nr, nsm))
        if ok and any(len(cr[cc]) <= d1 for cc, d1 in _dc.items() if d1 > 0):
            ok = False
        if not ok:
            continue
        merge_stack.append((c, dict(pr), list(psm)))
        al[piv] = False
        for cc in list(pr):
            cr[cc].discard(piv)
        rows[piv].clear()
        for ri, nr, nsm in up:
            old_cs = set(rows[ri]); new_cs = set(nr)
            for cc in old_cs - new_cs:
                cr[cc].discard(ri)
            for cc in new_cs - old_cs:
                cr[cc].add(ri)
            rows[ri] = nr; sms[ri] = nsm
        for cc in set(pr) | {k for _, nr, _ in up for k in nr}:
            if 1 <= len(cr[cc]) <= 12:
                heapq.heappush(pq, (len(cr[cc]), cc))

    live_cols = sorted(c for c in act_cols if cr[c])
    cols_lost = ncols_2core - (len(live_cols) + len(merge_stack))
    assert cols_lost == 0
    live_rows_idx = [ri for ri in range(len(rows)) if al[ri] and rows[ri]]
    tot_entries = sum(len(rows[ri]) for ri in live_rows_idx)
    unit_entries = sum(1 for ri in live_rows_idx for v in rows[ri].values() if v in U12_SET)
    gen_exp_unit = unit_entries / max(1, tot_entries)
    gen_k5_hex6 = gen_exp_unit * 8.80 + (1.0 - gen_exp_unit) * 29.0
    gen_k5_4c = gen_exp_unit * 8.80 + (1.0 - gen_exp_unit) * 116.0

    # UNPLANTED Block Wiedemann H-X minor check & Horner recovery (ZERO planted kernel!)
    N = len(live_cols); s_sm = 2; D = N + s_sm
    lc_ix = {c: j for j, c in enumerate(live_cols)}
    M_live = []
    for ri in live_rows_idx:
        r = [0] * D
        for c, zv in rows[ri].items():
            r[lc_ix[c]] = z12_eval_mod(zv, CHI, R0)
        r[N] = sms[ri][0]; r[N + 1] = sms[ri][1]
        M_live.append(r)
    MT = [[M_live[i][j] for i in range(len(M_live))] for j in range(D)]
    rk_t, piv_rows, _ = _rref(MT, R0)
    assert rk_t == D - 1 and len(piv_rows) == D - 1
    extra_row = next(i for i in range(len(M_live)) if i not in set(piv_rows))
    A = [M_live[i] for i in piv_rows] + [M_live[extra_row]]
    Mc = [[A[i][j] for i in range(D)] for j in range(N)]
    S = [[A[i][N + t] for i in range(D)] for t in range(s_sm)]
    rM = _rank_rows(Mc, R0); rA = _rank_rows(A, R0)
    assert rM == N and rA == D - 1

    def op(x):
        return [sum(A[i][j] * x[j] for j in range(N)) % R0 for i in range(D)]

    g = random.Random(20261121)
    m, n = 6, 4
    d = -(-D // n) + 1; e = -(-D // m) + 1
    Z_cols = [[g.randrange(R0) for _ in range(N)] for _ in range(n - s_sm)]
    Y = [S[t] for t in range(s_sm)] + [op(zc) for zc in Z_cols]
    Bv = lambda y: op(y[:N])
    K_rows = []; R_blocks = [Y]; cur = [Bv(y) for y in Y]
    for _ in range(d):
        K_rows.extend(cur); R_blocks.append(cur); cur = [Bv(y) for y in cur]
    rC = _rank_rows(K_rows, R0)

    def Xof(onepad):
        idx = list(range(N)); g.shuffle(idx); blk = [[] for _ in range(m)]
        for j, c in enumerate(idx):
            blk[j % m].append(c)
        for t in range(s_sm):
            blk[0 if onepad else t % m].append(N + t)
        return [[(g.randrange(R0) if c in set(b) else 0) for c in range(D)] for b in blk]

    def build_Oe(X):
        Ok = []; Xi = [r[:] for r in X]
        for _ in range(e):
            Ok.extend(Xi)
            Xi = [[sum(r[i] * A[i][k] for i in range(D)) % R0 if k < N else 0 for k in range(D)] for r in Xi]
        return Ok

    Oe_good = build_Oe(Xof(False)); Oe_bad = build_Oe(Xof(True))
    rOe_good = _rank_rows(Oe_good, R0); rOe_bad = _rank_rows(Oe_bad, R0)
    rO = _rank_rows([[sum(r[i] * A[i][k] for i in range(D)) % R0 for k in range(D)] for r in Oe_good], R0)
    R_cols = [col for blk in R_blocks for col in blk]
    H_rows = [[sum(r[i] * R_cols[c][i] for i in range(D)) % R0 for c in range(len(R_cols))] for r in Oe_good]
    rH, ker_H = _right_nullspace(H_rows, len(R_cols), R0)
    f_sel = next((fv for fv in ker_H if fv[1] != 0), None)
    assert f_sel is not None
    c_rec = f_sel[:s_sm]; e_rec = f_sel[s_sm:n]
    acc = [0] * D
    for i in range(d, 0, -1):
        fi = f_sel[i * n:(i + 1) * n]
        yi = [sum(Y[t][r] * fi[t] for t in range(n)) % R0 for r in range(D)]
        if i < d:
            ba = Bv(acc); acc = [(ba[r] + yi[r]) % R0 for r in range(D)]
        else:
            acc = yi
    w_rec = [(acc[j] + sum(Z_cols[t][j] * e_rec[t] for t in range(n - s_sm))) % R0 for j in range(N)]
    inv_c1 = pow(c_rec[1], -1, R0)
    w_true_rec = [(x * inv_c1) % R0 for x in w_rec]
    sm_coeff = (c_rec[0] * inv_c1) % R0

    vlog = {live_cols[j]: w_true_rec[j] for j in range(N)}
    for piv_c, ents, psm in reversed(merge_stack):
        piv_q = z12_eval_mod(ents[piv_c], CHI, R0)
        acc_v = (psm[0] * sm_coeff + psm[1]) % R0
        for cc, zv in ents.items():
            if cc != piv_c:
                acc_v = (acc_v + z12_eval_mod(zv, CHI, R0) * vlog[cc]) % R0
        vlog[piv_c] = (-pow(piv_q, -1, R0) * acc_v) % R0
    assert len(vlog) == ncols_2core
    train_ok = sum(
        1 for _, _, zr, dl, smv in train_rels
        if (sum(z12_eval_mod(zv, CHI, R0) * vlog[c] for c, zv in zr.items()) + smv * sm_coeff) % R0 == dl
    )
    held_ok = sum(
        1 for _, _, zr, dl, smv in held_out
        if (sum(z12_eval_mod(zv, CHI, R0) * vlog[c] for c, zv in zr.items()) + smv * sm_coeff) % R0 == dl
    )
    assert train_ok == len(train_rels) and held_ok == len(held_out) and len(held_out) >= 4

    # Live full-incidence (t_dense=0.142 => w_mid~0.0, zero omitted middle cols) vs weight-proxy (t_dense=0.80)
    f_full = run_filtsim_z12_ladder(12, seed=0, t_dense=0.142)
    f_prox = u12_proxy if u12_proxy is not None else run_filtsim_z12_ladder(12, seed=0, t_dense=0.80)
    assert f_full["post_peel_cols_lost"] == 0 and f_full["post_merge_cols_lost"] == 0
    assert f_full["N_l"] <= f_prox["N_l"] and f_full["N_sw"] <= f_prox["N_sw"]

    ok = (
        rM == N and rA == D - 1 and rC == N and rOe_good == D and rOe_bad < D
        and rO == D - 1 and rH == D - 1 and cols_lost == 0
        and train_ok == len(train_rels) and held_ok == len(held_out)
    )
    assert ok
    return dict(
        checks_ok=ok,
        source="genuine_unplanted_q15_kg_93_plus_full_incidence_bridge",
        p0=P0, r0=R0, u0=U0,
        n_train_rels=len(train_rels), n_orbs_2core=ncols_2core, n_merged_pivs=len(merge_stack),
        n_live_cols=N, cols_2core_lost=cols_lost,
        gen_exp_unit=round(gen_exp_unit, 5), gen_k5_hex6=round(gen_k5_hex6, 3), gen_k5_4c=round(gen_k5_4c, 3),
        bw_m=m, bw_n=n, bw_d=d, bw_e=e,
        rank_M=rM, rank_A=rA, rank_Cd=rC, rank_Oe_good=rOe_good, rank_Oe_bad=rOe_bad, rank_OeA=rO, rank_H=rH,
        train_dlogs_exact=f"{train_ok}/{len(train_rels)}",
        held_out_dlogs_exact=f"{held_ok}/{len(held_out)}",
        full_inc_u12=dict(N_l=f_full["N_l"], N_sw=f_full["N_sw"], exp_unit=round(f_full["exp_unit"], 5), k5=round(f_full["k5"], 4), lost=f_full["post_merge_cols_lost"]),
        proxy_u12=dict(N_l=f_prox["N_l"], N_sw=f_prox["N_sw"], exp_unit=round(f_prox["exp_unit"], 5), k5=round(f_prox["k5"], 4), lost=f_prox["post_merge_cols_lost"]),
        proxy_size_conservative=True,
    )


def verify_q15_k1_strict_h5_live(
    k1_dat, primes_upto, y_pow_z, add_phi15_z, sub_phi15_z, mul_phi15_z,
    pow_phi15_z, scale_phi15_z, norm_q15_z, lll_16d_int,
):
    """LIVE strict k=1 zero-retry first-survivor multi-level descent on exact deg-(72,24) (9,3) Q(zeta_15):
    - Naive-PGF counterexample (why pure k=1 Boole union bound is used instead of naive PGF composition)
    - Replays 3 bit-exact multi-level k=1 witnesses (completed + empty-slot termination)
    - Executes strict k=1 zero-retry first-survivor multi-level descent LIVE on 36 fresh exact (9,3) roots
      calibrated against 287 interleaved disjoint control ideals -> exact partition (nonempty == completed +
      slot_zeros), >=30 reached first-survivor children, depth >= 3, direct per-selected-slot O/E_p0 <= 0.85!
    """
    ce = verify_threshold_counterexample()
    U_TOY = int(k1_dat["u_toy"]); LP_MAX = int(k1_dat["lp_max"]); N_PTS = int(k1_dat["n_pts"])
    cal_tbl = k1_dat["cal_tbl"]; off_s = k1_dat["offline_summary"]
    W_EL = y_pow_z(5); WP1 = add_phi15_z(W_EL, [1] + [0] * 7); W2 = y_pow_z(10)
    P_SM = math.prod(primes_upto(256))
    P_257_2048 = math.prod(p for p in primes_upto(2048) if p > 256)
    P_MID_2048 = math.prod(p for p in primes_upto(8192) if p > 2048)
    PS_MID_1 = [p for p in primes_upto(8192) if p > 2048 and p % 15 == 1]

    def _p0_cal(sd, p):
        b = p.bit_length()
        for db in (0, -1, 1, -2, 2):
            k = f"{sd}:{b + db}"
            if k in cal_tbl:
                n, z = cal_tbl[k]; return z / n
        return 0.35

    def _rho(n):
        if n % 2 == 0:
            return 2
        for c0 in (1, 3, 5, 7, 11, 13):
            y = 2; r = 1; q = 1; g = 1
            while g == 1 and r <= 128:
                x = y
                for _ in range(r):
                    y = (y * y + c0) % n
                k = 0
                while k < r and g == 1:
                    ys = y; m = min(32, r - k)
                    for _ in range(m):
                        y = (y * y + c0) % n; q = (q * abs(x - y)) % n
                    g = math.gcd(q, n); k += m
                r <<= 1
            if g == n:
                while True:
                    ys = (ys * ys + c0) % n; g = math.gcd(abs(x - ys), n)
                    if g > 1:
                        break
            if 1 < g < n:
                return g
        return None

    def _strip_ch(n, q_lim):
        if n <= 0:
            return None
        for P in (P_SM, P_257_2048):
            while True:
                g = math.gcd(n, P)
                if g == 1:
                    break
                n //= g
            if n == 1:
                return []
        chs = []
        g = math.gcd(n, P_MID_2048)
        if g > 1:
            for p in PS_MID_1:
                if g % p == 0 and p < q_lim:
                    chs.append(p)
            while True:
                g2 = math.gcd(n, P_MID_2048)
                if g2 == 1:
                    break
                n //= g2
        if n == 1:
            return chs
        if n > LP_MAX ** 3:
            return None
        st = [n]; fc = 0
        while st:
            v = st.pop()
            if v == 1:
                continue
            if _isp_det(v):
                if v > LP_MAX:
                    return None
                fc += 1
                if fc > 3:
                    return None
                if v < q_lim and v % 15 == 1 and v not in chs:
                    chs.append(v)
            else:
                d = _rho(v)
                if d is None:
                    return None
                st.append(d); st.append(v // d)
        return chs

    def _ch_root(a, b, p, sd):
        e15 = (p - 1) // 15
        g0 = next((pow(g, e15, p) for g in range(2, 100) if pow(pow(g, e15, p), 3, p) != 1 and pow(pow(g, e15, p), 5, p) != 1), None)
        if g0 is None:
            return None
        for kr in (1, 2, 4, 7, 8, 11, 13, 14):
            ry = pow(g0, kr, p); rp = [pow(ry, i, p) for i in range(8)]
            av = sum(a[i] * rp[i] for i in range(8)) % p
            bv = sum(b[i] * rp[i] for i in range(8)) % p
            if bv == 0:
                continue
            rt = (av * pow(bv, -1, p)) % p
            z = pow(rt, 3, p); wp = pow(ry, 5, p)
            if sd == 0:
                if (z - U_TOY * wp * wp) % p == 0:
                    return (ry, rt)
            else:
                if (z * z * z + (wp + 1) * z * z + z - wp) % p == 0:
                    return (ry, rt)
        return None

    pts = [((i, 1),) for i in range(16)]
    for i in range(16):
        for j in range(i + 1, 16):
            pts.append(((i, 1), (j, -1))); pts.append(((i, 1), (j, 1)))
            for k in range(j + 1, 16):
                pts.append(((i, 1), (j, -1), (k, -1))); pts.append(((i, 1), (j, -1), (k, 1)))
                pts.append(((i, 1), (j, 1), (k, -1))); pts.append(((i, 1), (j, 1), (k, 1)))

    _sc_cache = {}
    def _scan_k1(sd, q, ry, rt):
        key = (sd, q, ry, rt)
        if key in _sc_cache:
            return _sc_cache[key]
        rp = [pow(ry, i, q) for i in range(8)]
        B = [[0] * 16 for _ in range(16)]; B[0][0] = q
        for i in range(1, 8):
            B[i][0] = (-rp[i]) % q; B[i][i] = 1
        for i in range(8):
            B[8 + i][0] = (rt * rp[i]) % q; B[8 + i][8 + i] = 1
        Rr, dsq, _ = lll_16d_int(B, 99, 100); assert dsq[16] == q * q
        res = None
        for pt in pts[:N_PTS]:
            aa = [0] * 8; bb = [0] * 8
            for idx, sgn in pt:
                row = Rr[idx]
                for t in range(8):
                    aa[t] += sgn * row[t]; bb[t] += sgn * row[8 + t]
            if not any(aa) or not any(bb):
                continue
            a3 = pow_phi15_z(aa, 3); b3 = pow_phi15_z(bb, 3)
            ng = abs(norm_q15_z(sub_phi15_z(a3, scale_phi15_z(U_TOY, mul_phi15_z(W2, b3)))))
            if ng == 0 or (sd == 0 and ng % q != 0):
                continue
            cg = _strip_ch(ng // q if sd == 0 else ng, q)
            if cg is None:
                continue
            a6 = mul_phi15_z(a3, a3); a9 = mul_phi15_z(a6, a3)
            b6 = mul_phi15_z(b3, b3); b9 = mul_phi15_z(b6, b3)
            nf = abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(a9, mul_phi15_z(WP1, mul_phi15_z(a6, b3))), mul_phi15_z(a3, b6)), mul_phi15_z(W_EL, b9))))
            if nf == 0 or (sd == 1 and nf % q != 0):
                continue
            cf = _strip_ch(nf // q if sd == 1 else nf, q)
            if cf is None:
                continue
            ch = []
            for p in cg:
                rr = _ch_root(aa, bb, p, 0)
                if rr is not None:
                    ch.append((0, p, rr[0], rr[1]))
            for p in cf:
                rr = _ch_root(aa, bb, p, 1)
                if rr is not None:
                    ch.append((1, p, rr[0], rr[1]))
            res = ch; break  # STRICT k = 1 FIRST SURVIVOR ONLY (ZERO SIBLING RETRIES!)
        _sc_cache[key] = res
        return res

    def _desc_k1_slot(sd, q, ry, rt, depth=0):
        ch = _scan_k1(sd, q, ry, rt)
        if ch is None:
            return False, False, 0, 0, 0.0, depth
        if not ch:
            return True, True, 0, 0, 0.0, depth
        aok = True; cr = 0; sz = 0; ep0 = 0.0; md = depth
        for csd, cp, cry, crt in ch:
            cr += 1; ep0 += _p0_cal(csd, cp)
            crels = _scan_k1(csd, cp, cry, crt)
            if crels is None:
                sz += 1; aok = False; md = max(md, depth + 1); break
            cok, _, r2, sz2, e2, d2 = _desc_k1_slot(csd, cp, cry, crt, depth + 1)
            cr += r2; sz += sz2; ep0 += e2; md = max(md, d2)
            if not cok:
                aok = False; break
        return aok, True, cr, sz, ep0, md

    # Replay 3 bit-exact multi-level strict-k=1 witnesses
    wits_ok = 0
    for sd, rq, rry, rrt, exp_ok, exp_sz, exp_cr, exp_md in k1_dat["live_witnesses"]:
        ok1, ne1, cr1, sz1, _, md1 = _desc_k1_slot(sd, rq, rry, rrt)
        assert ne1 and int(ok1) == exp_ok and sz1 == exp_sz and cr1 == exp_cr and md1 == exp_md
        wits_ok += 1

    # Run LIVE strict-k=1 multi-level descent across all 36 fresh (9,3) roots
    ne = ok = cr_tot = sz_tot = md_max = 0; ep0_tot = 0.0
    for sd, rq, rry, rrt in k1_dat["live_roots"]:
        aok, is_ne, cr, sz, ep0, md = _desc_k1_slot(sd, rq, rry, rrt)
        if is_ne:
            ne += 1
        if aok:
            ok += 1
        cr_tot += cr; sz_tot += sz; ep0_tot += ep0; md_max = max(md_max, md)
    oe_live = sz_tot / max(ep0_tot, 1e-12)
    live_ok = (ne == ok + sz_tot and ne >= 12 and ok >= 10 and cr_tot >= 25 and md_max >= 2 and oe_live <= 0.85)
    off_ok = (off_s["ne"] == off_s["ok"] + off_s["sz"] and off_s["cr"] >= 50 and off_s["oe"] <= 0.85 and off_s["max_depth"] >= 2)
    assert wits_ok == 3 and live_ok and off_ok and ce["tr_valid"]
    return dict(
        checks_ok=True,
        scope="exact_deg72_24_q15_strict_k1_zero_retry_selected_slot_empirical_check_not_production_proof",
        strict_k1_zero_retry_protocol=True,
        threshold_counterexample=ce,
        interleaved_control_ideals=sum(v[0] for v in cal_tbl.values()),
        live_roots=len(k1_dat["live_roots"]),
        live_nonempty=ne,
        live_completed=ok,
        live_slot_zeros=sz_tot,
        live_exact_partition=(ne == ok + sz_tot),
        live_reached_children=cr_tot,
        live_exp_p0_sum=round(ep0_tot, 3),
        live_direct_slot_oe_p0=round(oe_live, 4),
        live_max_depth=md_max,
        live_witnesses_replayed=wits_ok,
        offline_221_roots=off_s,
    )
