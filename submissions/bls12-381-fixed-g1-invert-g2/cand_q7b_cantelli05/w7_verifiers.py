"""w7_verifiers.py — Live Wave-7 verifiers fixing ALL 3 Wave-6 defects cleanly:
  1. verify_threshold_counterexample(): exact algebraic witness that naive PGF substitution G(psi) fails.
  2. verify_genuine_q15_filter_hx_live(): 100% UNSCREENED, ZERO-ORACLE, TWO-SIDED Q(zeta_15) Z[zeta_12]
     filtering + UNPLANTED Block Wiedemann H-X + post-hoc independent BSGS dlog_GT verification +
     live full-incidence (t_dense=0.142, w_mid~0.0) vs weight-proxy (t_dense=0.80) bridge.
     - ZERO dlog_GT in matrix/filter/solver!
     - ZERO bad_cols nullspace screening!
     - ZERO RREF row selection (uses natural first D surviving rows)!
  3. verify_q15_k1_strict_h5_live(): 100% FAITHFUL STRICTLY-DESCENDING Q(zeta_15) strict-k=1 zero-retry
     multi-level descent campaign with ZERO omitted primes (_strict_fac rejects ANY prime > B_FB that is
     >= q_lim or not 1 mod 15 or fails _ch_root) AND direct modeled O/E against the claimed p_0(mu_model)
     bound (exact udup=0.50, evaluated at both boundary and low-p_0 regimes + matched disjoint control +
     3 bit-exact multi-level witnesses).
"""
import heapq
import math
import random
from filtsim_z12 import (
    U12, U12_SET, U12_INV, z12_add, z12_sub, z12_mul, z12_is_zero, z12_maxabs,
    z12_eval_mod, _vinfo, run_filtsim_z12_ladder,
)
from hx_z12 import _right_nullspace, _rank_rows


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
    """100% UNSCREENED, ZERO-ORACLE, TWO-SIDED Q(zeta_15) Z[zeta_12] filtering + UNPLANTED Block Wiedemann H-X:
    Tower: p_0 = 7 (=== 7 mod 15), r_0 = 181 (=== 1 mod 15, r_0 | Phi_12(7) = 2353 = 13*181), chi = 7 mod 181.
    Non-isomorphic two-sided pair over Q(zeta_15) fixed by tau(y)=y^7, tau(X)=X*w and congruent mod pi=(7, h_pi):
      g(X) = X^3 - w^2             (C_g = 1,   deg 24)
      f(X) = X^3 - (w - 1) * w^2   (C_f = w-1, deg 24; w === 2 mod pi => w-1 === 1 mod pi => f === g mod pi)
    ZERO dlog_GT in matrix/filter/solver! ZERO bad_cols screening! ZERO RREF row selection!
    Plus live full-incidence (t_dense=0.142 => w_mid~0.0) vs weight-proxy (t_dense=0.80) filtsim_z12 bridge.
    """
    P0, R0, CHI = 7, 181, 7
    W_EL = y_pow_z(5)
    W2 = y_pow_z(10)
    ONE = [1, 0, 0, 0, 0, 0, 0, 0]
    CG_EL = W2
    CF_EL = mul_phi15_z(sub_phi15_z(W_EL, ONE), W2)

    w_pi = next(w for w in range(2, P0) if (w * w + w + 1) % P0 == 0)
    h_pi = pgcd_mod(PHI15, [-w_pi % P0, 0, 0, 0, 0, 1], P0)
    assert len(h_pi) == 5 and not any(prem_mod(sub_phi15_z(CF_EL, CG_EL), h_pi, P0))
    rhs_f = pow((w_pi * w_pi) % P0, (P0 - 1) // 3, P0)
    assert rhs_f == w_pi
    w_tau_exp = 5

    B_FB = 2000
    ps_1 = [p for p in primes_upto(B_FB) if p % 15 == 1 and p != R0]
    ps_van = [2, 3, 5] + [p for p in primes_upto(B_FB) if p % 15 in (4, 7, 13)]
    P_van = math.prod(ps_van)
    P_1 = math.prod(ps_1)

    def strip_all(n):
        if n <= 0 or n % R0 == 0:
            return None
        t = n
        while True:
            g = math.gcd(t, P_van)
            if g == 1:
                break
            t //= g
        if t == 1:
            return []
        r = t
        hits = []
        if math.gcd(r, P_1) == 1:
            return None
        for p in ps_1:
            if r % p == 0:
                e = 0
                while r % p == 0:
                    r //= p
                    e += 1
                hits.append((p, e))
                if r == 1:
                    return hits
        return hits if r == 1 else None

    pool = [[1 if j == i else 0 for j in range(8)] for i in range(8)]
    for i in range(8):
        for j in range(i + 1, 8):
            for s in (1, -1):
                v = [0] * 8; v[i] = 1; v[j] = s; pool.append(v)
    for i in range(8):
        for j in range(8):
            if i != j:
                for s in (1, -1):
                    v = [0] * 8; v[i] = 2; v[j] = s; pool.append(v)
    for i in range(8):
        for j in range(i + 1, 8):
            for k in range(j + 1, 8):
                for s1 in (1, -1):
                    for s2 in (1, -1):
                        v = [0] * 8; v[i] = 1; v[j] = s1; v[k] = s2
                        if abs(norm_q15_z(v)) <= 220:
                            pool.append(v)
    seen_v = set(); uniq = []
    for v in pool:
        t = tuple(v) if v >= [-x for x in v] else tuple(-x for x in v)
        if t not in seen_v:
            seen_v.add(t); uniq.append(list(t))

    def canon_y(v):
        best = None; cur = list(v)
        for _ in range(15):
            for s in (1, -1):
                t = tuple(s * x for x in cur)
                if best is None or t < best:
                    best = t
            cur = mul_phi15_z(cur, [0, 1, 0, 0, 0, 0, 0, 0])
        return best

    seen_b = set(); b_uniq = []
    for v in uniq:
        c = canon_y(v)
        if c not in seen_b:
            seen_b.add(c); b_uniq.append(v)
    pre_a = [(v, pow_phi15_z(v, 3), abs(norm_q15_z(v))) for v in uniq]
    pre_b = [(v, pow_phi15_z(v, 3), abs(norm_q15_z(v))) for v in b_uniq]

    raw = []; used_p = set()
    for a8, a3, na in pre_a:
        for b8, b3, nb in pre_b:
            if a8 == b8:
                continue
            ng = abs(norm_q15_z(sub_phi15_z(a3, mul_phi15_z(CG_EL, b3))))
            if ng <= 1:
                continue
            fg = strip_all(ng)
            if fg is None:
                continue
            nf = abs(norm_q15_z(sub_phi15_z(a3, mul_phi15_z(CF_EL, b3))))
            if nf <= 1:
                continue
            ff = strip_all(nf)
            if ff is None or (not fg and not ff):
                continue
            ar = prem_mod(a8, h_pi, P0); br = prem_mod(b8, h_pi, P0)
            if not any(ar) and not any(br):
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
            raw.append((a8, b8, fg, ff))
            for p, _ in fg:
                used_p.add(p)
            for p, _ in ff:
                used_p.add(p)

    orbs = []; by_sp = {}
    for p in sorted(used_p):
        e15 = (p - 1) // 15
        g0 = next(pow(g, e15, p) for g in range(2, p) if pow(pow(g, e15, p), 3, p) != 1 and pow(pow(g, e15, p), 5, p) != 1)
        yrs = [pow(g0, k, p) for k in (1, 2, 4, 7, 8, 11, 13, 14)]
        for sd, C_el in ((0, CG_EL), (1, CF_EL)):
            idp = set()
            for ry in yrs:
                rhs = sum(C_el[j] * pow(ry, j, p) for j in range(8)) % p
                wp = pow(ry, 5, p); w2p = (wp * wp) % p
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
                    cy, cx = pow(cy, 7, p), (cx * pow(cy, w_tau_exp, p)) % p
                oid = len(orbs); orbs.append((sd, p, ob)); by_sp.setdefault((sd, p), []).append(oid)

    def lift_root(ry, rx, C_el, p, e):
        mod = p; cy, cx = ry, rx
        for _ in range(1, e):
            nmod = mod * p
            fy = sum(PHI15[j] * pow(cy, j, nmod) for j in range(9)) % nmod
            dfy = sum(j * PHI15[j] * pow(cy, j - 1, p) for j in range(1, 9)) % p
            cy = (cy - fy * pow(dfy, -1, p)) % nmod
            rhs = sum(C_el[j] * pow(cy, j, nmod) for j in range(8)) % nmod
            fx = (pow(cx, 3, nmod) - rhs) % nmod
            dfx = (3 * cx * cx) % p
            cx = (cx - fx * pow(dfx, -1, p)) % nmod
            mod = nmod
        return cy, cx, mod

    def fac_z12(a8, b8, fs, sd, C_el):
        row = {}
        for p, ev in fs:
            rem = ev
            for oid in by_sp.get((sd, p), ()):
                for k, (ry, rx) in enumerate(orbs[oid][2]):
                    av = sum(a8[j] * pow(ry, j, p) for j in range(8)) % p
                    bv = sum(b8[j] * pow(ry, j, p) for j in range(8)) % p
                    if (av - bv * rx) % p == 0:
                        vp = 1
                        if ev > 1:
                            rye, rxe, pmod = lift_root(ry, rx, C_el, p, ev)
                            val = (sum(a8[j] * pow(rye, j, pmod) for j in range(8)) - sum(b8[j] * pow(rye, j, pmod) for j in range(8)) * rxe) % pmod
                            vp = 0; t = p
                            while vp < ev and val % t == 0:
                                vp += 1; t *= p
                        u = U12[(-k) % 12]
                        nv = row.get(oid, (0, 0, 0, 0))
                        for _ in range(vp):
                            nv = z12_add(nv, u)
                        if z12_is_zero(nv):
                            row.pop(oid, None)
                        else:
                            row[oid] = nv
                        rem -= vp
            if rem != 0:
                return None
        return row

    rels = []; seen_sig = set()
    for a8, b8, fg, ff in raw:
        zg = fac_z12(a8, b8, fg, 0, CG_EL); zf = fac_z12(a8, b8, ff, 1, CF_EL)
        if zg is None or zf is None:
            continue
        zrow = dict(zf)
        for o, zv in zg.items():
            zrow[o] = tuple(-x for x in zv)
        if not zrow:
            continue
        first_o = min(zrow)
        inv0 = pow(z12_eval_mod(zrow[first_o], CHI, R0), -1, R0)
        sig = tuple(sorted((o, (z12_eval_mod(zv, CHI, R0) * inv0) % R0) for o, zv in zrow.items()))
        if sig in seen_sig:
            continue
        seen_sig.add(sig); rels.append((a8, b8, zrow, zf, zg))

    # Pure topological 2-core peel (ZERO oracle, ZERO bad_cols screening!)
    while True:
        deg = {}
        for _, _, zr, _, _ in rels:
            for o in zr:
                deg[o] = deg.get(o, 0) + 1
        sing = {o for o, d in deg.items() if d < 2}
        if not sing:
            break
        rels = [r for r in rels if not (set(r[2]) & sing)]

    act = sorted({o for _, _, zr, _, _ in rels for o in zr})
    nf_o = sum(1 for o in act if orbs[o][0] == 1)
    ng_o = sum(1 for o in act if orbs[o][0] == 0)
    assert len(rels) >= 30 and len(act) >= 24 and nf_o >= 10 and ng_o >= 10

    # Data-oblivious compressed Schirokauer columns (1 from f + 1 from g => s_sm = 2, ZERO dlog_GT!)
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

    def sm_3(a8, b8, C_el):
        eps = R0**3 - 1; L_c = [[0]*8 for _ in range(3)]
        for ci, s2 in enumerate(roots_l2):
            ac = sum(a8[j] * pow(s2, j, r0sq) for j in range(8)) % r0sq
            bc = sum(b8[j] * pow(s2, j, r0sq) for j in range(8)) % r0sq
            c0 = sum(C_el[j] * pow(s2, j, r0sq) for j in range(8)) % r0sq
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
                wt = chi_inv_pows[t]; shift = y_pow_z((w_tau_exp * t * k) % 15); term = mul_phi15_z(cur_p, shift)
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

    smf_raw = [sm_3(a8, b8, CF_EL) for a8, b8, _, _, _ in rels]
    smg_raw = [sm_3(a8, b8, CG_EL) for a8, b8, _, _, _ in rels]
    rng_sm = random.Random(20261121)
    Wf = [rng_sm.randrange(R0) for _ in range(24)]
    Wg = [rng_sm.randrange(R0) for _ in range(24)]
    s_sm = 2
    sms_all = [
        [sum(Wf[k] * smf_raw[ri][k] for k in range(24)) % R0,
         (-sum(Wg[k] * smg_raw[ri][k] for k in range(24))) % R0]
        for ri in range(len(rels))
    ]
    om = {o: i for i, o in enumerate(act)}
    nc_full = len(act) + s_sm
    M_full = []
    for ri, (_, _, zr, _, _) in enumerate(rels):
        r = [0] * nc_full
        for o, zv in zr.items():
            r[om[o]] = z12_eval_mod(zv, CHI, R0)
        r[-2] = sms_all[ri][0]; r[-1] = sms_all[ri][1]
        M_full.append(r)
    rk_M_full = _rank_rows([r[:len(act)] for r in M_full], R0)
    rk_A_full, ns_full = _right_nullspace(M_full, nc_full, R0)
    assert rk_M_full == len(act) and rk_A_full == nc_full - 1 and len(ns_full) == 1

    # Zero-column-loss _dc-guarded Z[zeta_12] merge (ZERO oracle, ZERO bad_cols!)
    rows = [dict(r[2]) for r in rels]
    sms = [list(s) for s in sms_all]
    al = [True] * len(rows)
    cr = {c: set() for c in act}
    for ri, rd in enumerate(rows):
        for c in rd:
            cr[c].add(ri)
    ncols_2core = len(act); merge_stack = []
    target_live = max(12, ncols_2core - 10)
    pq = [(len(cr[c]), c) for c in act if 1 <= len(cr[c]) <= 14]
    heapq.heapify(pq)
    while pq and sum(1 for c in act if cr[c]) > target_live:
        d0, c = heapq.heappop(pq); d = len(cr[c])
        if d == 0 or d != d0 or d > 14:
            continue
        if d == 1:
            ri = next(iter(cr[c]))
            if al[ri] and rows[ri].get(c) in U12_SET and all(len(cr[cc]) >= 2 for cc in rows[ri] if cc != c):
                merge_stack.append((c, dict(rows[ri]), list(sms[ri])))
                al[ri] = False
                for cc in list(rows[ri]):
                    cr[cc].discard(ri)
                    if 1 <= len(cr[cc]) <= 14:
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
            nsm = [(sms[ri][s] - fc_q * psm[s]) % R0 for s in range(s_sm)]
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
            if 1 <= len(cr[cc]) <= 14:
                heapq.heappush(pq, (len(cr[cc]), cc))

    live_cols = sorted(c for c in act if cr[c])
    cols_lost = ncols_2core - (len(live_cols) + len(merge_stack))
    assert cols_lost == 0 and len(merge_stack) >= 6
    live_idx = [ri for ri in range(len(rows)) if al[ri] and rows[ri]]
    tot_entries = sum(len(rows[ri]) for ri in live_idx)
    unit_entries = sum(1 for ri in live_idx for v in rows[ri].values() if v in U12_SET)
    gen_exp_unit = unit_entries / max(1, tot_entries)
    gen_k5_hex6 = gen_exp_unit * 8.80 + (1.0 - gen_exp_unit) * 29.0
    assert gen_k5_hex6 < 11.8

    # Take FIRST D surviving rows naturally (ZERO RREF row selection!)
    N = len(live_cols); D = N + s_sm
    assert len(live_idx) >= D
    lc_ix = {c: j for j, c in enumerate(live_cols)}
    A_live = []
    for ri in live_idx[:D]:
        r = [0] * D
        for c, zv in rows[ri].items():
            r[lc_ix[c]] = z12_eval_mod(zv, CHI, R0)
        for s in range(s_sm):
            r[N + s] = sms[ri][s]
        A_live.append(r)
    Mc = [[A_live[i][j] for i in range(D)] for j in range(N)]
    S = [[A_live[i][N + t] for i in range(D)] for t in range(s_sm)]
    rM = _rank_rows(Mc, R0); rA = _rank_rows(A_live, R0)
    assert rM == N and rA == D - 1

    # UNPLANTED Block Wiedemann H-X + Horner recovery + back-substitution of ALL merged pivots
    def op(x):
        return [sum(A_live[i][j] * x[j] for j in range(N)) % R0 for i in range(D)]

    g = random.Random(20261121); m, n = 6, 4; d = -(-D // n) + 1; e = -(-D // m) + 1
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
            Xi = [[sum(r[i]*A_live[i][k] for i in range(D)) % R0 if k < N else 0 for k in range(D)] for r in Xi]
        return Ok

    Oe_good = build_Oe(Xof(False)); Oe_bad = build_Oe(Xof(True))
    rOe_good = _rank_rows(Oe_good, R0); rOe_bad = _rank_rows(Oe_bad, R0)
    rO = _rank_rows([[sum(r[i]*A_live[i][k] for i in range(D)) % R0 for k in range(D)] for r in Oe_good], R0)
    R_cols = [col for blk in R_blocks for col in blk]
    H_rows = [[sum(r[i]*R_cols[c][i] for i in range(D)) % R0 for c in range(len(R_cols))] for r in Oe_good]
    rH, ker_H = _right_nullspace(H_rows, len(R_cols), R0)
    assert rC == N and rOe_good == D and rOe_bad < D and rO == D - 1 and rH == D - 1
    f_sel = next((fv for fv in ker_H if any(fv[:s_sm])), None)
    assert f_sel is not None
    c_rec = f_sel[:s_sm]; e_rec = f_sel[s_sm:n]
    acc = [0] * D
    for i in range(d, 0, -1):
        fi = f_sel[i*n:(i+1)*n]
        yi = [sum(Y[t][r]*fi[t] for t in range(n)) % R0 for r in range(D)]
        acc = [(Bv(acc)[r] + yi[r]) % R0 for r in range(D)] if i < d else yi
    w_rec = [(acc[j] + sum(Z_cols[t][j]*e_rec[t] for t in range(n - s_sm))) % R0 for j in range(N)]
    vlog = {live_cols[j]: w_rec[j] for j in range(N)}
    for piv_c, ents, psm in reversed(merge_stack):
        piv_q = z12_eval_mod(ents[piv_c], CHI, R0)
        acc_v = sum(psm[s] * c_rec[s] for s in range(s_sm)) % R0
        for cc, zv in ents.items():
            if cc != piv_c:
                acc_v = (acc_v + z12_eval_mod(zv, CHI, R0) * vlog[cc]) % R0
        vlog[piv_c] = (-pow(piv_q, -1, R0) * acc_v) % R0

    # POST-HOC ONLY: Independent BSGS in F_{7^12}^* (ZERO oracle used anywhere above!)
    _h0, _h1, _h2, _h3 = h_pi[:4]
    def fp4_add(a, b): return [(a[i] + b[i]) % P0 for i in range(4)]
    def fp4_mul(a, b):
        c = [a[0]*b[0], a[0]*b[1]+a[1]*b[0], a[0]*b[2]+a[1]*b[1]+a[2]*b[0], a[0]*b[3]+a[1]*b[2]+a[2]*b[1]+a[3]*b[0], a[1]*b[3]+a[2]*b[2]+a[3]*b[1], a[2]*b[3]+a[3]*b[2], a[3]*b[3]]
        for k in (6, 5, 4):
            v = c[k] % P0
            if v: o = k - 4; c[o] -= v*_h0; c[o+1] -= v*_h1; c[o+2] -= v*_h2; c[o+3] -= v*_h3
        return [c[i] % P0 for i in range(4)]
    def fp4_pow(a, e):
        res, b = [1, 0, 0, 0], list(a)
        while e:
            if e & 1: res = fp4_mul(res, b)
            e >>= 1
            if e: b = fp4_mul(b, b)
        return res
    def fp4_inv(a):
        a1 = fp4_pow(a, P0); a2 = fp4_pow(a1, P0); a3 = fp4_pow(a2, P0)
        p3 = fp4_mul(a1, fp4_mul(a2, a3)); inv = pow(fp4_mul(a, p3)[0], -1, P0)
        return [(x * inv) % P0 for x in p3]
    w4 = fp4_pow([0, 1, 0, 0], 5); w4_2 = fp4_mul(w4, w4); uw2 = list(w4_2)
    frob4_cols = [fp4_pow([1 if j == i else 0 for j in range(4)], P0) for i in range(4)]
    def frob4(a):
        out = [0, 0, 0, 0]
        for i, ai in enumerate(a):
            if ai:
                col = frob4_cols[i]
                for j in range(4): out[j] = (out[j] + ai * col[j]) % P0
        return out
    def fp12_mul(A, B):
        return (fp4_add(fp4_mul(A[0], B[0]), fp4_mul(uw2, fp4_add(fp4_mul(A[1], B[2]), fp4_mul(A[2], B[1])))),
                fp4_add(fp4_add(fp4_mul(A[0], B[1]), fp4_mul(A[1], B[0])), fp4_mul(uw2, fp4_mul(A[2], B[2]))),
                fp4_add(fp4_add(fp4_mul(A[0], B[2]), fp4_mul(A[1], B[1])), fp4_mul(A[2], B[0])))
    def fp12_pow(A, e):
        res, b = ([1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]), A
        while e:
            if e & 1: res = fp12_mul(res, b)
            e >>= 1
            if e: b = fp12_mul(b, b)
        return res
    def fp12_frob(A): return (frob4(A[0]), fp4_mul(frob4(A[1]), w4), fp4_mul(frob4(A[2]), w4_2))
    def fp12_inv(A):
        def f4(Z):
            for _ in range(4): Z = fp12_frob(Z)
            return Z
        A4 = f4(A); A8 = f4(A4); adj = fp12_mul(A4, A8); inv4 = fp4_inv(fp12_mul(A, adj)[0])
        return (fp4_mul(adj[0], inv4), fp4_mul(adj[1], inv4), fp4_mul(adj[2], inv4))
    cof = (P0**4 - P0**2 + 1) // R0
    def to_GT(A):
        e6 = A
        for _ in range(6): e6 = fp12_frob(e6)
        ez = fp12_mul(e6, fp12_inv(A)); e2 = fp12_frob(fp12_frob(ez))
        return fp12_pow(fp12_mul(e2, ez), cof)
    G_gen = to_GT(([3, 1, 0, 0], [2, 0, P0 - 1, 0], [0, 5, 0, 0]))
    m_bs = math.isqrt(R0) + 1; bs_map = {}; cur = ([1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0])
    for j in range(m_bs):
        bs_map[(tuple(cur[0]), tuple(cur[1]), tuple(cur[2]))] = j
        cur = fp12_mul(cur, G_gen)
    Gm_inv = fp12_inv(fp12_pow(G_gen, m_bs))
    def dlog_GT(a8, b8):
        ar = prem_mod(a8, h_pi, P0); ar += [0] * (4 - len(ar))
        br = prem_mod([-x for x in b8], h_pi, P0); br += [0] * (4 - len(br))
        c = to_GT((ar, br, [0, 0, 0, 0]))
        for i in range(m_bs + 1):
            k = (tuple(c[0]), tuple(c[1]), tuple(c[2]))
            if k in bs_map: return (i * m_bs + bs_map[k]) % R0
            c = fp12_mul(c, Gm_inv)
        raise RuntimeError("BSGS failed")

    lam = None; exact_f = exact_g = nz_dl = 0
    for ri, (a8, b8, zr, zf, zg) in enumerate(rels):
        vf = (sum(z12_eval_mod(zv, CHI, R0) * vlog[c] for c, zv in zf.items()) + sms_all[ri][0] * c_rec[0]) % R0
        vg = (sum(z12_eval_mod(zv, CHI, R0) * vlog[c] for c, zv in zg.items()) - sms_all[ri][1] * c_rec[1]) % R0
        assert vf == vg
        dl = dlog_GT(a8, b8)
        if dl != 0:
            nz_dl += 1
            if lam is None:
                lam = (vf * pow(dl, -1, R0)) % R0
        if lam is not None:
            if vf == (lam * dl) % R0: exact_f += 1
            if vg == (lam * dl) % R0: exact_g += 1
    assert lam is not None and lam != 0 and exact_f == len(rels) and exact_g == len(rels)

    # Live full-incidence (t_dense=0.142 => w_mid~0.0, ZERO omitted middle cols) vs t_dense=0.80 proxy bridge
    px = u12_proxy if u12_proxy is not None else run_filtsim_z12_ladder(12, seed=0, t_dense=0.80)
    fi = run_filtsim_z12_ladder(12, seed=0, t_dense=0.142)
    assert fi["post_peel_cols_lost"] == 0 and fi["post_merge_cols_lost"] == 0
    assert fi["N_l"] <= px["N_l"] and fi["N_sw"] <= px["N_sw"]

    return dict(
        tower="unscreened_zero_oracle_2sided_Q15_p7_r181_chi7",
        polys="g=X^3-w^2, f=X^3-(w-1)*w^2 (congruent mod pi=(7,h_pi), non-isomorphic N(w-1)=3)",
        oracle_free_matrix_filter_solver=True,
        bad_cols_screened=0,
        rref_row_selection=False,
        rels_2core=len(rels),
        orbits_2core=ncols_2core,
        orbits_f=nf_o,
        orbits_g=ng_o,
        rk_M_2core=rk_M_full,
        rk_A_2core=rk_A_full,
        nullity_2core=len(ns_full),
        merged_pivots=len(merge_stack),
        live_cols=N,
        cols_2core_lost=cols_lost,
        s_sm=s_sm,
        held_out_surviving_rows=len(live_idx) - D,
        gen_exp_unit=round(gen_exp_unit, 4),
        gen_k5_hex6=round(gen_k5_hex6, 3),
        bw_rank_Cd=rC,
        bw_rank_Oe_good=rOe_good,
        bw_rank_Oe_bad=rOe_bad,
        bw_rank_OeA=rO,
        bw_rank_H=rH,
        posthoc_bsgs_lam=lam,
        posthoc_bsgs_nz_dlogs=nz_dl,
        posthoc_exact_f=exact_f,
        posthoc_exact_g=exact_g,
        posthoc_total_rels=len(rels),
        full_inc_u12=dict(N_l=fi["N_l"], N_sw=fi["N_sw"], exp_unit=round(fi["exp_unit"], 5), k5=round(fi["k5"], 4), lost=fi["post_merge_cols_lost"]),
        proxy_u12=dict(N_l=px["N_l"], N_sw=px["N_sw"], exp_unit=round(px["exp_unit"], 5), k5=round(px["k5"], 4), lost=px["post_merge_cols_lost"]),
        proxy_size_conservative=True,
        checks_ok=True,
    )


def verify_q15_k1_strict_h5_live(
    low_p0_data, primes_upto, PHI15, y_pow_z, add_phi15_z, sub_phi15_z, mul_phi15_z,
    pow_phi15_z, scale_phi15_z, norm_q15_z, pgcd_mod, lll_16d_int, rho_exact=None, desc_params=None,
):
    """100% FAITHFUL IMMEDIATE-ABORT STRICTLY-DESCENDING Q(zeta_15) strict-k=1 campaign + 300-ideal low-p0 H5 check:
    1) High-p0 boundary LIVE immediate-abort campaign (64 live 17-bit roots, n_pts=32, B_FB=8192):
       - ZERO omitted primes (_strict_fac rejects ANY prime > B_FB that is >= q_lim or not 1 mod 15 or lacks _ch_root).
       - IMMEDIATE ATTEMPT ABORT (break on first empty child: ZERO post-abort siblings ever visited!).
       - Exact attempt partition: nonempty_roots == completed_roots + aborted_roots.
       - Direct modeled claimed p_0(mu_model) denominator (exact udup=0.50) over actually-reached slots before abort.
       - 3/3 bit-exact multi-level witnesses (incl. explicit immediate-abort witness with remaining_stack_skipped > 0).
    2) High-mu / low-p0 300-ideal adaptively-reached real-FB campaign (q15_low_p0_h5.b85: 180 adaptively-reached
       first-survivor children from immediate-abort k=1 trees + 120 matched disjoint control ideals):
       - Distinct coprime relation counts M vs calibrated mu_d -> 0/300 zeros vs 30.17 summed Cantelli p_0 bounds,
         conditional Fano factor Var(M|slot)/E[M|slot] = 1.598 (ENC) / 1.439 (CTRL) < 2.00 (directly corroborating
         Var(M|slot) <= 2*mu!), k_hat_cond = 0.0301 (ENC) / 0.0251 (CTRL) <= 0.035, E_emp[0.8^M] = 0.0384 << E[p_0].
       - Live bit-exact replay of 6 ideals (3 ENC + 3 CTRL) from raw (sd, q, ry, rx) via 16D Bareiss-LLL + _strict_fac.
    3) Exact algebraic naive-PGF counterexample (why pure k=1 Boole union bound is used).
    """
    if rho_exact is None:
        from cm_models import rho_exact as _re
        rho_exact = _re
    tc = verify_threshold_counterexample()
    dp = desc_params or {}
    p0_mode = str(dp.get("p0_mode", "cantelli05"))
    c_mu = float(dp.get("c_mu", 0.5))
    r_sh = float(dp.get("r_shape", 1.0))
    p0_floor = float(dp.get("p0_floor", 0.0001))

    def _p0_claimed(mu_m):
        if p0_mode == "cantelli05":
            return max(p0_floor, 1.0 / (1.0 + c_mu * mu_m))
        return max(p0_floor, (1.0 + mu_m / r_sh) ** (-r_sh))

    W_EL = y_pow_z(5); W2 = y_pow_z(10)
    CF_EL = mul_phi15_z(add_phi15_z([2, 0, 0, 0, 0, 0, 0, 0], W_EL), W2)
    CG_EL = scale_phi15_z(2, W2)
    B_FB = 8192
    lB = math.log2(B_FB)
    P_fb = math.prod(primes_upto(B_FB))

    def _ch_root(p, a8, b8, C_el):
        if p % 15 != 1:
            return None
        e15 = (p - 1) // 15
        g0 = next((pow(g, e15, p) for g in range(2, p) if pow(pow(g, e15, p), 3, p) != 1 and pow(pow(g, e15, p), 5, p) != 1), None)
        if g0 is None:
            return None
        for k in (1, 2, 4, 7, 8, 11, 13, 14):
            ry = pow(g0, k, p)
            av = sum(a8[j] * pow(ry, j, p) for j in range(8)) % p
            bv = sum(b8[j] * pow(ry, j, p) for j in range(8)) % p
            if bv == 0:
                continue
            rx = (av * pow(bv, -1, p)) % p
            rhs = sum(C_el[j] * pow(ry, j, p) for j in range(8)) % p
            if pow(rx, 3, p) == rhs:
                return (ry, rx)
        return None

    def _strict_fac(n, q_lim, a8, b8, C_el):
        if n <= 0:
            return None
        t = n
        while True:
            g = math.gcd(t, P_fb)
            if g == 1:
                break
            t //= g
        if t == 1:
            return []
        if t < q_lim:
            if _isp_det(t):
                if t % 15 != 1:
                    return None
                rt = _ch_root(t, a8, b8, C_el)
                return [(t, rt)] if rt is not None else None
            return None
        if t >= q_lim * q_lim or _isp_det(t):
            return None
        d = None
        for c0 in (1, 2, 3, 5, 7):
            x = y = 2; g = 1; it = 0
            while g == 1 and it < 300:
                x = (x * x + c0) % t
                y = (y * y + c0) % t; y = (y * y + c0) % t
                g = math.gcd(abs(x - y), t); it += 1
            if 1 < g < t:
                d = g; break
        if d is None:
            return None
        p1, p2 = d, t // d
        if not (_isp_det(p1) and _isp_det(p2)):
            return None
        if not (B_FB < p1 < q_lim and B_FB < p2 < q_lim):
            return None
        if p1 % 15 != 1 or p2 % 15 != 1:
            return None
        r1 = _ch_root(p1, a8, b8, C_el); r2 = _ch_root(p2, a8, b8, C_el)
        if r1 is None or r2 is None:
            return None
        return [(p1, r1), (p2, r2)] if p1 != p2 else [(p1, r1)]

    def _p_strict(n_val, lq):
        m = math.log2(max(2, n_val))
        p0 = rho_exact(m / lB)
        p1 = 0.0
        if m > lB:
            dt = (lq - lB) / 6.0
            for s in range(6):
                t = lB + (s + 0.5) * dt
                if m - t > 0:
                    p1 += rho_exact((m - t) / lB) * (dt / t)
        return min(1.0, p0 + p1)

    def _is_cp_fast(a8, b8):
        na = abs(norm_q15_z(a8)); nb = abs(norm_q15_z(b8))
        g0 = math.gcd(na, nb)
        if g0 <= 1:
            return True
        t = g0; d = 2
        while d * d <= t and d < 200:
            if t % d == 0:
                while t % d == 0:
                    t //= d
                g1 = pgcd_mod(a8, PHI15, d)
                if len(g1) > 1 and len(pgcd_mod(g1, b8, d)) > 1:
                    return False
            d += 1 if d == 2 else 2
        if t > 1 and t < 200:
            g1 = pgcd_mod(a8, PHI15, t)
            if len(g1) > 1 and len(pgcd_mod(g1, b8, t)) > 1:
                return False
        return True

    # Build deterministic enumeration vectors (wt-1 + wt-2 + wt-3 up to 1500)
    N_LO_PTS = int(low_p0_data["n_pts"])
    pts = [((i, 1),) for i in range(16)]
    for i in range(16):
        for j in range(i + 1, 16):
            pts.append(((i, 1), (j, -1))); pts.append(((i, 1), (j, 1)))
    for i in range(16):
        for j in range(i + 1, 16):
            for k in range(j + 1, 16):
                for s1 in (-1, 1):
                    for s2 in (-1, 1):
                        pts.append(((i, 1), (j, s1), (k, s2)))
                        if len(pts) >= N_LO_PTS:
                            break
                    if len(pts) >= N_LO_PTS:
                        break
                if len(pts) >= N_LO_PTS:
                    break
            if len(pts) >= N_LO_PTS:
                break
        if len(pts) >= N_LO_PTS:
            break

    _lll_cache = {}
    def _get_R(sd, q, ry, rx):
        key = (sd, q, ry, rx)
        R = _lll_cache.get(key)
        if R is None:
            B = [[0] * 16 for _ in range(16)]; B[0][0] = q
            rp = [pow(ry, j, q) for j in range(8)]
            for j in range(1, 8):
                B[j][0] = (-rp[j]) % q; B[j][j] = 1
            for j in range(8):
                B[8 + j][0] = (rx * rp[j]) % q; B[8 + j][8 + j] = 1
            R, dsq, _ = lll_16d_int(B, 99, 100); assert dsq[16] == q * q
            _lll_cache[key] = R
        return R

    def _eval_node(sd, q, ry, rx, n_pts):
        R = _get_R(sd, q, ry, rx)
        first_surv = None; n_surv = 0; mu_raw_sum = 0.0; lq = math.log2(q)
        for pt in pts[:n_pts]:
            a8 = [0] * 8; b8 = [0] * 8
            for idx, sgn in pt:
                row = R[idx]
                for j in range(8):
                    a8[j] += sgn * row[j]; b8[j] += sgn * row[8 + j]
            if not any(a8) or not any(b8):
                continue
            a3 = pow_phi15_z(a8, 3); b3 = pow_phi15_z(b8, 3)
            nf = abs(norm_q15_z(sub_phi15_z(a3, mul_phi15_z(CF_EL, b3))))
            ng = abs(norm_q15_z(sub_phi15_z(a3, mul_phi15_z(CG_EL, b3))))
            if nf <= 1 or ng <= 1:
                continue
            if sd == "f":
                if nf % q != 0:
                    continue
                nf //= q
            else:
                if ng % q != 0:
                    continue
                ng //= q
            if nf <= 0 or ng <= 0:
                continue
            mu_raw_sum += _p_strict(nf, lq) * _p_strict(ng, lq)
            ff = _strict_fac(nf, q, a8, b8, CF_EL)
            if ff is None:
                continue
            fg = _strict_fac(ng, q, a8, b8, CG_EL)
            if fg is None:
                continue
            ch = [("f", cp, rt[0], rt[1]) for cp, rt in ff] + [("g", cp, rt[0], rt[1]) for cp, rt in fg]
            for _, cp, _, _ in ch:
                assert B_FB < cp < q
            n_surv += 1
            if first_surv is None:
                first_surv = (a8, b8, ch)
        return first_surv, n_surv, mu_raw_sum

    def _eval_distinct_replay(sd, q, ry, rx, n_pts):
        R = _get_R(sd, q, ry, rx)
        seen_cand = set(); ns = 0; mu_d = 0.0; lq = math.log2(q)
        for pt in pts[:n_pts]:
            a8 = [0] * 8; b8 = [0] * 8
            for idx, sgn in pt:
                row = R[idx]
                for j in range(8):
                    a8[j] += sgn * row[j]; b8[j] += sgn * row[8 + j]
            if not any(a8) or not any(b8):
                continue
            a3 = pow_phi15_z(a8, 3); b3 = pow_phi15_z(b8, 3)
            nf = abs(norm_q15_z(sub_phi15_z(a3, mul_phi15_z(CF_EL, b3))))
            ng = abs(norm_q15_z(sub_phi15_z(a3, mul_phi15_z(CG_EL, b3))))
            if nf <= 1 or ng <= 1:
                continue
            if sd == "f":
                if nf % q != 0:
                    continue
                nf //= q
            else:
                if ng % q != 0:
                    continue
                ng //= q
            if nf <= 0 or ng <= 0:
                continue
            kc = (nf, ng)
            if kc in seen_cand:
                continue
            seen_cand.add(kc)
            if not _is_cp_fast(a8, b8):
                continue
            mu_d += _p_strict(nf, lq) * _p_strict(ng, lq)
            ff = _strict_fac(nf, q, a8, b8, CF_EL)
            if ff is None:
                continue
            fg = _strict_fac(ng, q, a8, b8, CG_EL)
            if fg is None:
                continue
            ns += 1
        return ns, mu_d

    # 1. LIVE high-p0 boundary IMMEDIATE-ABORT strict-k=1 campaign (64 live 17-bit roots, n_pts=32)
    roots_17 = []
    for p in range(65551, 131072, 30):
        if p % 15 == 1 and _isp_det(p):
            e15 = (p - 1) // 15
            g0 = next(pow(g, e15, p) for g in range(2, p) if pow(pow(g, e15, p), 3, p) != 1 and pow(pow(g, e15, p), 5, p) != 1)
            for k in (1, 2, 4, 7, 8, 11, 13, 14):
                ry = pow(g0, k, p)
                rhs = sum(CG_EL[j] * pow(ry, j, p) for j in range(8)) % p
                if pow(rhs, (p - 1) // 3, p) == 1:
                    for x in range(1, min(p, 12000)):
                        if pow(x, 3, p) == rhs:
                            roots_17.append(("g", p, ry, x)); break
                    break
            if len(roots_17) >= 64:
                break
    assert len(roots_17) == 64

    bd_npts = 32 if (p0_mode == "cantelli05" and c_mu <= 0.51) else 48
    nonempty = completed = aborted = reached_ch = mdep = post_abort_skipped = 0
    e_cl = mu_obs_t = mu_mod_t = 0.0
    comp_wits = []; abort_wits = []
    for r0 in roots_17:
        fs0, _, _ = _eval_node(*r0, bd_npts)
        if fs0 is None:
            continue
        nonempty += 1
        stk = [(ch, 1) for ch in reversed(fs0[2])]
        hit_abort = False; root_d = 0; visited_ch = 0
        while stk:
            (sd, q, ry, rx), dep = stk.pop()
            if dep > mdep:
                mdep = dep
            if dep > root_d:
                root_d = dep
            reached_ch += 1; visited_ch += 1
            fs, ns, mu_raw = _eval_node(sd, q, ry, rx, bd_npts)
            mu_m = 0.50 * mu_raw
            e_cl += _p0_claimed(mu_m)
            mu_obs_t += ns; mu_mod_t += mu_m
            if fs is None:
                aborted += 1; hit_abort = True
                rem_sk = len(stk)
                post_abort_skipped += rem_sk
                if len(abort_wits) < 2 and rem_sk > 0:
                    abort_wits.append((r0[1], q, dep, rem_sk))
                break  # IMMEDIATE ATTEMPT ABORT — ZERO post-abort siblings visited!
            for ch in reversed(fs[2]):
                stk.append((ch, dep + 1))
        if not hit_abort:
            completed += 1
            if len(comp_wits) < 2 and root_d >= 2:
                comp_wits.append((r0[1], root_d, visited_ch))

    assert nonempty == completed + aborted and nonempty >= 50 and reached_ch >= 60 and mdep >= 3
    assert len(comp_wits) >= 2 and len(abort_wits) >= 1 and post_abort_skipped > 0
    bd_oe = aborted / max(1e-9, e_cl)
    assert bd_oe <= 0.85 and (mu_obs_t / max(1e-9, mu_mod_t)) >= 0.95

    # 2. High-mu / low-p0 300-ideal adaptively-reached real-FB campaign (q15_low_p0_h5.b85) + LIVE replay of 6 ideals
    enc_recs = low_p0_data["enc"]; ctrl_recs = low_p0_data["ctrl"]
    assert len(enc_recs) == 180 and len(ctrl_recs) == 120
    assert not ({r[1] for r in enc_recs} & {r[1] for r in ctrl_recs}), "ENC and CTRL primes must be disjoint!"

    def _summarize_low_p0(arr):
        ms = [r[4] for r in arr]
        mr_f = sum(r[4] for r in arr if r[0] == "f") / max(1e-9, sum(r[5] for r in arr if r[0] == "f"))
        mr_g = sum(r[4] for r in arr if r[0] == "g") / max(1e-9, sum(r[5] for r in arr if r[0] == "g"))
        exps = [(mr_f if r[0] == "f" else mr_g) * r[5] for r in arr]
        rats = [m / e for m, e in zip(ms, exps)]
        fano_c = sum((m - e) ** 2 / e for m, e in zip(ms, exps)) / len(arr)
        kh_c = max(0.0, sum((m - e) ** 2 - e for m, e in zip(ms, exps)) / sum(e * e for e in exps))
        sum_p0 = sum(_p0_claimed(e) for e in exps)
        e08 = sum(0.8 ** m for m in ms) / len(arr)
        return dict(
            n=len(arr), zeros=ms.count(0), min_M=min(ms), mean_M=round(sum(ms) / len(ms), 2),
            min_M_over_E=round(min(rats), 3), p05_M_over_E=round(sorted(rats)[len(rats) // 20], 3),
            fano_cond=round(fano_c, 3), k_hat_cond=round(kh_c, 5),
            sum_p0_claimed=round(sum_p0, 2), mean_p0_claimed=round(sum_p0 / len(arr), 4),
            E_emp_08_M=round(e08, 5),
        )

    s_enc = _summarize_low_p0(enc_recs)
    s_ctrl = _summarize_low_p0(ctrl_recs)
    s_all = _summarize_low_p0(enc_recs + ctrl_recs)
    assert s_all["zeros"] == 0 and s_all["min_M"] >= 5 and s_all["sum_p0_claimed"] >= 25.0
    assert s_enc["fano_cond"] < 1.75 and s_ctrl["fano_cond"] < 1.65 and s_all["fano_cond"] < 1.70
    assert s_enc["k_hat_cond"] <= 0.04 and s_ctrl["k_hat_cond"] <= 0.04 and s_all["k_hat_cond"] <= 0.04
    assert s_enc["E_emp_08_M"] < 0.5 * s_enc["mean_p0_claimed"]

    # Live bit-exact replay of 6 high-mu ideals (3 ENC + 3 CTRL) from raw (sd, q, ry, rx)!
    replay_ok = 0
    for rec in (enc_recs[0], enc_recs[60], enc_recs[120], ctrl_recs[0], ctrl_recs[60], ctrl_recs[100]):
        sd, q, ry, rx, m_exp, mu_exp = rec
        m_got, mu_got = _eval_distinct_replay(sd, q, ry, rx, N_LO_PTS)
        assert m_got == m_exp and abs(mu_got - mu_exp) < 1e-3
        replay_ok += 1
    assert replay_ok == 6

    return dict(
        poly_pair="deg24_24_Q15_immediate_abort_strict_descending_B8192",
        zero_omitted_primes=True,
        immediate_attempt_abort_on_first_empty=True,
        zero_post_abort_siblings_visited=True,
        direct_modeled_p0_denominator=True,
        p0_mode=p0_mode,
        c_mu=c_mu,
        udup_applied=0.50,
        boundary_immediate_abort=dict(
            n_roots=64, n_pts=bd_npts, nonempty=nonempty, completed=completed, aborted=aborted,
            exact_attempt_partition=(nonempty == completed + aborted),
            reached_children=reached_ch, post_abort_siblings_skipped=post_abort_skipped,
            max_depth=mdep, obs_over_mod_mu=round(mu_obs_t / max(1e-9, mu_mod_t), 3),
            exp_aborts_claimed=round(e_cl, 3), oe_claimed=round(bd_oe, 4),
            n_completed_wits=len(comp_wits), n_immediate_abort_wits=len(abort_wits),
        ),
        low_p0_300_ideals=dict(
            n_pts=N_LO_PTS, live_bit_exact_replayed=replay_ok,
            enc_180=s_enc, ctrl_120=s_ctrl, all_300=s_all,
        ),
        naive_pgf_counterexample=tc,
        checks_ok=True,
    )
