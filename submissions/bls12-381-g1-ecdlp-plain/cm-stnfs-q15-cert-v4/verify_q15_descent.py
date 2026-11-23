#!/usr/bin/env python3
"""verify_q15_descent.py - pure-Python stdlib verifier for the Q(zeta_15) descent (certified k = 1 Boole + Wald boot).

Checks:
  [7]  16D ball half-counts (R = 2, 4, 256) and the exact 384-mesh evaluation of PLAN16/BOOT16 with one certified trial
       per slot (K_TRIALS = 1, P_NODE = exp(-MU_CERT) + DELTA_CERT, certificate cost charged per draw).
  [7c] Yield certificate (q15_cert16.py): archived production campaign q15_cert_prod.b85 (76 entries x 6 ideals x 2
       p_rand x 256 cells) re-checked record by record, plus live recomputation of one lattice per side from its seed.
  [7p] Live exact count of the 4,740 known-log f_9 ideals in [2^14, 2^16) (the p_rand pool).
  [11] Live integer BOUND_CELL on 12 sampled cells, (135, g) short-subspace rejection, multi-bin UDUP check.
  [11b] Live q15_retain16 prime-splitting + p = 2 (K = 140) check + j >= 74 tail + counted_fraction.
  [12] Live q15_ecm16 L_min = 2 branch-tree DP + 12-bit geometric ladder + Suyama/Newton/singular-N-1 check.
"""
import base64
import bisect
import json
import lzma
import math
import random
import sys
import time
from math import ceil, log2
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import q15_descent16 as DS
import q15_counts16 as CT
import q15_ecm16 as EC
import q15_retain16 as RT
import q15_bound_cell as BC
import q15_cert16 as QC

_SUMPROD = getattr(math, "sumprod", lambda a, b: sum(x * y for x, y in zip(a, b)))

U_BLS = -0xd201000000010000
LOG2U8 = 8.0 * log2(abs(U_BLS))
PHI15 = [1, -1, 0, 1, -1, 1, 0, -1, 1]


def _mul_z(a, b):
    t = [0] * 15
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                t[i + j] += x * y
    for k in range(14, 7, -1):
        c = t[k]
        if c:
            o = k - 8
            t[o] -= c; t[o + 1] += c; t[o + 3] -= c
            t[o + 4] += c; t[o + 5] -= c; t[o + 7] += c
            t[k] = 0
    return t[:8]


_YPOW = [[1, 0, 0, 0, 0, 0, 0, 0]]
for _i in range(1, 15):
    _YPOW.append(_mul_z(_YPOW[-1], [0, 1, 0, 0, 0, 0, 0, 0]))
_W = _YPOW[5]
_WP1 = [(1 if i == 0 else 0) + _W[i] for i in range(8)]
_W2 = _YPOW[10]


def _norm_q15(a):
    if not any(a):
        return 0
    def sig(v, k):
        res = [0] * 8
        for j, cj in enumerate(v):
            if cj:
                pw = _YPOW[(k * j) % 15]
                for i in range(8):
                    res[i] += cj * pw[i]
        return res
    p = _mul_z(a, sig(a, 4))
    p = _mul_z(p, sig(p, 2))
    return _mul_z(p, sig(p, 11))[0]


def _norms_fg(aa, bb):
    a2, b2 = _mul_z(aa, aa), _mul_z(bb, bb)
    a3, b3 = _mul_z(a2, aa), _mul_z(b2, bb)
    w2b3 = _mul_z(_W2, b3)
    ng = abs(_norm_q15([a3[i] - U_BLS * w2b3[i] for i in range(8)]))
    a6, b6 = _mul_z(a3, a3), _mul_z(b3, b3)
    a9, b9 = _mul_z(a6, a3), _mul_z(b6, b3)
    t1 = _mul_z(_WP1, _mul_z(a6, b3))
    t2 = _mul_z(a3, b6)
    t3 = _mul_z(_W, b9)
    nf = abs(_norm_q15([a9[i] + t1[i] + t2[i] - t3[i] for i in range(8)]))
    return nf, ng


def _is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n == p:
            return True
        if n % p == 0:
            return False
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if a >= n:
            continue
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


def _tr15(d):
    g = math.gcd(15, d % 15)
    return 8 if g == 15 else (-4 if g == 5 else (-2 if g == 3 else 1))


_G8 = [[_tr15(i - j) for j in range(8)] for i in range(8)]


def _lll_16d_qcm(B_in, num=99, den=100):
    def _bareiss(B0, use_g8, p_num=num, p_den=den):
        n = 16
        M = [list(r) for r in B0]
        sp = _SUMPROD
        if use_g8:
            def dot(u, v):
                s = 0
                for i in range(8):
                    ui, ub = u[i], u[8 + i]
                    gi = _G8[i]
                    if ui:
                        s += ui * sp(gi, v[:8])
                    if ub:
                        s += ub * sp(gi, v[8:])
                return s
        else:
            dot = sp
        d = [1] * (n + 1)
        lam = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1):
                u = dot(M[i], M[j])
                for k in range(j):
                    u = (d[k + 1] * u - lam[i][k] * lam[j][k]) // d[k]
                if j < i:
                    lam[i][j] = u
                else:
                    d[i + 1] = u
        k = 1
        while k < n:
            for l in range(k - 1, -1, -1):
                if 2 * abs(lam[k][l]) > d[l + 1]:
                    q = (2 * lam[k][l] + d[l + 1]) // (2 * d[l + 1]) if lam[k][l] > 0 else -((-2 * lam[k][l] + d[l + 1]) // (2 * d[l + 1]))
                    if q:
                        M[k] = [x - q * y for x, y in zip(M[k], M[l])]
                        lam[k][l] -= q * d[l + 1]
                        for j in range(l):
                            lam[k][j] -= q * lam[l][j]
            if p_den * (d[k - 1] * d[k + 1] + lam[k][k - 1] ** 2) < p_num * d[k] * d[k]:
                M[k], M[k - 1] = M[k - 1], M[k]
                for j in range(k - 1):
                    lam[k][j], lam[k - 1][j] = lam[k - 1][j], lam[k][j]
                for i in range(k + 1, n):
                    t = lam[i][k - 1]
                    lam[i][k - 1] = (lam[k][k - 1] * t + d[k - 1] * lam[i][k]) // d[k]
                    lam[i][k] = (d[k + 1] * t - lam[k][k - 1] * lam[i][k]) // d[k]
                d[k] = (d[k - 1] * d[k + 1] + lam[k][k - 1] ** 2) // d[k]
                if k > 1:
                    k -= 1
            else:
                k += 1
        return M, d

    B1, d1 = _bareiss(B_in, False, 11, 40)
    B2, d2 = _bareiss(B1, True, num, den)
    return B2, d1[16], d2[16]


def _sq_hnf(q, s, t):
    rp = [pow(s, i, q) for i in range(8)]
    B = [[0] * 16 for _ in range(16)]
    B[0][0] = q
    for i in range(1, 8):
        B[i][0] = (-rp[i]) % q; B[i][i] = 1
    for i in range(8):
        B[8 + i][0] = (t * rp[i]) % q; B[8 + i][8 + i] = 1
    return B


def _roots_mod_q(poly, q, rng):
    f = RT._trim([c % q for c in poly])
    if len(f) <= 1:
        return []
    xp = RT.ppowmod([0, 1], q, f, q)
    if len(xp) < 2:
        xp += [0] * (2 - len(xp))
    xp[1] = (xp[1] - 1) % q
    g = RT.pgcd(f, RT._trim(xp), q)
    out = []

    def split(h):
        deg = len(h) - 1
        if deg <= 0:
            return
        if deg == 1:
            out.append((-h[0] * pow(h[1], -1, q)) % q)
            return
        while True:
            a = rng.randrange(q)
            t = RT.ppowmod([a, 1], (q - 1) // 2, h, q)
            if not t:
                continue
            t[0] = (t[0] - 1) % q
            d = RT.pgcd(h, RT._trim(t), q)
            if 0 < len(d) - 1 < deg:
                quo, _ = RT.pdiv(h, d, q)
                split(d)
                split(quo)
                return

    split(g)
    return sorted(set(out))


def _sample_ideal(bits, side, rng):
    while True:
        q = rng.randrange(1 << (bits - 1), 1 << bits)
        q = q - (q % 15) + 1
        if q.bit_length() != bits or not _is_prime(q):
            continue
        e15 = (q - 1) // 15
        g0 = next((pow(g, e15, q) for g in range(2, 100) if pow(pow(g, e15, q), 3, q) != 1 and pow(pow(g, e15, q), 5, q) != 1), None)
        if g0 is None:
            continue
        s = pow(g0, rng.choice((1, 2, 4, 7, 8, 11, 13, 14)), q)
        w = pow(s, 5, q)
        poly = [(-U_BLS * w * w) % q, 0, 0, 1] if side == "g" else [(-w) % q, 0, 0, 1, 0, 0, (w + 1) % q, 0, 0, 1]
        rts = _roots_mod_q(poly, q, rng)
        if rts:
            return q, s, rng.choice(rts)


def _count_f9_ideals_p1416():
    """Exact live count of degree-1 f_9 ideals with prime norm p in [2^14, 2^16) (p === 1 mod 15)."""
    cnt_f9 = 0
    n_split = 0
    for p in range(1 << 14, 1 << 16):
        if p % 15 != 1 or not _is_prime(p):
            continue
        n_split += 1
        w = next(pow(g, (p - 1) // 3, p) for g in range(2, 100) if pow(g, (p - 1) // 3, p) != 1)
        w2 = (w * w) % p
        for wi in (w, w2):
            poly = [(-wi) % p, 0, 0, 1, 0, 0, (wi + 1) % p, 0, 0, 1]
            xp = RT.ppowmod([0, 1], p, poly, p)
            if len(xp) < 2:
                xp += [0] * (2 - len(xp))
            xp[1] = (xp[1] - 1) % p
            g = RT.pgcd(poly, RT._trim(xp), p)
            cnt_f9 += 4 * (len(g) - 1)
    return cnt_f9, n_split


def _build_dp_ball(rmax=32):
    nmax = rmax * rmax
    rd = [[0] * (nmax + 1) for _ in range(17)]
    rd[0][0] = 1
    for d in range(1, 17):
        prev = rd[d - 1]; cur = rd[d]
        for s in range(nmax + 1):
            v = prev[s]
            if not v:
                continue
            cur[s] += v
            k = 1
            while s + k * k <= nmax:
                cur[s + k * k] += 2 * v
                k += 1
    pref = [0] * (nmax + 1)
    acc = 0
    for s in range(1, nmax + 1):
        acc += rd[16][s]
        pref[s] = acc
    return rd, pref


def _sample_prim_ball_dp(rd, pref, R, rng):
    r2 = R * R
    tot = pref[r2]
    while True:
        u = rng.randrange(tot)
        s = bisect.bisect_right(pref, u, 1, r2 + 1)
        c = [0] * 16
        for d in range(16, 0, -1):
            rem = s
            v = rng.randrange(rd[d][rem])
            w0 = rd[d - 1][rem]
            if v < w0:
                c[d - 1] = 0
            else:
                v -= w0
                k = 1
                while True:
                    wk = rd[d - 1][rem - k * k]
                    if v < wk:
                        c[d - 1] = k; s -= k * k; break
                    v -= wk
                    if v < wk:
                        c[d - 1] = -k; s -= k * k; break
                    v -= wk
                    k += 1
        if math.gcd(*c) == 1:
            for x in reversed(c):
                if x:
                    if x < 0:
                        c = [-y for y in c]
                    break
            return c


def run_all_checks(root_dir=None, seed=1, mesh=384):
    t0 = time.perf_counter()
    rdir = Path(root_dir) if root_dir else HERE
    st_b85 = (rdir / "q15_stats.b85").read_text()
    st = json.loads(lzma.decompress(base64.b85decode("".join(st_b85.split()))).decode())
    c7e_b85 = (rdir / "q15_cond_7e.b85").read_text()
    c7e = json.loads(lzma.decompress(base64.b85decode("".join(c7e_b85.split()))).decode())

    # [7] Exact 384-mesh evaluation of PLAN16/BOOT16 (certified single trial per slot, Boole, Wald boot)
    Pm = DS.P()
    DS.load_stats(st, Pm)
    assert Pm.LB == 68.5 and Pm.LX == 68.5 and Pm.BP_LO == 14.0 and Pm.BP_HI == 16.0
    assert Pm.K_TRIALS == 1 and Pm.MU_MIN == 15.0 and Pm.MU_CERT == 15.0 and Pm.DELTA_CERT == 1.0e-9
    assert Pm.P_EMPTY == 1.0e-4 and Pm.P_NODE == Pm.P_EMPTY + Pm.DELTA_CERT and math.exp(-Pm.MU_CERT) < Pm.P_EMPTY / 300.0 and Pm.CERT_DRAW == 2.0 ** 34 and Pm.CERT_PASS_LO == 0.5 and Pm.D_MAX == 64 and Pm.K_ATT == 6 and Pm.P_BAD == 0.03 and Pm.P_BAD_BOOT == 0.01 and Pm.P_LOW == 0.01
    assert Pm.LMAX == 260 and Pm.BALL_PAIR == 0.0 and Pm.BOOT_KEEP_FRAC == 0.90
    assert Pm.BOUND_CELL == 5000.0 and Pm.SIEVE_HIT == 2600.0 and Pm.UDUP == 0.40
    Pm.MESH = int(mesh); Pm.TRUNC = 0.0
    Pm.MAIN_LB = 71.19; Pm.MAIN_UF = 546.0 / 71.19; Pm.MAIN_UG = 670.0 / 71.19

    assert DS.CELLS[2] == 17072 and DS.CELLS[4] == 642328192 and DS.CELLS[8] == 35245626789232
    vol256 = 0.5 * (math.pi ** 8) / math.factorial(8) * (256 ** 16)
    assert abs(DS.CELLS[256] / vol256 - 1.0) < 1e-4

    # Strict +-16-ULP Dickman enclosure
    assert DS.rho_lo(2.0) < (1.0 - math.log(2.0)) < DS.rho_hi(2.0)
    assert DS.rho_lo(1032.0 / 1024.0) < (1.0 - math.log(1032.0 / 1024.0)) < DS.rho_hi(1032.0 / 1024.0)
    assert DS.rho_lo(3080.0 / 1024.0) < DS.rho_series(3080.0 / 1024.0) < DS.rho_hi(3080.0 / 1024.0)
    assert DS.rho_lo(3306.0 / 1024.0) < DS.rho_series(3306.0 / 1024.0) < DS.rho_hi(3306.0 / 1024.0)

    # Collection-cut-conditioned miss_side strict monotonicity at LX = 68.5
    _mf_seq = [DS.miss_side(b_m, Pm.LOG2_RELS, Pm.MAIN_UF, Pm.LX, Pm.MAIN_LB, Pm.ORBIT, Pm.P_SIEVE_F, Pm.RET_MULT_F) for b_m in range(40, 71)]
    _mg_seq = [DS.miss_side(b_m, Pm.LOG2_RELS, Pm.MAIN_UG, Pm.LX, Pm.MAIN_LB, Pm.ORBIT, Pm.P_SIEVE_G, Pm.RET_MULT_G) for b_m in range(40, 71)]
    _mg_uncond = DS.miss_side(40, Pm.LOG2_RELS, Pm.MAIN_UG, Pm.LX, Pm.MAIN_LB, Pm.ORBIT, Pm.LX, 1.0)
    assert _mf_seq[0] > 0.0 and _mg_seq[0] > _mg_uncond > 0.0
    assert all(_mf_seq[i] < _mf_seq[i + 1] and _mg_seq[i] < _mg_seq[i + 1] for i in range(len(_mf_seq) - 1))

    t_ev = time.perf_counter()
    res = DS.evaluate(Pm, DS.PLAN16, DS.BOOT16)
    dt_ev = time.perf_counter() - t_ev
    rows = DS.report(res)
    bt = res["boot"]
    EN, _ = CT.expected_nodes(Pm, res)
    two_log = res["two_log_units"]
    min_mu = min(r["mu"] for r in rows)
    min_ret = min(r["ret"] for r in rows)
    max_loss = max(max(r["loss_f"], r["loss_g"]) for r in rows)
    max_miss_node = max(d["miss"] for d in res["DET"].values())

    assert len(rows) == 76
    assert min_mu >= 15.0 and bt["phi_tree"] < 0.37 and bt["phi"] == 0.0 and two_log <= 107.52
    assert min_ret >= 0.90 and max_loss < 5e-3
    assert all(r["how"] in ("sieve", "ecm") for r in rows)
    _det_sieve = [v for v in res["DET"].values() if v["how"] == "sieve"]
    _det_ecm = [v for v in res["DET"].values() if v["how"] == "ecm"]
    assert len(_det_sieve) == 71 and len(_det_ecm) == 5
    assert all(v["cofac_frac"] == 1.0 and v["eager"] is True for v in res["DET"].values())
    assert all(abs(v["lin_keep"] - max(0.0, 1.0 - v["miss"])) < 1e-12 for v in res["DET"].values())
    assert abs(bt["lin_keep"] - max(0.0, 1.0 - bt["miss"])) < 1e-12

    # [7c] Yield certificate: production campaign (q15_cert_prod.b85), the algorithm's exact trial law (uniform p_rand from
    # the 4,740-ideal pool, exact orbit sizes, admission rules), first draws at delta_1 = DELTA_CERT / 2.
    cp = json.loads(lzma.decompress(base64.b85decode("".join((rdir / "q15_cert_prod.b85").read_text().split()))).decode())
    assert cp["S"] == QC.S_CERT and cp["delta"] == Pm.DELTA_CERT
    keys = sorted(res["DET"].keys())
    assert abs(QC.draw_delta(1, Pm.DELTA_CERT) - Pm.DELTA_CERT / 2) < 1e-24 and sum(QC.draw_delta(i, 1.0) for i in range(1, 100001)) < 1.0
    assert QC.prand_pool(sys.modules[__name__], RT)[1] == 4740
    cpb = json.loads(lzma.decompress(base64.b85decode("".join((rdir / "q15_cert_boot.b85").read_text().split()))).decode())
    cp["seed_boot"] = cpb["seed_boot"]
    ent, bent = cp["entries"], cpb["boot_entries"]
    assert sorted(int(k) for k in ent) == list(range(76)) and all(len(e["draws"]) == 200 for e in ent.values())
    boot_keys = [i for i, k in enumerate(keys) if k[1] == "g" and k[0][0] < float(DS.BOOT16[0])]
    assert sorted(int(k) for k in bent) == boot_keys and all(len(e["draws"]) == 400 for e in bent.values())
    fails = sum(1 for e in list(ent.values()) + list(bent.values()) for s_ok, cv, _ in e["draws"] if not (s_ok == 1 and cv >= Pm.MU_CERT))
    assert fails == 0
    cp95 = lambda n_: 2.0 * (1.0 - 0.05 ** (1.0 / n_))
    p_bad_entry = cp95(200)
    p_bad_boot = cp95(600)
    assert p_bad_entry <= Pm.P_BAD and p_bad_boot <= Pm.P_BAD_BOOT and 2.0 ** -Pm.D_MAX <= 1e-18
    ratios = sorted(cv / res["DET"][keys[int(k)]]["mu"] for k, e in ent.items() for _, cv, _ in e["draws"])
    mratios = [mh / res["DET"][keys[int(k)]]["mu"] for k, e in ent.items() for _, _, mh in e["draws"]]
    cert_min = min(cv for e in list(ent.values()) + list(bent.values()) for _, cv, _ in e["draws"])
    assert cert_min >= 4.0 * Pm.MU_CERT
    live_cert = []
    for ei, sd, src in ((50, cp["seed"], ent), (53, cp["seed"], ent), (boot_keys[0], cp["seed_boot"], bent)):
        lr = QC.certify_entry(DS, sys.modules[__name__], BC, RT, Pm, res, keys[ei], 1, 1, sd, first_only=True)[0]
        ar = src[str(ei)]["first"]
        assert lr[:5] == ar[:5] and lr[8:] == ar[8:] and all(abs(lr[k] - ar[k]) <= 1e-3 * max(1.0, abs(ar[k])) for k in (5, 6, 7))
        live_cert.append(lr)
    # [7t] Toy lattice-conditional law: per-lattice [mu_L, M_L] (q15_toy_cond.b85)
    tc = json.loads(lzma.decompress(base64.b85decode("".join((rdir / "q15_toy_cond.b85").read_text().split()))).decode())["cohorts"]
    toy = {}
    pois3 = lambda u: math.exp(-u) * (1.0 + u + u * u / 2.0 + u ** 3 / 6.0)
    for nm, rows_t in tc.items():
        mus = [r[0] for r in rows_t]; Ms = [r[1] for r in rows_t]
        D_t = sum((m - u) ** 2 for u, m in zip(mus, Ms)) / sum(mus)
        hi10 = [(u, m) for u, m in zip(mus, Ms) if u >= 10.0]
        hi15 = [(u, m) for u, m in zip(mus, Ms) if u >= 15.0]
        toy[nm] = dict(n=len(rows_t), D=round(D_t, 3), zeros=sum(1 for m in Ms if m == 0), exp_zeros_poisson=round(sum(math.exp(-u) for u in mus), 2),
                       n_mu_ge10=len(hi10), zeros_mu_ge10=sum(1 for _, m in hi10 if m == 0), le3_mu_ge10=sum(1 for _, m in hi10 if m <= 3),
                       exp_le3_mu_ge10=round(sum(pois3(u) for u, _ in hi10), 2), n_mu_ge15=len(hi15), le3_mu_ge15=sum(1 for _, m in hi15 if m <= 3))
    assert toy["c2Q"]["n"] == toy["c1A"]["n"] == toy["c1C"]["n"] == 1536
    assert all(t["D"] < 1.4 and t["zeros_mu_ge10"] == 0 and t["le3_mu_ge15"] == 0 for t in toy.values())
    assert sum(t["n_mu_ge10"] for t in toy.values()) == 2646 and sum(t["n_mu_ge15"] for t in toy.values()) == 1163

    # [7p] Live exact count of the 4,740 f_9 ideals in [2^14, 2^16) (p_rand pool)
    n_pool_f9, n_split_primes = _count_f9_ideals_p1416()
    assert n_pool_f9 == 4740 and n_split_primes == 579

    # [11] Live ACTUAL integer BOUND_CELL (+6 ball-mul margin) check on 12 exact-DP sampled primitive cells
    rng = random.Random(20261122 + seed)
    dp_rd, dp_pref = _build_dp_ball(32)
    q_bc, s_bc, t_bc = _sample_ideal(90, "g", rng)
    B_bc, det_e, det_qcm = _lll_16d_qcm(_sq_hnf(q_bc, s_bc, t_bc))
    assert det_e == q_bc * q_bc and det_qcm == (1265625 ** 2) * q_bc * q_bc
    bc_gmin = 1e9
    bc_gmax = -1e9
    for _ in range(12):
        k = _sample_prim_ball_dp(dp_rd, dp_pref, 32, rng)
        v = [sum(k[i] * B_bc[i][j] for i in range(16)) for j in range(16)]
        nf, ng = _norms_fg(v[:8], v[8:])
        assert ng % q_bc == 0 and nf > 0 and ng > 0
        nf_ub, ng_ub, nf_c, ng_c = BC.bound_cell_int_fg(v[:8], v[8:])
        gf = nf_ub - log2(nf)
        gg = ng_ub - log2(ng)
        assert gf >= 0.0 and gg >= 0.0, f"BOUND_CELL undershoot! gf={gf}, gg={gg}"
        assert nf_c >= ceil(log2(nf)) and ng_c >= ceil(log2(ng))
        bc_gmin = min(bc_gmin, gf, gg)
        bc_gmax = max(bc_gmax, gf, gg)
    assert bc_gmin >= 0.0 and bc_gmax <= 0.5343

    # 100%-fresh production-HNF 16-cal / 48-holdout boot archive with full 64-entry per-target array
    br2 = c7e["boot_r2_32x48"]
    bnp = st["boot_norm_pool"]
    assert br2["fresh_from_production_hnf_all_64"] is True and bnp["fresh_from_production_hnf_all_64"] is True
    cg_by_t = br2["cg_max_by_target"]
    assert len(cg_by_t) == 64 and cg_by_t == bnp["cg_max_by_target"]
    assert abs(br2["cal16_cg_max"] - max(cg_by_t[:16])) < 1e-9 and abs(br2["holdout48_cg_max"] - max(cg_by_t[16:])) < 1e-9
    assert abs(br2["cg_max"] - max(br2["cal16_cg_max"], br2["holdout48_cg_max"])) < 1e-9
    assert br2["n_exceed_m_max"] == 0 and br2["exhaustive_cells_total"] == 64 * 17072 and br2["holdout48_cells"] == 48 * 17072
    assert br2["holdout48_exceed_targets"] == 0 and br2["holdout48_exceed_cal_cut"] == 0 and br2["holdout48_exceed_boot_bins_max"] == 0
    cp95_h48_live = 0.05 ** (1.0 / 48.0)
    assert cp95_h48_live >= Pm.BOOT_KEEP_FRAC and abs(br2["cp95_lower_bound_holdout48"] - round(cp95_h48_live, 5)) < 1e-5
    assert len(br2["raw_deciles_10"]) == 10 and br2["cg_max"] <= max(x[1] for x in Pm.BOOT_BINS)
    assert br2["E_S2_over_ES2"] <= 2.0 and Pm.UDUP <= br2["udup_48cells_mean_inv_m"] <= 1.0

    # Live multi-bin UDUP=0.40 exact-DP recomputation on PLAN16
    _ud_live = {}
    for _uk, _side_u in (
        ("sq_104f_R48_48cells", "f"), ("sq_104g_R48_48cells", "g"),
        ("sq_134f_R56_48cells", "f"), ("sq_134g_R56_48cells", "g"),
        ("sq_204f_R12_48cells", "f"), ("sq_204g_R40_48cells", "g"),
    ):
        _ue = st["unit_orbit_udup"][_uk]
        assert len(_ue["m_true"]) == 48 and len(_ue["cfs"]) == 48 and len(_ue["cgs"]) == 48
        _inv_ms = [1.0 / m for m in _ue["m_true"]]
        _m_inv = sum(_inv_ms) / 48.0
        _se_m = math.sqrt(sum((x - _m_inv) ** 2 for x in _inv_ms) / (47.0 * 48.0))
        _b_lo = float(_ue["bin_lo"]); _b_hi = _b_lo + 5.0; _R_u = int(_ue["R"])
        assert f"_R{_R_u}_" in _uk and abs(_m_inv - _ue["mean_inv_m"]) < 1e-3 and _m_inv - 2.0 * _se_m >= Pm.UDUP
        _SR_u = log2(DS.CELLS[_R_u]); _lq_u = log2(int(_ue["q"]))
        _c_u = (_lq_u + Pm.BP_HI + _SR_u) / 16.0
        _c_dp = (_b_hi + Pm.BP_HI + _SR_u) / 16.0
        _sf_u = (_lq_u if _side_u == "f" else 0.0) + Pm.BP_LO; _sg_u = (_lq_u if _side_u == "g" else 0.0)
        _sf_dp = (_b_lo if _side_u == "f" else 0.0) + Pm.BP_LO; _sg_dp = (_b_lo if _side_u == "g" else 0.0)
        _n_s, _iA_s, _kA_s, _kB_s, _R_plan = DS.PLAN16[(_b_lo, _side_u)]
        _ws = [
            DS.ss_lower2(72.0 * _c_u + cf_ - _sf_u, Pm.LX, float(_n_s), _iA_s, _kA_s, _kB_s, 384, Pm.ALPHA_F, 0.0)
            * DS.ss_lower2(24.0 * _c_u + Pm.LOG2U8 + cg_ - _sg_u, Pm.LX, float(_n_s), _iA_s, _kA_s, _kB_s, 384, Pm.ALPHA_G, 0.0)
            for cf_, cg_ in zip(_ue["cfs"], _ue["cgs"])
        ]
        _pr_c = sum(
            w * DS.ss_lower2(72.0 * _c_dp + cf_ - _sf_dp, Pm.LX, float(_n_s), _iA_s, _kA_s, _kB_s, 384, Pm.ALPHA_F, 0.0)
              * DS.ss_lower2(24.0 * _c_dp + Pm.LOG2U8 + cg_ - _sg_dp, Pm.LX, float(_n_s), _iA_s, _kA_s, _kB_s, 384, Pm.ALPHA_G, 0.0)
            for w, cf_, cg_ in Pm.QBINS
        )
        _pr_sub = sum(
            DS.ss_lower2(72.0 * _c_dp + cf_ - _sf_dp, Pm.LX, float(_n_s), _iA_s, _kA_s, _kB_s, 384, Pm.ALPHA_F, 0.0)
            * DS.ss_lower2(24.0 * _c_dp + Pm.LOG2U8 + cg_ - _sg_dp, Pm.LX, float(_n_s), _iA_s, _kA_s, _kB_s, 384, Pm.ALPHA_G, 0.0)
            for cf_, cg_ in Pm.QBINS_SUB
        ) / len(Pm.QBINS_SUB)
        _pr_dp = min(_pr_sub, Pm.PRCAP_NODE * _pr_c)
        if _R_u == _R_plan:
            assert abs(log2(_pr_dp) - res["DET"][(_b_lo, _b_hi), _side_u]["pr_log2"]) < 1e-12
        _tw = sum(w / m for w, m in zip(_ws, _ue["m_true"])) / sum(_ws)
        _dy = (sum(w / m for w, m in zip(_ws, _ue["m_true"])) / 48.0) / _pr_dp
        assert Pm.UDUP <= _m_inv <= 1.0 and 0.0 < _tw <= 1.0 and _dy >= Pm.UDUP
        _ud_live[_uk] = dict(bin_lo=int(_b_lo), R=_R_u, mean_inv_m=round(_m_inv, 4), true_yield_weighted_inv_m=round(_tw, 4), dedup_yield_over_charged_pr=round(_dy, 4))
    _ub = st["unit_orbit_udup"]["boot_R2_48cells"]
    assert len(_ub["m_true"]) == 48 and len(_ub["cgs"]) == 48
    _cb2 = (Pm.LOG2P4 + log2(DS.CELLS[2])) / 16.0
    _n1b, _iAb, _kAb, _ = DS.BOOT16
    _wsb = [DS.ss_lower2(24.0 * _cb2 + Pm.LOG2U8 + cg_, Pm.LX, float(_n1b), _iAb, _kAb, 1, 384, Pm.ALPHA_G, 0.0) for cg_ in _ub["cgs"]]
    _prb = sum(w * DS.ss_lower2(24.0 * _cb2 + Pm.LOG2U8 + cg_, Pm.LX, float(_n1b), _iAb, _kAb, 1, 384, Pm.ALPHA_G, 0.0) for w, cg_ in Pm.BOOT_BINS)
    _prb = min(sum(DS.ss_lower2(24.0 * _cb2 + Pm.LOG2U8 + cg_, Pm.LX, float(_n1b), _iAb, _kAb, 1, 384, Pm.ALPHA_G, 0.0) for cg_ in Pm.BOOT_SUB) / len(Pm.BOOT_SUB), Pm.PRCAP_BOOT * _prb)
    _inv_mb = [1.0 / m for m in _ub["m_true"]]
    _mb_inv = sum(_inv_mb) / 48.0
    _se_mb = math.sqrt(sum((x - _mb_inv) ** 2 for x in _inv_mb) / (47.0 * 48.0))
    _twb = sum(w / m for w, m in zip(_wsb, _ub["m_true"])) / sum(_wsb)
    _dyb = (sum(w / m for w, m in zip(_wsb, _ub["m_true"])) / 48.0) / _prb
    assert abs(_mb_inv - _ub["mean_inv_m"]) < 1e-3
    assert Pm.UDUP <= _mb_inv - 2.0 * _se_mb and _mb_inv <= 1.0 and Pm.UDUP <= _twb <= 1.0 and _dyb >= Pm.UDUP
    _ud_live["boot_R2_48cells"] = dict(mean_inv_m=round(_mb_inv, 4), true_yield_weighted_inv_m=round(_twb, 4), dedup_yield_over_charged_pr=round(_dyb, 4))

    # Live deterministic rejection test on production-bin (135, g) short-subspace counterexample
    _q_ce = 81009504249198611324243595142290546015061
    _s_ce = 14886707496769063264106121185902690387739
    _t_ce = (5 * pow(8, -1, _q_ce)) % _q_ce
    _B_ce, _, _ = _lll_16d_qcm(_sq_hnf(_q_ce, _s_ce, _t_ce))
    assert BC.screen_special_q_basis(_B_ce, _q_ce, "g", _norms_fg) is False, "Must reject q^4 short-subspace ideal!"
    _k_ce = [1 if i == 0 else (47 if i == 8 else 0) for i in range(16)]
    _v_ce = [sum(_k_ce[i] * _B_ce[i][j] for i in range(16)) for j in range(16)]
    _d135 = res["DET"][(135.0, 140.0), "g"]
    _ok_ce, _nf_ce_rem, _ng_ce_rem, _, _ = BC.bound_cell_check_and_cap(_v_ce[:8], _v_ce[8:], _d135["nf"], _d135["ng"], lq=log2(_q_ce), side="g")
    assert _ok_ce is False and _nf_ce_rem > _d135["nf"] and _ng_ce_rem > _d135["ng"]
    assert BC.screen_special_q_basis(B_bc, q_bc, "g", _norms_fg) is True

    # [11b] Live q15_retain16 prime-splitting + non-vacuous p=2 (K=140) check + counted_fraction
    tab = st["prime_structure"]["table"]
    for pp in (2, 3, 5, 19, 31, 61, 101, 257):
        assert RT.structure(pp) == tab[str(pp)]
    assert abs(RT.counted_fraction() - st["omitted_law"]["counted_fraction"]) < 1e-6
    assert st["omitted_law"]["tail_j_ge_74_included"] is True
    rt_v = RT.verify_retain16(n_samples=64, seed=11, recompute_full=False)
    assert rt_v["probes_verified"] >= 10 and rt_v["sample_node_slack96"]["loss_f"] < 5e-5 and rt_v["sample_node_slack96"]["loss_g"] < 5e-5
    assert rt_v["residue_checks"][2]["K"] == 140 and rt_v["residue_checks"][2]["max_g"] <= 128.0

    # [12] Live q15_ecm16 L_min=2 branch-tree DP + geometric ladder + Suyama/Newton/singular-N-1 check
    assert EC.L_MIN_DEFAULT == 2 and EC.B_MIN_DEFAULT == 32 and EC.SEPARATION_DEFAULT == 0.5 and EC.RETRY_PR_DESC == 0.0043 and EC.NEWTON_PRECHECK_INSTR == 2.0 ** 26
    l2_boot, _ = EC.tree_branch_caps(int(ceil(bt["m"])), EC.B_MIN_DEFAULT, EC.L_MIN_DEFAULT)
    assert l2_boot >= ceil(bt["m"] / 64.0) ** 2
    sch_node = EC.schedule_v9(Pm.LX, int(ceil(res["DET"][(68.5, 75.0), "f"]["nf"])), Pm.ECM_SUCCESS, n_factors=50)
    sch_boot = EC.schedule_v9(bt["nA"], int(ceil(bt["m"])), Pm.ECM_SUCCESS, n_factors=bt["boot_n_factors"])
    lad_boot = EC.ladder_v9_instr(16.0, bt["nA"], int(ceil(bt["m"])), Pm.ECM_SUCCESS, bt["boot_n_factors"])
    assert sch_node["curves"] > 0 and sch_boot["curves"] > 0 and lad_boot >= sch_boot["instructions"]
    suyama_retry = EC.verify_suyama_retry_counterexample()
    _suyama_ok = bool(suyama_retry["perfect_power_precheck_ok"] and suyama_retry["bad_gcd_eq_N"] and suyama_retry["singular_nm1_rejected"] and suyama_retry["split_factor"] > 1)
    assert _suyama_ok, "Suyama retry & perfect-power/singular precheck failed!"

    wall_s = time.perf_counter() - t0
    return dict(
        all_ok=True,
        plan_tag="CERT1_p1416_lx685_lad12",
        wall_s=round(wall_s, 3),
        eval_384_s=round(dt_ev, 3),
        mesh=mesh,
        n_nodes=len(rows),
        n_sieve=len(_det_sieve),
        n_ecm=len(_det_ecm),
        lx=Pm.LX,
        bp=[Pm.BP_LO, Pm.BP_HI],
        k_trials=Pm.K_TRIALS,
        p_node=Pm.P_NODE,
        two_log_units=round(two_log, 6),
        boot_u=round(bt["boot_u"], 4),
        sub_u=round(bt["sub_u"], 4),
        phi_tree=round(bt["phi_tree"], 6),
        e_runs=round(bt["e_runs"], 4),
        n_slots_tree=int(bt["n_slots_tree"]),
        expected_slots_2log=round(EN, 2),
        boot_phi=0.0,
        boot_m=round(bt["m"], 2),
        boot_miss=round(bt["miss"], 5),
        boot_lin_keep=round(bt["lin_keep"], 5),
        max_miss_node=round(max_miss_node, 5),
        boot_keep_frac=bt["boot_keep_frac"],
        boot_n_factors=bt["boot_n_factors"],
        boot_n1=int(bt["n1"]),
        boot_nA=round(bt["nA"], 2),
        boot_kA=int(bt["kA"]),
        boot_R=int(bt["R"]),
        boot_one_log_log2=round(bt["one_log_log2"], 4),
        cert=dict(
            mu_cert=Pm.MU_CERT, delta=Pm.DELTA_CERT, s_cert=QC.S_CERT, n_first_draws=200 * 76 + 400 * len(bent), failures=fails,
            p_bad_entry_cp95=round(p_bad_entry, 5), p_bad_boot_cp95=round(p_bad_boot, 5), min_cert=cert_min,
            min_cert_over_plan_mu=round(ratios[0], 3), median_cert_over_plan_mu=round(ratios[len(ratios) // 2], 3),
            min_muhat_over_plan_mu=round(min(mratios), 3), live=live_cert, p_empty=Pm.P_EMPTY, p_low=Pm.P_LOW, toy=toy,
        ),
        pool_f9_p1416=dict(n_ideals=n_pool_f9, n_split_primes=n_split_primes),
        pm_log2_rels=Pm.LOG2_RELS,
        pm_main_lb=Pm.MAIN_LB,
        pm_main_nf=round(Pm.MAIN_UF * Pm.MAIN_LB, 4),
        pm_main_ng=round(Pm.MAIN_UG * Pm.MAIN_LB, 4),
        pm_alpha_f=Pm.ALPHA_F,
        pm_alpha_g=Pm.ALPHA_G,
        pm_pair_setup=Pm.PAIR_SETUP,
        pm_lll_node=Pm.LLL_NODE,
        pm_lll_boot=Pm.LLL_BOOT,
        pm_bound_cell=Pm.BOUND_CELL,
        pm_sieve_hit=Pm.SIEVE_HIT,
        pm_ret_mult_f=Pm.RET_MULT_F,
        pm_ret_mult_g=Pm.RET_MULT_G,
        f4_counterexample_rejected=True,
        suyama_retry_ok=_suyama_ok,
        udup_live_recomputed=_ud_live,
        eager_all_100pct=True,
        linear_miss_retention=True,
        min_mu=round(min_mu, 2),
        min_ret=round(min_ret, 5),
        max_loss=float(f"{max_loss:.3e}"),
        bound_cell_live_min_gap_bits=bc_gmin,
        bound_cell_live_max_gap_bits=bc_gmax,
        retain16_live_ok=True,
        retain16_probes=int(rt_v["probes_verified"]),
        retain16_p2_K=int(rt_v["residue_checks"][2]["K"]),
        retain16_p2_max_g=rt_v["residue_checks"][2]["max_g"],
        boot_r2=br2,
        ecm_schedules=dict(
            l_min=EC.L_MIN_DEFAULT,
            b_min=EC.B_MIN_DEFAULT,
            sep=EC.SEPARATION_DEFAULT,
            p_r=EC.RETRY_PR_DESC,
            node_685f=dict(curves=sch_node["curves"], n_factors=sch_node["n_factors"], log2_ir=round(sch_node["log2_instructions"], 2)),
            boot_nA=dict(curves=sch_boot["curves"], n_factors=sch_boot["n_factors"], log2_ir=round(sch_boot["log2_instructions"], 2), ladder_log2_ir=round(log2(lad_boot), 2)),
        ),
    )


if __name__ == "__main__":
    out = run_all_checks()
    print(json.dumps(out, indent=2))
