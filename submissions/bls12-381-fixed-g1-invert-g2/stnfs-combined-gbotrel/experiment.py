#!/usr/bin/env python3
"""Entrypoint experiment script for the combined STNFS submission: CADO
filtering, 11M/bit Suyama Montgomery ECM validation, |u|^(1/2)-balanced 12D
boot lattice & both-side special-q witnesses, toy fixed-ball descent, and
Markov-capped two-window descent recurrence."""

import bisect
import cmath
import functools
import heapq
import json
import math
import os
import random
import sys
import time
from decimal import Decimal as D, localcontext
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stnfs_core import (
    F_TERMS,
    G_TERMS,
    H_POLY,
    LN2,
    LOG2U6,
    MONT5,
    P_MOD,
    R_ORD,
    RHO_TAB,
    RSTEP,
    U_ABS,
    UNIT,
    certify_ellipsoid_counts,
    compute_prime_structure_and_alphas,
    evaluate_main_phase_live,
    exact_norms,
    kappa_rho,
    log2_li,
    p_div,
    p_gcd,
    p_powmod,
    p_sub,
    p_trim,
    primes_range,
    rho_hi,
    rho_lo,
    structure_of_prime,
    verify_tower_isomorphism_and_fields,
    zy_mul,
    zy_norm,
)


def build_filt_system(L, us, excess, seed, log2U, t_dense=0.8, nt=2000):
    rng = random.Random(seed)
    LiB = (2.0**L) / (L * LN2)
    Rprod = excess * LiB
    Up = 2**log2U
    Rp = int(round(excess * Up))
    lsmall = max(8.0, math.log2(Rp / 25.0))
    ts = [t_dense + (1.0 - t_dense) * (i + 0.5) / nt for i in range(nt)]
    dens = [2.0 ** (L * (t - 1.0)) / t for t in ts]
    cum, tot = [], 0.0
    for d in dens:
        tot += d
        cum.append(tot)
    cum = [c / tot for c in cum]
    cols, w = [], []
    for s_, u in enumerate(us):
        for _ in range(Up // 2):
            i = bisect.bisect_left(cum, rng.random())
            t = ts[min(i, nt - 1)] + (rng.random() - 0.5) * (1.0 - t_dense) / nt
            lam = Rprod * 2.0 * kappa_rho(u, t) / (2.0 ** (L * t))
            cols.append((s_, t))
            w.append(min(1.0, lam / Rp))
    small = []
    for s_, u in enumerate(us):
        for q in primes_range(2, int(2**lsmall) + 1):
            t = math.log(q) / (L * LN2)
            pr = min(1.0, 2.0 * kappa_rho(u, t) / q)
            cols.append((s_, t))
            w.append(pr)
            small.append((len(cols) - 1, pr, q))
    tmin = lsmall / L
    w_mid = sum(
        kappa_rho(u, tmin + (t_dense - tmin) * (i + 0.5) / 1000)
        * ((t_dense - tmin) / 1000)
        / (tmin + (t_dense - tmin) * (i + 0.5) / 1000)
        for u in us
        for i in range(1000)
    )
    nsamp = (Up // 2) * len(us)
    dense = [i for i in range(nsamp) if w[i] >= 0.2]
    sparse = [i for i in range(nsamp) if w[i] < 0.2]
    scum, stot = [], 0.0
    for i in sparse:
        stot += w[i]
        scum.append(stot)

    def r_pois(lam):
        if lam < 30.0:
            lim = math.exp(-lam)
            k, p = 0, 1.0
            while True:
                p *= rng.random()
                if p <= lim:
                    return k
                k += 1
        return max(0, int(round(rng.gauss(lam, math.sqrt(lam)))))

    rows = []
    for _ in range(Rp):
        row = {}
        for i in dense:
            if rng.random() < w[i]:
                row[i] = 1
        for _ in range(r_pois(stot)):
            j = sparse[bisect.bisect_left(scum, rng.random() * stot)]
            row[j] = row.get(j, 0) + 1
        for j, pr, q in small:
            if rng.random() < pr:
                v = 1
                while rng.random() < 1.0 / q:
                    v += 1
                row[j] = v
        row[-1] = w_mid
        rows.append(row)
    return cols, w, rows


def run_filtering_and_missing_audit(seed=1, log2U=12):
    cols, w, rows = build_filt_system(69.0, (750.0 / 69.0, 498.5 / 69.0), 1.3, seed, log2U)
    U = len(cols)
    Usamp = 2 * (2**log2U // 2)
    Rp = len(rows)
    occ = [0] * U
    for r in rows:
        for c in r:
            if c >= 0:
                occ[c] += 1

    colw = occ[:]
    alive = [True] * len(rows)
    nalive = len(rows)

    def purge_sing():
        nonlocal nalive
        ch = True
        while ch:
            ch = False
            for i, r in enumerate(rows):
                if alive[i] and any(c >= 0 and colw[c] == 1 for c in r):
                    alive[i] = False
                    nalive -= 1
                    ch = True
                    for cc in r:
                        if cc >= 0:
                            colw[cc] -= 1

    purge_sing()
    ncols = sum(1 for x in colw if x > 0)
    target_exc = max(2, int(0.002 * Usamp))
    while nalive - ncols > target_exc:
        parent = list(range(len(rows)))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        c_rows = {}
        for i, r in enumerate(rows):
            if alive[i]:
                for c in r:
                    if c >= 0 and colw[c] == 2:
                        c_rows.setdefault(c, []).append(i)
        for rs in c_rows.values():
            if len(rs) == 2:
                a, b = find(rs[0]), find(rs[1])
                if a != b:
                    parent[a] = b
        comp = {}
        for i in range(len(rows)):
            if alive[i]:
                comp.setdefault(find(i), []).append(i)
        comps = sorted(comp.values(), key=lambda rs: -sum((len(rows[i]) - 1) + rows[i].get(-1, 0.0) for i in rs))
        need = nalive - ncols - target_exc
        nb = max(1, min(need // 2 + 1, max(1, len(comps) // 100)))
        for rs in comps[:nb]:
            for i in rs:
                if alive[i]:
                    alive[i] = False
                    nalive -= 1
                    for c in rows[i]:
                        if c >= 0:
                            colw[c] -= 1
        purge_sing()
        ncols = sum(1 for x in colw if x > 0)

    kept = [dict(rows[i]) for i in range(len(rows)) if alive[i]]
    known = [False] * U
    for r in kept:
        for c in r:
            if c >= 0:
                known[c] = True
    ch = True
    while ch:
        ch = False
        for r in rows:
            unk = [c for c in r if c >= 0 and not known[c]]
            if len(unk) == 1:
                known[unk[0]] = True
                ch = True

    bands = []
    for lo, hi in ((1, 2), (2, 3), (3, 5), (5, 8), (8, 1e9)):
        idx = [i for i in range(Usamp) if lo <= w[i] * Rp < hi and occ[i] > 0]
        if idx:
            meas = sum(1 for i in idx if not known[i]) / len(idx)
            mod = sum(min(1.0, math.exp(-w[i] * Rp) * (1.0 + w[i] * Rp + 0.5 * ((w[i] * Rp) ** 2))) for i in idx) / len(idx)
            assert meas <= mod
            bands.append(dict(band=f"[{lo},{hi})", ideals=len(idx), measured_missing=round(meas, 4), model_poisson_le2=round(mod, 4)))

    m_rows = [dict(r) for r in kept]
    colrows = {}
    for i, r in enumerate(m_rows):
        for c in r:
            if c >= 0:
                colrows.setdefault(c, set()).add(i)
    nrows = len(m_rows)
    totw = sum((len(r) - 1) + r.get(-1, 0.0) for r in m_rows)
    heap = [(len(s), c) for c, s in colrows.items()]
    heapq.heapify(heap)
    while heap and totw / max(nrows, 1) < 100.0:
        wgt, c = heapq.heappop(heap)
        if c not in colrows or len(colrows[c]) != wgt:
            if c in colrows:
                heapq.heappush(heap, (len(colrows[c]), c))
            continue
        s = colrows[c]
        if wgt == 0:
            del colrows[c]
            continue
        if wgt > 64:
            break
        if wgt == 1:
            i = next(iter(s))
            nrows -= 1
            totw -= (len(m_rows[i]) - 1) + m_rows[i].get(-1, 0.0)
            for cc in m_rows[i]:
                if cc >= 0:
                    colrows[cc].discard(i)
                    if colrows[cc]:
                        heapq.heappush(heap, (len(colrows[cc]), cc))
                    else:
                        del colrows[cc]
            m_rows[i] = {}
            continue
        pivs = [i for i in s if abs(m_rows[i][c]) == 1]
        if not pivs:
            continue
        p = min(pivs, key=lambda i: (len(m_rows[i]) - 1) + m_rows[i].get(-1, 0.0))
        pr, pc = m_rows[p], m_rows[p][c]
        if any(abs(m_rows[i].get(cc, 0) - m_rows[i][c] * pc * v) > 255 for i in s if i != p for cc, v in pr.items() if cc >= 0):
            continue
        for i in list(s):
            if i == p:
                continue
            r = m_rows[i]
            m = r[c] * pc
            r[-1] = r.get(-1, 0.0) + pr.get(-1, 0.0)
            totw += pr.get(-1, 0.0)
            for cc, v in pr.items():
                if cc < 0:
                    continue
                nv = r.get(cc, 0) - m * v
                if nv == 0:
                    if cc in r:
                        del r[cc]
                        colrows[cc].discard(i)
                        totw -= 1
                        if colrows[cc]:
                            heapq.heappush(heap, (len(colrows[cc]), cc))
                        else:
                            del colrows[cc]
                else:
                    if cc not in r:
                        colrows[cc].add(i)
                        totw += 1
                        heapq.heappush(heap, (len(colrows[cc]), cc))
                    r[cc] = nv
        nrows -= 1
        totw -= (len(pr) - 1) + pr.get(-1, 0.0)
        for cc in pr:
            if cc != c and cc >= 0:
                colrows[cc].discard(p)
                if colrows[cc]:
                    heapq.heappush(heap, (len(colrows[cc]), cc))
                else:
                    del colrows[cc]
        m_rows[p] = {}
        del colrows[c]

    ratio_2_12 = nrows / Usamp
    ratio_2_17 = 1.0 / 18.600
    pad_n = 1.60
    assert ratio_2_17 * pad_n > ratio_2_12
    log2_U_prod = math.log2(2.0**log2_li(69.0) + 2.0**31)
    log2_N = log2_U_prod + math.log2(ratio_2_17 * pad_n)
    log2_la_units = math.log2(3.0) + 2.0 * log2_N + math.log2(100.0 * 96.0 / UNIT)

    return {
        "scale_2_12_N_over_U": round(ratio_2_12, 4),
        "scale_2_17_N_over_U": round(ratio_2_17, 4),
        "pad_N": pad_n,
        "charged_N_over_U": round(ratio_2_17 * pad_n, 4),
        "log2_N_filtered": round(log2_N, 4),
        "log2_la_units": round(log2_la_units, 4),
        "missing_log_bands": bands,
    }


@functools.lru_cache(None)
def ecm_semismooth(fbits, b1, b2, steps=256):
    m = fbits - math.log2(12.0)
    base, top = math.log2(b1), math.log2(b2)
    acc = rho_lo(m / base)
    if top > base:
        for i in range(steps):
            lo = base + (top - base) * i / steps
            hi = base + (top - base) * (i + 1) / steps
            acc += rho_lo((m - lo) / base) * math.log(hi / lo)
    return min(acc, 1.0 - 1e-15)


@functools.lru_cache(None)
def ecm_power_corr(fbits, b1):
    m = fbits - math.log2(12.0)
    u = m / math.log2(b1)
    tilt = rho_hi(max(u - 2.0, 1.0)) / max(rho_lo(u), 1e-300)
    y = math.sqrt(b1)
    return max(0.0, 1.0 - tilt * 2.0 / (y * math.log(y)))


@functools.lru_cache(None)
def ecm_schedule(fbits, nbits, success=0.999):
    nbits = int(math.ceil(nbits))
    L = max(5, math.ceil(nbits / 64))
    mul_i = MONT5 * ((L / 5.0) ** 2)
    gcd_i = 2.0 * nbits * (12 * L + 16)
    best, wheel, J = None, 30030, 2880
    for k in range(10, 37, 2):
        b1 = 2**k
        s1 = math.ceil(1.5 * b1) * 11
        pi_b1 = 1.26 * b1 / math.log(b1)
        for ratio in (1, 10, 30, 100, 300):
            b2 = b1 * ratio
            p1 = ecm_semismooth(fbits, b1, b2) * ecm_power_corr(fbits, b1)
            if p1 <= 1e-12:
                continue
            curves = math.ceil(math.log1p(-success) / math.log1p(-p1))
            G = math.ceil(b2 / wheel) + 1 if ratio > 1 else 0
            s2 = 6 * (wheel // 2 + G) + 3 * J * G if ratio > 1 else 0
            ops = 1.1 * (s1 + s2)
            per_curve = ops * mul_i + 2.0 * gcd_i
            redo = p1 * (ops * mul_i + (pi_b1 + G) * gcd_i)
            cost = curves * (per_curve + redo)
            if best is None or cost < best["instructions"]:
                best = dict(fbits=fbits, nbits=nbits, B1=b1, B2=b2, curves=curves, p_curve=p1, stage1_M=s1, stage2_M=s2, instructions=cost, log2_instructions=round(math.log2(cost), 3))
    return best


def ecm_instr(fbits, nbits, success=0.999):
    return ecm_schedule(fbits, nbits, success)["instructions"]


def xdbl_mod(P_, a24, n):
    x, z = P_
    u = ((x + z) * (x + z)) % n
    v = ((x - z) * (x - z)) % n
    w = (u - v) % n
    return (u * v) % n, (w * (v + a24 * w)) % n


def xadd_mod(P_, Q, diff, n):
    x, z = P_
    y, w = Q
    dx, dz = diff
    u = ((x + z) * (y - w)) % n
    v = ((x - z) * (y + w)) % n
    return (dz * ((u + v) ** 2)) % n, (dx * ((u - v) ** 2)) % n


def ladder_mod(P_, k, a24, n):
    if k == 0:
        return (1, 0)
    R0, R1 = P_, xdbl_mod(P_, a24, n)
    for bit in bin(k)[3:]:
        if bit == "0":
            R1 = xadd_mod(R0, R1, P_, n)
            R0 = xdbl_mod(R0, a24, n)
        else:
            R0 = xadd_mod(R0, R1, P_, n)
            R1 = xdbl_mod(R1, a24, n)
    return R0


def mr_probable_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if a % n == 0:
            continue
        x = pow(a % n, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True


def validate_ecm_stage1_stage2(trials=50, seed=1):
    rng = random.Random(seed)
    p_b1 = primes_range(2, 2049)
    out = []
    for bits, b1, ratio in ((36, 1024, 50), (40, 2048, 50)):
        b2, wheel = b1 * ratio, 210
        s1_hits = s2_hits = valid = 0
        pred_sum = 0.0
        for _ in range(trials):
            p = rng.randrange(1 << (bits - 1), 1 << bits) | 1
            while not mr_probable_prime(p):
                p += 2
            sigma = rng.randrange(6, p)
            u = (sigma * sigma - 5) % p
            v = (4 * sigma) % p
            den = (16 * pow(u, 3, p) * v) % p
            if not den:
                continue
            a24 = (pow(v - u, 3, p) * (3 * u + v) * pow(den, -1, p)) % p
            Q = (pow(u, 3, p), pow(v, 3, p))
            for q in p_b1:
                if q > b1:
                    break
                pw = q
                while pw * q <= b1:
                    pw *= q
                Q = ladder_mod(Q, pw, a24, p)
            valid += 1
            lp = math.log2(p)
            pred_sum += ecm_semismooth(lp, b1, b2) * ecm_power_corr(lp, b1)
            if Q[1] == 0:
                s1_hits += 1
                continue
            baby = {j: ladder_mod(Q, j, a24, p) for j in range(1, wheel // 2 + 1) if math.gcd(j, wheel) == 1}
            step = ladder_mod(Q, wheel, a24, p)
            prev, cur = (1, 0), step
            found2 = False
            for i in range(1, math.ceil(b2 / wheel) + 1):
                x, z = cur
                for j, (y, w) in baby.items():
                    if (b1 < i * wheel - j <= b2 or b1 < i * wheel + j <= b2) and (x * w - y * z) % p == 0:
                        found2 = True
                        break
                if found2:
                    s2_hits += 1
                    break
                nxt = xdbl_mod(step, a24, p) if i == 1 else xadd_mod(cur, step, prev, p)
                prev, cur = cur, nxt
        emp = (s1_hits + s2_hits) / valid
        mod_c = pred_sum / valid
        se = math.sqrt(mod_c * (1.0 - mod_c) / valid)
        assert emp >= mod_c - 2.5 * se
        out.append(dict(bits=bits, B1=b1, B2=b2, trials=valid, stage1=s1_hits, stage2=s2_hits, empirical=round(emp, 4), model_corrected=round(mod_c, 4)))
    return out


CELLS = {
    4: 13315520, 6: 1573359464, 8: 47982017060, 10: 687193400292,
    12: 6073556677752, 14: 38414493101628, 16: 190058478539564, 20: 2754504316057176,
    24: 24505531366261200, 28: 155610653675916932, 32: 771909602591773000, 40: 11221331971811657968,
    48: 99995402971660463340, 56: 635616373364584622852, 64: 3155035383949694198996, 80: 45900097950907999286464,
    96: 409193506988359583551012, 112: 2601677688965738307448548, 128: 12916156661164168796080328,
    160: 187942845067274337246511324, 192: 1675659512501707896389090868, 224: 10654623198171176305341831664,
    256: 52897551248440914796084637156, 288: 217401149762189752879411886484, 320: 769748406668707279529221545272,
    384: 6863095027560050329435265332068,
}
WB_SHIFT = 32
WB_SCALE = 1 << WB_SHIFT


def lll_12d(B_in, prec=260):
    B = [v[:] for v in B_in]
    n = len(B)
    swaps = 0
    with localcontext() as ctx:
        ctx.prec = prec
        mu, norm, star = [], [], []
        for v in B:
            s = list(map(D, v))
            row = []
            for j, b in enumerate(star):
                x = sum(D(a) * c for a, c in zip(v, b)) / norm[j]
                row.append(x)
                s = [a - x * c for a, c in zip(s, b)]
            star.append(s)
            mu.append(row)
            norm.append(sum(x * x for x in s))
        k = 1
        while k < n:
            for j in range(k - 1, -1, -1):
                q = int(mu[k][j].to_integral_value(rounding="ROUND_HALF_EVEN"))
                if q:
                    B[k] = [x - q * y for x, y in zip(B[k], B[j])]
                    for i in range(j):
                        mu[k][i] -= q * mu[j][i]
                    mu[k][j] -= q
            m = mu[k][k - 1]
            if norm[k] >= (D("0.99") - m * m) * norm[k - 1]:
                k += 1
                continue
            A = norm[k] + m * m * norm[k - 1]
            newm = m * norm[k - 1] / A
            norm[k] = norm[k] * norm[k - 1] / A
            norm[k - 1] = A
            mu[k][k - 1] = newm
            B[k], B[k - 1] = B[k - 1], B[k]
            for j in range(k - 1):
                mu[k][j], mu[k - 1][j] = mu[k - 1][j], mu[k][j]
            for i in range(k + 1, n):
                t = mu[i][k]
                mu[i][k] = mu[i][k - 1] - m * t
                mu[i][k - 1] = t + newm * mu[i][k]
            k = max(1, k - 1)
            swaps += 1
    return B, swaps


def gs_exact_12d(B):
    stars, mu, norm = [], [], []
    for v in B:
        row, s = [], list(map(F, v))
        for j, b in enumerate(stars):
            m = sum(x * y for x, y in zip(v, b)) / norm[j]
            row.append(m)
            s = [x - m * y for x, y in zip(s, b)]
        stars.append(s)
        mu.append(row)
        norm.append(sum(x * x for x in s))
    return mu, norm


def ball_point_12d(rng, Rr):
    while True:
        g = [rng.gauss(0, 1) for _ in range(12)]
        nrm = math.sqrt(sum(x * x for x in g))
        r = Rr * (rng.random() ** (1.0 / 12.0))
        c = [int(round(x * r / nrm)) for x in g]
        if any(c) and sum(x * x for x in c) <= Rr * Rr:
            for x in reversed(c):
                if x:
                    if x < 0:
                        c = [-y for y in c]
                    break
            return c


def f6_add(a, b):
    return [(x + y) % P_MOD for x, y in zip(a, b)]


def f6_neg(a):
    return [(-x) % P_MOD for x in a]


def f6_mul(a, b):
    return [x % P_MOD for x in zy_mul(a, b)]


def f6_pow(a, n):
    v = [1, 0, 0, 0, 0, 0]
    while n:
        if n & 1:
            v = f6_mul(v, a)
        a = f6_mul(a, a)
        n >>= 1
    return v


def f12_mul(a, b):
    ac = f6_mul(a[0], b[0])
    bd = f6_mul(a[1], b[1])
    cross = f6_mul(f6_add(a[0], a[1]), f6_add(b[0], b[1]))
    return f6_add(ac, f6_neg(f6_mul([0, U_ABS, 0, 0, 0, 0], bd))), f6_add(cross, f6_neg(f6_add(ac, bd)))


def f12_pow(a, n):
    v = ([1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0])
    while n:
        if n & 1:
            v = f12_mul(v, a)
        a = f12_mul(a, a)
        n >>= 1
    return v


def roots_mod_q(f, q, rng):
    f = p_trim([x % q for x in f])
    f = p_gcd(f, p_sub(p_powmod([0, 1], q, f, q), [0, 1], q), q)

    def split(g):
        if len(g) <= 1:
            return []
        if len(g) == 2:
            return [(-g[0] * pow(g[1], -1, q)) % q]
        while True:
            a = [rng.randrange(q) for _ in range(len(g) - 1)]
            h = p_gcd(g, p_sub(p_powmod(a, (q - 1) // 2, g, q), [1], q), q)
            if 1 < len(h) < len(g):
                q_poly, _ = p_div(g, h, q)
                return split(h) + split(q_poly)

    return sorted(split(f))


def build_specialq_lattice(bits, side, rng):
    while True:
        q = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if not mr_probable_prime(q):
            continue
        ss = roots_mod_q(H_POLY, q, rng)
        opts = []
        for s in ss:
            ff = [U_ABS * s, 0, 1] if side == "g" else [1 + s - s * s - s**3, 0, s**5, 0, 0, 0, 2 * (s**3), 0, 0, 0, -2 * s, 0, 1]
            opts.extend((s, t) for t in roots_mod_q(ff, q, rng))
        if opts:
            break
    s, t = rng.choice(opts)
    w = [pow(s, j, q) for j in range(6)] + [(-t * pow(s, j, q)) % q for j in range(6)]
    B0 = [[q] + [0] * 11] + [[(-w[j]) % q] + [int(k == j) for k in range(1, 12)] for j in range(1, 12)]
    B, swaps = lll_12d(B0)
    assert all(sum(a * b for a, b in zip(v, w)) % q == 0 for v in B)
    _, nn = gs_exact_12d(B)
    assert math.prod(nn) == q * q
    return q, s, t, B, swaps


@functools.lru_cache(None)
def complex_forms():
    def croots(c):
        n = len(c) - 1
        c = [complex(x) / c[-1] for x in c]
        rad = 1.0 + max(abs(x) for x in c[:-1])
        z = [rad * cmath.exp(2j * math.pi * (k + 0.25) / n) for k in range(n)]
        for _ in range(2500):
            nz = [z[i] - sum(c[m] * (z[i] ** m) for m in range(n + 1)) / math.prod(z[i] - z[j] for j in range(n) if j != i) for i in range(n)]
            if max(abs(a - b) for a, b in zip(z, nz)) < 1e-15 * rad:
                return nz
            z = nz
        return z

    F_, G_ = [], []
    for e in croots(H_POLY):
        F_ += [(e, x) for x in croots([1 + e - e * e - e**3, 0, e**5, 0, 0, 0, 2 * (e**3), 0, 0, 0, -2 * e, 0, 1])]
        x = cmath.sqrt(-U_ABS * e)
        G_ += [(e, x), (e, -x)]
    return F_, G_


def cell_bound_forms(v, fs):
    t = 0.0
    for e, x in fs:
        p = [e**j for j in range(6)]
        L = abs(sum(v[j] * p[j] for j in range(6)) - x * sum(v[6 + j] * p[j] for j in range(6)))
        t += math.log2(L + (2.0**-46) * (sum(abs(v[j]) * abs(p[j]) for j in range(6)) + abs(x) * sum(abs(v[6 + j]) * abs(p[j]) for j in range(6))))
    return t


def audit_balanced_boot_and_both_side_specialq(seed=31001):
    rng = random.Random(seed)
    x = ([rng.randrange(P_MOD) for _ in range(6)], [rng.randrange(1, P_MOD) for _ in range(6)])
    z = f12_pow(f12_pow(x, (P_MOD**12 - 1) // R_ORD), rng.randrange(1, R_ORD))
    ratio = f6_mul(f6_neg(z[0]), f6_pow(z[1], P_MOD**6 - 2))
    raw_basis = [[P_MOD if i == j else 0 for j in range(12)] for i in range(6)] + [
        f6_mul(ratio, [0] * j + [1]) + [int(i == j) for i in range(6)] for j in range(6)
    ]

    scaled_in = [row[:6] + [c * WB_SCALE for c in row[6:]] for row in raw_basis]
    B_scaled, swaps_boot = lll_12d(scaled_in)
    mu_b, norm_b = gs_exact_12d(B_scaled)
    for i in range(12):
        assert all(abs(m) <= F(1, 2) for m in mu_b[i])
        if i:
            assert norm_b[i] >= (F(98, 100) - mu_b[i][-1] ** 2) * norm_b[i - 1]
    assert math.prod(norm_b) == (P_MOD**12) * (WB_SCALE**12)

    B_bal = [row[:6] + [c // WB_SCALE for c in row[6:]] for row in B_scaled]
    for v in B_bal:
        assert f6_mul(v[:6], z[1]) == f6_neg(f6_mul(v[6:], z[0]))

    HALF_LOG2U6 = 0.5 * LOG2U6
    boot_stats = []
    for Rr in (4, 96):
        S = math.log2(CELLS[Rr])
        cb = (6.0 * math.log2(P_MOD) + S) / 12.0
        cgs = []
        for _ in range(12):
            k = ball_point_12d(rng, Rr)
            v = [sum(k[i] * B_bal[i][j] for i in range(12)) for j in range(12)]
            assert f6_mul(v[:6], z[1]) == f6_neg(f6_mul(v[6:], z[0]))
            ng = zy_norm([a + b for a, b in zip(zy_mul(v[:6], v[:6]), zy_mul([0, U_ABS], zy_mul(v[6:], v[6:])))])
            assert ng > 0
            cgs.append(math.log2(ng) - 12.0 * cb - HALF_LOG2U6)
        boot_stats.append(dict(R=Rr, cg_bal_mean=round(sum(cgs) / len(cgs), 2), cg_bal_max=round(max(cgs), 2)))

    sq_audit = []
    FF, GG = complex_forms()
    for bits, side, s_seed in ((150, "g", 12201), (200, "f", 13300)):
        s_rng = random.Random(s_seed)
        q, s, t, B_sq, sw = build_specialq_lattice(bits, side, s_rng)
        lq = math.log2(q)
        Rr = 96
        S = math.log2(CELLS[Rr])
        c = (lq + S) / 12.0
        cfs, cgs = [], []
        max_gap = 0.0
        for _ in range(10):
            k = ball_point_12d(s_rng, Rr)
            while math.gcd(*k) > 1:
                k = ball_point_12d(s_rng, Rr)
            v = [sum(k[i] * B_sq[i][j] for i in range(12)) for j in range(12)]
            nf, ng = exact_norms(v[:6], v[6:])
            assert (nf if side == "f" else ng) % q == 0
            cfs.append(math.log2(nf) - 72.0 * c)
            cgs.append(math.log2(ng) - 12.0 * c - LOG2U6)
            for N_val, fs in ((nf, FF), (ng, GG)):
                gap = cell_bound_forms(v, fs) - math.log2(N_val)
                assert gap >= -1e-9
                max_gap = max(max_gap, gap)
        sq_audit.append(dict(bits=bits, side=side, swaps=sw, cf_max=round(max(cfs), 2), cg_max=round(max(cgs), 2), certified_bound_max_gap_bits=round(max_gap, 12)))

    return {
        "balanced_boot_wb_shift": WB_SHIFT,
        "boot_norm_reduction_bits": round(HALF_LOG2U6, 2),
        "boot_swaps": swaps_boot,
        "boot_ball_stats": boot_stats,
        "specialq_both_sides": sq_audit,
    }


GRID_RET = 4


def _convolve_pdf(a, b):
    out = {}
    for x, px in a.items():
        for y, py in b.items():
            out[x + y] = out.get(x + y, 0.0) + px * py
    return out


def _prune_pdf(pdf):
    keep = {x: v for x, v in pdf.items() if v > 1e-60}
    lost = sum(v for x, v in pdf.items() if x not in keep)
    if lost:
        top = max(pdf)
        keep[top] = keep.get(top, 0.0) + lost
    return keep


def local_omitted_pdf(p, st, side, sieved, up):
    lp = math.log2(p)
    rnd = math.ceil if up else math.floor
    pdf = {0: 1.0}
    for d1, e, mf, mg in st:
        mults = mf if side == "f" else mg
        if not mults or (sieved and d1 != 1):
            continue
        q = p**d1
        tail = lambda j: (q / (q + 1.0)) * sum(min(1.0 / q, mu * (q ** (-((j + mu - 1) // mu)))) for mu in mults)
        comp = {}
        prev = tail(1)
        comp[0] = 1.0 - prev
        for j in range(1, 1000):
            nxt = tail(j + 1)
            m_u = rnd(d1 * (j - 1 if sieved else j) * lp * GRID_RET)
            comp[m_u] = comp.get(m_u, 0.0) + (prev - nxt)
            prev = nxt
            if nxt < 1e-60:
                m_u2 = rnd(d1 * (j + (0 if sieved else 1)) * lp * GRID_RET)
                comp[m_u2] = comp.get(m_u2, 0.0) + nxt
                break
        pdf = _convolve_pdf(pdf, comp)
    return pdf


@functools.lru_cache(None)
def base_primitive_fraction(pmax_exact=512):
    f = 1.0
    for p in primes_range(2, pmax_exact):
        for d1, e, mf, mg in structure_of_prime(p):
            f *= 1.0 - (p ** (-2.0 * d1))
    for j in range(int(math.log2(pmax_exact)), 200):
        f *= 1.0 - 6.0 * 1.26 * (2.0 ** (j + 1)) / ((j + 1) * LN2) * (2.0 ** (-2 * j))
    return f


@functools.lru_cache(None)
def omitted_mass_pdf(side, up=True, pmax_exact=512):
    pdf = {0: 1.0}
    ordinary = 1.0
    for p in primes_range(2, pmax_exact):
        st = structure_of_prime(p)
        if p >= 256:
            for d1, e, mf, mg in st:
                if d1 > 1:
                    ordinary -= len(mf if side == "f" else mg) * (p ** (-d1))
        pdf = _prune_pdf(_convolve_pdf(pdf, local_omitted_pdf(p, st, side, p >= 256, up)))
    for j in range(int(math.log2(pmax_exact)), 69):
        lam = 72.0 * 1.26 * (2.0 ** (j + 1)) / ((j + 1) * LN2) * (2.0 ** (-2 * j))
        comp, tot = {}, 0.0
        mass = (j + 1) if up else j
        rnd = math.ceil if up else math.floor
        for t in range(1, 12):
            pr = lam * (2.0 ** (-j * (t - 1))) * (1.0 - 2.0 ** (-j))
            u_m = rnd(mass * t * GRID_RET)
            comp[u_m] = comp.get(u_m, 0.0) + pr
            tot += pr
        for t in range(2, 8):
            pr = 0.5 * lam * lam * (2.0 ** (-j * (t - 2))) * (t - 1)
            u_m = rnd(mass * t * GRID_RET)
            comp[u_m] = comp.get(u_m, 0.0) + pr
            tot += pr
        comp[0] = max(0.0, 1.0 - tot)
        pdf = _prune_pdf(_convolve_pdf(pdf, comp))
        ordinary -= (36.0 if side == "f" else 6.0) * 1.26 * (2.0 ** (j + 1)) / ((j + 1) * LN2) * (2.0 ** (-2 * j))
    return pdf, ordinary


KA_MAX = 16


def _conv_list(a, b):
    out = [0.0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return out


@functools.lru_cache(None)
def window_rows(n, base, iA, mesh=48):
    h = (n - base) / mesh
    wA = [math.log1p(h / (base + i * h)) if i < iA else 0.0 for i in range(mesh)]
    wB = [math.log1p(h / (base + i * h)) if i >= iA else 0.0 for i in range(mesh)]
    PA = [[1.0]]
    for j in range(1, KA_MAX + 1):
        PA.append([x / j for x in _conv_list(PA[-1], wA)])
    PB = [[1.0], wB]
    rows = {}
    for jA in range(KA_MAX + 1):
        for jB in range(2):
            rows[(jA, jB)] = PA[jA] if jB == 0 else _conv_list(PA[jA], PB[1])
    return rows, h


@functools.lru_cache(None)
def ss_cum2(m, base, n, iA, mesh, alpha):
    rows, h = window_rows(n, base, iA, mesh)
    out = [[0.0, 0.0] for _ in range(KA_MAX + 1)]
    tab, ntab = RHO_TAB, len(RHO_TAB)
    for (jA, jB), row in rows.items():
        j = jA + jB
        imax = len(row) - 1
        ifeas = min(imax, int(math.floor((m - j * h - j * base) / h))) if h > 0 else imax
        if ifeas < 0:
            continue
        k0 = (m - j * base + alpha) / (base * RSTEP)
        hk = h / (base * RSTEP)
        s = 0.0
        for i in range(ifeas + 1):
            w = row[i]
            if w:
                k = math.ceil(k0 - i * hk)
                s += w * (1.0 if k <= 0 else (tab[k] if k < ntab else 0.0))
        out[jA][jB] = s
    cum = [[0.0, 0.0] for _ in range(KA_MAX + 1)]
    for jB in range(2):
        acc = 0.0
        for jA in range(KA_MAX + 1):
            acc += out[jA][jB]
            cum[jA][jB] = acc
    return tuple(tuple(x) for x in cum)


def ss_lower2(m, base, n, iA, kA, kB, mesh, alpha):
    m = math.ceil(m)
    if n <= base:
        return rho_lo((m + alpha) / base)
    cum = ss_cum2(m, base, n, iA, mesh, alpha)
    return cum[kA][0] + (cum[kA][1] if kB else 0.0)


@functools.lru_cache(None)
def ss_upper2(m, base, n, iA, kA, kB, mesh, alpha):
    m = math.floor(m)
    if n <= base:
        return rho_hi((m + alpha) / base)
    rows, h = window_rows(n, base, iA, mesh)
    tot = 0.0
    for (jA, jB), row in rows.items():
        if jA > kA or (jB and not kB):
            continue
        j = jA + jB
        for i, w in enumerate(row):
            if not w:
                continue
            left = j * base + i * h
            if left > m:
                continue
            u = (m - left - j * h + alpha) / base
            tot += w * (1.0 if u <= 1.0 else rho_hi(u))
    return min(tot, 1.0)


@functools.lru_cache(None)
def node_retention_loss(side, m, n, iA, kA, kB, alpha, slack=96.0, step=64):
    s_units = int(math.floor((slack - 1.0 - m / 128.0) * GRID_RET))
    pdf_up, _ = omitted_mass_pdf(side, True)
    base = 69.0
    num = sum(pr * (ss_upper2(max(m - math.ceil((x / GRID_RET) / step) * step, 1.0), base, n, iA, kA, kB, 48, alpha) if m - math.ceil((x / GRID_RET) / step) * step > 1.0 else 1.0) for x, pr in pdf_up.items() if x > s_units)
    den = sum(pr * (ss_lower2(max(m - math.floor((x / GRID_RET) / step) * step, 1.0), base, n, iA, kA, kB, 48, alpha) if m - math.floor((x / GRID_RET) / step) * step > 1.0 else 1.0) for x, pr in pdf_up.items() if x <= s_units)
    return num / den if den > 0.0 else 1.0


def node_retention(nf, ng, n, iA, kA, kB, slack=96.0):
    lf = node_retention_loss("f", math.ceil(nf), n, iA, kA, kB, 5.0, slack)
    lg = node_retention_loss("g", math.ceil(ng), n, iA, kA, kB, 1.0, slack)
    _, ordf = omitted_mass_pdf("f", True)
    _, ordg = omitted_mass_pdf("g", True)
    return max(0.0, 1.0 - lf - lg) * max(0.0, ordf + ordg - 1.0) * base_primitive_fraction(), lf, lg


@functools.lru_cache(None)
def miss_side_ideals(u4, log2_rels=64.162, u_main=10.869):
    u = u4 / 4.0
    tot = 0.0
    for b10 in range(400, 690):
        b = b10 / 10.0
        lam = (2.0 ** (log2_rels - b)) * 2.0 * kappa_rho(u_main, b / 69.0)
        tot += kappa_rho(u, b / 69.0) * math.log((b + 0.1) / b) * math.exp(-lam) * (1.0 + lam + 0.5 * lam * lam)
    return tot


def missing_ideals(nf, ng):
    return (miss_side_ideals(math.ceil(4.0 * nf / 69.0), 64.162, 750.0 / 69.0) if nf else 0.0) + (
        miss_side_ideals(math.ceil(4.0 * ng / 69.0), 64.162, 498.5 / 69.0) if ng else 0.0
    )


def surv_upper_prob(m, base, T):
    if T <= 0.0:
        return 0.0
    if T >= m:
        return 1.0
    return min(1.0, max(1.0, T / base) * rho_hi((math.floor(m) - T) / base))


QBINS = [
    (1 / 6, -26.4 + 23, 14.5 + 6), (1 / 6, -11.2 + 23, 14.5 + 6), (1 / 6, -0.1 + 23, 14.5 + 6),
    (1 / 6, 10.5 + 23, 14.5 + 6), (1 / 6, 23.2 + 23, 14.5 + 6), (1 / 6, 86.6 + 23, 14.5 + 6),
]
QBINS_LO = [
    (-140.05 - 23, -54.56 - 6), (-26.4 - 23, -54.56 - 6), (-11.2 - 23, -54.56 - 6),
    (-0.1 - 23, -54.56 - 6), (10.5 - 23, -54.56 - 6), (23.2 - 23, -54.56 - 6),
]
BOOT_BINS = [(0.5, 1.4 + 4), (0.4, 6.8 + 4), (0.1, 12.7 + 4)]

_PLAN_TABLE = {
    69: ((69, 48, 1, 0, 288, 1), (69, 48, 1, 0, 384, 1)),
    75: ((75, 42, 1, 0, 224, 1), (75, 48, 1, 0, 256, 1)),
    80: ((75, 42, 1, 0, 320, 1), (75, 42, 1, 0, 384, 1)),
    85: ((80, 24, 2, 0, 192, 1), (75, 48, 2, 0, 224, 1)),
    90: ((75, 48, 2, 0, 256, 1), (75, 48, 2, 0, 320, 0)),
    95: ((80, 24, 2, 0, 384, 1), (80, 18, 1, 1, 384, 1)),
    100: ((85, 12, 1, 1, 384, 1), (80, 24, 2, 1, 384, 0)),
    105: ((85, 18, 2, 1, 320, 1), (85, 18, 2, 1, 384, 0)),
    110: ((85, 18, 2, 1, 384, 0), (90, 24, 2, 1, 384, 1)),
    115: ((90, 24, 2, 1, 384, 1), (95, 24, 2, 1, 384, 0)),
    120: ((95, 24, 2, 1, 384, 1), (95, 48, 3, 0, 384, 1)),
    125: ((100, 24, 2, 1, 384, 0), (100, 24, 3, 1, 384, 0)),
    130: ((100, 24, 3, 1, 384, 1), (100, 24, 4, 1, 384, 0)),
    135: ((100, 30, 3, 1, 384, 1), (110, 30, 3, 1, 384, 1)),
    140: ((100, 36, 3, 1, 384, 0), (110, 36, 3, 1, 384, 0)),
    145: ((110, 30, 3, 1, 384, 1), (110, 36, 4, 1, 384, 1)),
    150: ((110, 36, 3, 1, 384, 0), (110, 42, 4, 1, 384, 0)),
    155: ((110, 30, 4, 1, 384, 0), (125, 30, 4, 1, 384, 1)),
    160: ((115, 30, 4, 1, 384, 0), (130, 30, 4, 1, 384, 0)),
    165: ((125, 24, 4, 1, 384, 0), (130, 36, 4, 1, 384, 0)),
    170: ((125, 30, 4, 1, 384, 1), (135, 36, 4, 1, 384, 0)),
    175: ((130, 30, 4, 1, 384, 1), (150, 30, 4, 1, 384, 0)),
    180: ((130, 36, 4, 1, 320, 0), (130, 36, 6, 1, 320, 0)),
    185: ((140, 30, 4, 1, 384, 0), (140, 30, 6, 1, 384, 0)),
    190: ((145, 30, 4, 1, 384, 0), (160, 24, 6, 1, 384, 1)),
    195: ((150, 30, 4, 1, 384, 0), (160, 42, 4, 1, 384, 0)),
    200: ((150, 36, 4, 1, 320, 0), (150, 36, 6, 1, 288, 0)),
    205: ((165, 30, 4, 1, 384, 1), (165, 30, 6, 1, 384, 1)),
    210: ((160, 36, 4, 1, 320, 0), (170, 30, 6, 1, 384, 0)),
}


def evaluate_descent_plan():
    lmax = 215
    bins = [(69, 75)] + [(l, l + 5) for l in range(75, lmax, 5)]
    W, PHI, DET = {}, {}, {}
    pairs = 2.0 * (2.0 ** log2_li(69.0))
    cfmax, cgmax = max(b[1] for b in QBINS), max(b[2] for b in QBINS)
    hits_per_cell = 2.0 * math.log(69.0 / 8.0)

    def wmax(side, n):
        return max([W[(b, side)] for b in bins if b[0] < n - 1e-9] + [0.0])

    def phimax(side, n):
        return max([PHI[(b, side)] for b in bins if b[0] < n - 1e-9] + [0.0])

    def miss_terms():
        low = [k for k in W if k[0][0] <= 75]
        return max([PHI[k] for k in low] + [0.0]), max([W[k] for k in low] + [0.0])

    for b in bins:
        for s_idx, side in enumerate(("f", "g")):
            n, iA, kA, kB, R, thin = _PLAN_TABLE[b[0]][s_idx]
            n = float(n)
            S = math.log2(CELLS[R])
            c, c_lo = (b[1] + S) / 12.0, (b[0] + S) / 12.0
            h = (n - 69.0) / 48.0 if n > 69.0 else 0.0
            nA = 69.0 + iA * h if n > 69.0 else 69.0
            T = kA * nA + kB * n + 96.0
            pr = sv = 0.0
            for (w, cf, cg), (cf_lo, cg_lo) in zip(QBINS, QBINS_LO):
                nf = 72.0 * c + cf - (b[0] if side == "f" else 0.0)
                ng = 12.0 * c + LOG2U6 + cg - (b[0] if side == "g" else 0.0)
                pr += w * ss_lower2(nf, 69.0, n, iA, kA, kB, 48, 5.0) * ss_lower2(ng, 69.0, n, iA, kA, kB, 48, 1.0)
                nf_lo = 72.0 * c_lo + cf_lo - (b[1] if side == "f" else 0.0)
                ng_lo = 12.0 * c_lo + LOG2U6 + cg_lo - (b[1] if side == "g" else 0.0)
                sv += w * surv_upper_prob(nf_lo, 69.0, T) * surv_upper_prob(ng_lo, 69.0, T)

            nf_max = 72.0 * c + cfmax - (b[0] if side == "f" else 0.0)
            ng_max = 12.0 * c + LOG2U6 + cgmax - (b[0] if side == "g" else 0.0)
            det = 1.0 - 0.001 * ((nf_max + ng_max) / 20.0 + 2 * (kA + kB))
            miss = missing_ideals(nf_max, ng_max)
            thin = bool(thin or b[0] <= 75)
            ret, lf, lg = node_retention(nf_max, ng_max, n, iA, kA, kB, 96.0)
            mu = CELLS[R] * pr * 0.25 * det * (math.exp(-miss) if thin else 1.0) * ret
            assert mu >= 16.0

            surv_cost = (2.0**26) + (2.0**24.5) + ecm_instr(69.0, nf_max, 0.999) + ecm_instr(69.0, ng_max, 0.999) + 2.0 * (ecm_instr(nA, T, 0.999) + 2.0**24) + 2.0**34
            survivors = CELLS[R] * (sv + pr)
            sieve = CELLS[R] * (hits_per_cell * 1100.0 + 10.0 + 2.0 + 2400.0) + pairs * (2.0**25.1)
            own = (2.0**36) + sieve + survivors * surv_cost

            pm, wm = miss_terms()
            miss_ch = 0.0 if thin else miss
            psi = kA * (phimax("f", nA) + phimax("g", nA)) + kB * (phimax("f", n) + phimax("g", n)) + miss_ch * pm
            assert psi <= 0.80
            retry = 1.0 / (1.0 - psi)
            sub = kA * (wmax("f", nA) + wmax("g", nA)) + kB * (wmax("f", n) + wmax("g", n)) + miss_ch * wm
            W[(b, side)] = own + retry * sub
            win = (kA if nA > 69.0 + 1e-9 else 0) + (kB if n > nA + 1e-9 else 0)
            r_nb = 0.5 if win <= 1 else 1.0
            PHI[(b, side)] = max((1.0 + mu * (1.0 - psi) / r_nb) ** (-r_nb), 0.02)
            assert PHI[(b, side)] <= 0.05
            DET[(b, side)] = dict(mu=round(mu, 2), own_log2=round(math.log2(own), 3), W_log2=round(math.log2(W[(b, side)]), 3), psi=round(psi, 4), phi=round(PHI[(b, side)], 5), ret=round(ret, 4))

    n1, iA, kA, Rb = (215, 42, 12, 4)
    Sb = math.log2(CELLS[Rb])
    cb = (6.0 * 380.99700 + Sb) / 12.0
    nA_b = 69.0 + iA * ((n1 - 69.0) / 48.0)
    HALF_LOG2U6 = 0.5 * LOG2U6
    m_max = 12.0 * cb + HALF_LOG2U6 + max(x[1] for x in BOOT_BINS)
    pr_b = sum(w * ss_lower2(12.0 * cb + HALF_LOG2U6 + cg, 69.0, float(n1), iA, kA, 1, 48, 1.0) for w, cg in BOOT_BINS)
    det_b = 1.0 - 0.001 * (m_max / 20.0 + kA + 1)
    _, ordf = omitted_mass_pdf("f", True)
    _, ordg = omitted_mass_pdf("g", True)
    cnt_frac = max(0.0, ordf + ordg - 1.0) * base_primitive_fraction()
    trials = 16.0 / (pr_b * det_b * 0.25 * cnt_frac)
    balls = 1.0 + trials / CELLS[Rb]
    trial_cost = (2.0**26) + ecm_instr(nA_b, m_max, 0.999) + (2.0**24) + (2.0**34)
    pm, wm = miss_terms()
    mb = missing_ideals(0, m_max)
    psi_b = kA * phimax("g", nA_b) + phimax("g", float(n1)) + mb * pm
    assert psi_b <= 0.80
    retry_b = 1.0 / (1.0 - psi_b)
    sub_b = kA * wmax("g", nA_b) + wmax("g", float(n1)) + mb * wm
    one_log = (2.0**36) * balls + trials * trial_cost + retry_b * sub_b
    phi_boot = psi_b**16
    two_log_units = math.log2(2.0 * one_log / UNIT)

    return {
        "bins_count": len(bins) * 2,
        "min_mu": min(d["mu"] for d in DET.values()),
        "max_psi": max(max(d["psi"] for d in DET.values()), round(psi_b, 4)),
        "max_phi": max(d["phi"] for d in DET.values()),
        "min_retention": min(d["ret"] for d in DET.values()),
        "balanced_boot": {
            "n1": n1,
            "nA": round(nA_b, 2),
            "kA": kA,
            "R": Rb,
            "m_max_bits": round(m_max, 2),
            "pr_log2": round(math.log2(pr_b), 3),
            "trials_log2": round(math.log2(trials), 3),
            "trial_cost_log2": round(math.log2(trial_cost), 3),
            "own_work_log2_instr": round(math.log2(trials * trial_cost), 3),
            "subtree_work_log2_instr": round(math.log2(sub_b), 3),
            "psi_boot": round(psi_b, 4),
            "phi_boot_16": phi_boot,
            "one_log_log2_instr": round(math.log2(one_log), 4),
        },
        "two_log_expected_units_log2": round(two_log_units, 4),
    }


def run_toy_fixed_ball(seed=1000):
    rng = random.Random(seed)
    small = primes_range(2, 1 << 10)
    while True:
        q = rng.randrange(1 << 13, 1 << 14) | 1
        if not mr_probable_prime(q):
            continue
        opts = [(s, t) for s in roots_mod_q(H_POLY, q, rng) for t in roots_mod_q([4 * s % q, 0, 1], q, rng)]
        if opts:
            s, t = rng.choice(opts)
            break
    w = [pow(s, j, q) for j in range(6)] + [(-t * pow(s, j, q)) % q for j in range(6)]
    B0 = [[q] + [0] * 11] + [[(-w[j]) % q] + [int(k == j) for k in range(1, 12)] for j in range(1, 12)]
    B, _ = lll_12d(B0, prec=120)

    cells, c = [], [0] * 12

    def rec(d, rem, pos):
        if d < 0:
            if pos:
                cells.append(tuple(c))
            return
        r = math.isqrt(rem)
        for x in range(0 if not pos else -r, r + 1):
            if not pos and x == 0:
                c[d] = 0
                rec(d - 1, rem, False)
            elif pos or x > 0:
                c[d] = x
                rec(d - 1, rem - x * x, True)

    rec(11, 4, False)
    raw, distinct, model_sum = 0, set(), 0.0
    f_toy_terms = [(4, [1]), (2, [0, 1]), (0, [0, 0, 1])]
    win_primes = primes_range((1 << 10) + 1, 1 << 13)

    def fac_win(n_val):
        fs = []
        for p in small:
            if p * p > n_val:
                break
            while n_val % p == 0:
                fs.append(p)
                n_val //= p
        if n_val == 1:
            return True, tuple(fs)
        if n_val <= (1 << 10):
            fs.append(n_val)
            return True, tuple(fs)
        if n_val >> 26:
            return False, ()
        if mr_probable_prime(n_val):
            return (True, tuple(fs + [n_val])) if n_val <= (1 << 13) else (False, ())
        w_cnt = 0
        for p in win_primes:
            if p * p > n_val:
                break
            while n_val % p == 0:
                fs.append(p)
                n_val //= p
                w_cnt += 1
                if w_cnt > 2:
                    return False, ()
        if n_val > 1:
            if n_val > (1 << 13) or not mr_probable_prime(n_val):
                return False, ()
            fs.append(n_val)
            w_cnt += 1
        return (w_cnt <= 2), tuple(sorted(fs))

    for cell in cells:
        v = [sum(cell[i] * B[i][j] for i in range(12)) for j in range(12)]
        a, b = v[:6], v[6:]
        ap, bp = [[1]], [[1]]
        for _ in range(4):
            ap.append(zy_mul(ap[-1], a))
            bp.append(zy_mul(bp[-1], b))
        fs_el = [0] * 6
        for d, cf in f_toy_terms:
            term = zy_mul(cf, zy_mul(ap[d], bp[4 - d]))
            fs_el = [x + y for x, y in zip(fs_el, term)]
        gs_el = [x + y for x, y in zip(ap[2], zy_mul([0, 4], bp[2]))]
        nf, ng = zy_norm(fs_el), zy_norm(gs_el)
        if nf == 0 or ng == 0:
            continue
        assert ng % q == 0
        ng //= q
        model_sum += ss_lower2(math.log2(nf), 10.0, 13.0, 48, 2, 0, 48, 0.0) * ss_lower2(math.log2(ng), 10.0, 13.0, 48, 2, 0, 48, 0.0)
        okf, ff = fac_win(nf)
        if not okf:
            continue
        okg, gf = fac_win(ng)
        if not okg:
            continue
        raw += 1
        distinct.add((ff, gf))

    assert len(distinct) > 0
    return {
        "q": q,
        "cells": len(cells),
        "raw_relations": raw,
        "distinct_relations": len(distinct),
        "model_mean_times_udup": round(model_sum * 0.25, 2),
    }


def parse_input():
    params = {"n_draws": 2500, "seed": 1}
    if len(sys.argv) > 1 and sys.argv[1].startswith("{"):
        params.update(json.loads(sys.argv[1]))
        return params
    if not sys.stdin.isatty():
        raw = sys.stdin.read().strip()
        if raw:
            data = json.loads(raw)
            if isinstance(data, list) and data:
                data = data[0]
            elif isinstance(data, dict) and "inputs" in data and data["inputs"]:
                data = data["inputs"][0]
            if isinstance(data, dict):
                params.update(data)
    return params


def main():
    t0 = time.time()
    cfg_in = parse_input()
    n_draws = int(cfg_in.get("n_draws", 2500))
    seed = int(cfg_in.get("seed", 1))

    tower_info = verify_tower_isomorphism_and_fields()
    alpha_info = compute_prime_structure_and_alphas(512)
    ellipsoid_cnts = certify_ellipsoid_counts()
    main_phase = evaluate_main_phase_live(n_draws=n_draws, n_bound_checks=24)
    filt_info = run_filtering_and_missing_audit(seed=seed, log2U=12)
    ecm_audit = validate_ecm_stage1_stage2(trials=50, seed=seed)
    lattice_audit = audit_balanced_boot_and_both_side_specialq(seed=31001 + seed)
    descent_info = evaluate_descent_plan()
    toy_info = run_toy_fixed_ball(seed=1000 + seed)

    log2_relcol = main_phase["log2_relcol_units"]
    log2_la = filt_info["log2_la_units"]
    two_log_exp = descent_info["two_log_expected_units_log2"]
    descent_cap = 120.5

    markov_abort = 2.0 ** (two_log_exp - descent_cap)
    phi_boot = descent_info["balanced_boot"]["phi_boot_16"]
    p_main = 1.0 - (2.0**-190)
    success_prob = p_main * (1.0 - 2.0 * phi_boot - markov_abort)
    assert success_prob >= 0.39

    total_units = (
        (2.0**log2_relcol)
        + (2.0**log2_la)
        + (2.0**descent_cap)
        + (2.0**90.0)
        + (2.0**88.0)
        + (2.0**72.0)
        + (2.0**80.0)
        + (2.0**11.0)
    )
    total_log2_units = math.log2(total_units)
    claimed_score = 123.0
    assert total_log2_units <= claimed_score

    out = {
        "status": "ok",
        "elapsed_seconds": round(time.time() - t0, 2),
        "tower_and_field_certificates": tower_info,
        "murphy_alpha_and_prime_structure": alpha_info,
        "ellipsoid_lattice_count_certificates": ellipsoid_cnts,
        "main_phase_relation_collection": main_phase,
        "filtering_and_linear_algebra": filt_info,
        "ecm_stage1_stage2_validation": ecm_audit,
        "balanced_boot_and_specialq_audit": lattice_audit,
        "descent_evaluation": descent_info,
        "toy_fixed_ball_validation": toy_info,
        "ledger_summary": {
            "log2_relcol_units": log2_relcol,
            "log2_la_units": log2_la,
            "descent_expected_two_logs_units": two_log_exp,
            "descent_cap_units": descent_cap,
            "markov_abort_prob": round(markov_abort, 4),
            "modelled_success_probability": round(success_prob, 4),
            "total_log2W_full_ledger": round(total_log2_units, 6),
            "submitted_claim_bound": claimed_score,
        },
    }
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()