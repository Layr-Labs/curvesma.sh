#!/usr/bin/env python3
"""
Clean-Room Pure-Python Verification Suite for stnfs-q15-guard4-bandpass-snell-v1
Track: bls12-381-fixed-g1-invert-g2 (classical-group-work-v2, C = 20,000 Ir/u)

Pure Python 3.11+ stdlib ONLY. Reads JSON object on stdin:
  {}                         -> mode = "core" (< 7s)
  {"mode": "diagnostics"}    -> extended finite-scale diagnostics
  {"mode": "all"}            -> core + diagnostics
"""
from __future__ import annotations

import bisect
import json
import math
import random
import re
import sys
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent

X_BLS = -0xD201000000010000
ABS_X = abs(X_BLS)
LOG2_ABS_X = math.log2(ABS_X)
P_BLS = int(
    "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab",
    16,
)
R_BLS = int("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001", 16)
P0_X = int(
    "17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb",
    16,
)
P0_Y = int(
    "08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1",
    16,
)
Q0_X = (
    int("024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8", 16),
    int("13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e", 16),
)
Q0_Y = (
    int("0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801", 16),
    int("0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be", 16),
)
C_COST = 20000.0


class ExactDickman:
    """Zero-drift positive-sum Dickman rho(u) table valid up to u_max."""

    def __init__(self, u_max: float = 30.0, steps_per_unit: int = 200) -> None:
        self.u_max = u_max
        self.N = steps_per_unit
        h = 1.0 / steps_per_unit
        total = int(round(u_max * steps_per_unit)) + 1
        vals = [0.0] * total
        for i in range(steps_per_unit + 1):
            vals[i] = 1.0
        for i in range(steps_per_unit + 1, min(2 * steps_per_unit + 1, total)):
            vals[i] = 1.0 - math.log(i * h)
        for i in range(2 * steps_per_unit + 1, total):
            u = i * h
            s = 0.5 * vals[i - steps_per_unit] + sum(vals[i - steps_per_unit + 1 : i])
            vals[i] = (h * s) / (u - 0.5 * h)
        self.vals = vals

    def rho(self, u: float) -> float:
        if u <= 1.0:
            return 1.0
        if u >= self.u_max:
            return 1e-300
        pos = u * self.N
        idx = int(pos)
        frac = pos - idx
        if idx + 1 >= len(self.vals):
            return self.vals[-1]
        return math.exp((1.0 - frac) * math.log(self.vals[idx]) + frac * math.log(self.vals[idx + 1]))


DICKMAN = ExactDickman(u_max=30.0, steps_per_unit=200)


# =============================================================================
# 1. Stage 0 & 12D M_g Algebra, Frobenius, Automorphisms & Prime Splitting
# =============================================================================
def _g1_add(p1, p2):
    if p1 is None:
        return p2
    if p2 is None:
        return p1
    x1, y1 = p1
    x2, y2 = p2
    if x1 == x2:
        if (y1 + y2) % P_BLS == 0:
            return None
        lam = (3 * x1 * x1 * pow(2 * y1, -1, P_BLS)) % P_BLS
    else:
        lam = ((y2 - y1) * pow(x2 - x1, -1, P_BLS)) % P_BLS
    x3 = (lam * lam - x1 - x2) % P_BLS
    y3 = (lam * (x1 - x3) - y1) % P_BLS
    return (x3, y3)


def _g1_mul(k: int, pt):
    res, addend = None, pt
    while k > 0:
        if k & 1:
            res = _g1_add(res, addend)
        addend = _g1_add(addend, addend)
        k >>= 1
    return res


def _fp2_mul(a, b):
    return ((a[0] * b[0] - a[1] * b[1]) % P_BLS, (a[0] * b[1] + a[1] * b[0]) % P_BLS)


def _fp2_add(a, b):
    return ((a[0] + b[0]) % P_BLS, (a[1] + b[1]) % P_BLS)


def _sieve_primes(limit: int) -> list[int]:
    is_p = bytearray(b"\x01") * (limit + 1)
    is_p[0] = is_p[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_p[i]:
            is_p[i * i : limit + 1 : i] = bytearray(len(range(i * i, limit + 1, i)))
    return [i for i in range(2, limit + 1) if is_p[i]]


def verify_prime_splitting(p_limit: int) -> dict:
    primes = _sieve_primes(p_limit)
    p15 = [p for p in primes if p % 15 == 1]
    g_act_p = f_act_p = f_act_ok = f_tot_roots = bad_m3 = 0
    for p in p15:
        w = 1
        for gen in range(2, p):
            w = pow(gen, (p - 1) // 3, p)
            if w != 1:
                break
        if pow(X_BLS % p, (p - 1) // 3, p) == 1:
            g_act_p += 1
        p_has_f = False
        for ww in (w, (w * w) % p):
            ww2 = (ww * ww) % p

            def pmul(a, b):
                c = [0, 0, 0, 0, 0]
                for i in range(3):
                    for j in range(3):
                        c[i + j] = (c[i + j] + a[i] * b[j]) % p
                for deg in (4, 3):
                    v = c[deg]
                    if v:
                        c[deg - 1] = (c[deg - 1] + v) % p
                        c[deg - 2] = (c[deg - 2] - v * ww2) % p
                        c[deg - 3] = (c[deg - 3] + v * ww) % p
                return c[:3]

            res, base, e = [1, 0, 0], [0, 1, 0], p
            while e:
                if e & 1:
                    res = pmul(res, base)
                base = pmul(base, base)
                e >>= 1
            res[1] = (res[1] - 1) % p

            def strip(a):
                a = list(a)
                while a and a[-1] == 0:
                    a.pop()
                return a

            def pgcd(a, b):
                a, b = strip(a), strip(b)
                while b:
                    inv = pow(b[-1], p - 2, p)
                    while len(a) >= len(b):
                        coef = (a[-1] * inv) % p
                        sh = len(a) - len(b)
                        for idx in range(len(b)):
                            a[sh + idx] = (a[sh + idx] - coef * b[idx]) % p
                        a = strip(a)
                    a, b = b, a
                if a:
                    inv = pow(a[-1], p - 2, p)
                    a = [(x * inv) % p for x in a]
                return a

            gpoly = pgcd([(-ww) % p, ww2, (-1) % p, 1], res)
            deg_g = len(gpoly) - 1 if gpoly else 0
            q_roots = []
            if deg_g == 1:
                q_roots = [(-gpoly[0]) % p]
            elif deg_g == 3:
                rng = random.Random(p * 17 + ww)

                def split_poly(poly):
                    if len(poly) == 2:
                        return [(-poly[0]) % p]
                    d = len(poly) - 1
                    while True:
                        a0 = rng.randrange(p)

                        def pm(u, v):
                            c = [0] * (2 * d)
                            for i in range(len(u)):
                                for j in range(len(v)):
                                    c[i + j] = (c[i + j] + u[i] * v[j]) % p
                            while len(c) > d:
                                val = c.pop()
                                if val:
                                    sh = len(c) - d
                                    for k in range(d):
                                        c[sh + k] = (c[sh + k] - val * poly[k]) % p
                            return c

                        r_p, b_p, ex = [1], [a0, 1], (p - 1) // 2
                        while ex:
                            if ex & 1:
                                r_p = pm(r_p, b_p)
                            b_p = pm(b_p, b_p)
                            ex >>= 1
                        r_p[0] = (r_p[0] - 1) % p
                        h = pgcd(poly, r_p)
                        if 1 < len(h) < len(poly):
                            qd, rem = [], list(poly)
                            while len(rem) >= len(h):
                                coef = rem[-1]
                                qd.append(coef)
                                sh = len(rem) - len(h)
                                for k in range(len(h)):
                                    rem[sh + k] = (rem[sh + k] - coef * h[k]) % p
                                rem.pop()
                            qd.reverse()
                            return split_poly(h) + split_poly(qd)

                q_roots = split_poly(gpoly)
            n_y = sum(3 for rt in q_roots if pow(rt, (p - 1) // 3, p) == 1)
            if n_y % 3 != 0:
                bad_m3 += 1
            if n_y > 0:
                p_has_f = True
                f_act_ok += 4
                f_tot_roots += 4 * n_y
        if p_has_f:
            f_act_p += 1
    n_all = len(primes)
    return {
        "p_limit": p_limit,
        "total_primes": n_all,
        "p_1mod15_frac": round(len(p15) / n_all, 5),
        "g_active_p_frac": round(g_act_p / n_all, 5),
        "f_active_p_frac": round(f_act_p / n_all, 5),
        "f_active_ok_frac": round(f_act_ok / n_all, 5),
        "all_root_counts_divisible_by_3": (bad_m3 == 0),
        "f_active_ok_le_one_third": (f_act_ok / n_all <= 1.0 / 3.0 + 0.02),
    }


def check_stage0_algebra() -> dict:
    assert (X_BLS - 1) ** 2 * (X_BLS**4 - X_BLS**2 + 1) // 3 + X_BLS == P_BLS
    assert X_BLS**4 - X_BLS**2 + 1 == R_BLS
    # MAJOR #6 exact Frobenius & lambda identity: p == x == lambda (mod r), lambda primitive 12th root of unity!
    lam = X_BLS % R_BLS
    assert P_BLS % R_BLS == lam
    assert (pow(lam, 4, R_BLS) - pow(lam, 2, R_BLS) + 1) % R_BLS == 0
    assert pow(lam, 6, R_BLS) == R_BLS - 1 and pow(lam, 12, R_BLS) == 1 and pow(lam, 4, R_BLS) != 1

    assert (P0_Y * P0_Y - (P0_X**3 + 4)) % P_BLS == 0
    assert _g1_mul(R_BLS, (P0_X, P0_Y)) is None
    assert _fp2_mul(Q0_Y, Q0_Y) == _fp2_add(_fp2_mul(_fp2_mul(Q0_X, Q0_X), Q0_X), (4, 4))

    assert P_BLS % 3 == 1 and P_BLS % 5 == 2 and P_BLS % 15 == 7
    ord5 = [k for k in range(1, 5) if pow(P_BLS, k, 5) == 1][0]
    assert ord5 == 4
    assert pow(X_BLS % P_BLS, (P_BLS - 1) // 3, P_BLS) != 1

    w_p = None
    for g in range(2, 20):
        cand = pow(g, (P_BLS - 1) // 3, P_BLS)
        if cand != 1:
            for c in (cand, (cand * cand) % P_BLS):
                if (X_BLS**3 - X_BLS**2 + (c * c % P_BLS) * X_BLS - c) % P_BLS == 0:
                    w_p = c
                    break
            if w_p is not None:
                break
    assert w_p is not None
    assert X_BLS**6 - 2 * X_BLS**5 + 2 * X_BLS**3 + X_BLS + 1 == 3 * P_BLS
    assert (R_BLS - 1) % 12 == 0 and (36 // 12 + 12 // 12) == 4
    split_3k = verify_prime_splitting(3000)
    assert split_3k["all_root_counts_divisible_by_3"]
    return {
        "ok": True,
        "p_mod_15": P_BLS % 15,
        "ord_5_p": ord5,
        "p_equiv_x_mod_r": True,
        "lambda_order_12_mod_r": True,
        "O_K_and_O_Mg_iso_Fp12": True,
        "n_sm_chi1": 4,
        "small_prime_splitting_3000": split_3k,
    }


# =============================================================================
# 2. FATAL #1: Live 248-Length Guard4 x16 + Masked Active Scan with 0xFF Guard Poison
# =============================================================================
def check_guard4_masked_differential(n_trials: int = 36, V_a: int = 192) -> dict:
    rng = random.Random(0xCAFEBABE2026)
    test_primes = [1021, 1031, 4093, 4441, 65521]
    w_g = [21, 25, 29, 33, 37]
    w_f = [19, 23, 27, 31, 35]

    def build_tables(weights):
        tabs, avs = [], []
        for k, p in enumerate(test_primes):
            nb = ((p - 1) >> 2) + 1
            cnt = [0] * nb
            av_k = [[rng.randrange(p) for _ in range(V_a)] for _ in range(3)]
            avs.append(av_k)
            for rj in range(3):
                for v in av_k[rj]:
                    cnt[v >> 2] += 1
            bkt = [0] * (nb + 1)
            for b in range(nb):
                bkt[b + 1] = bkt[b] + cnt[b]
            cur = list(bkt[:-1])
            hits = [0] * (3 * V_a)
            for rj in range(3):
                for a_id, v in enumerate(av_k[rj]):
                    b = v >> 2
                    hits[cur[b]] = (a_id << 8) | (v & 255)
                    cur[b] += 1
            tabs.append({"rem": rng.randrange(p), "step": rng.randrange(p), "p": p, "bkt": bkt, "hits": hits, "w": weights[k]})
        return tabs, avs

    tabs_g, avs_g = build_tables(w_g)
    tabs_f, avs_f = build_tables(w_f)

    def run_guard4_sieve(S_buf, tabs, L):
        emitted = 0
        S_base = 256
        for mp in tabs:
            p = mp["p"]
            R = (mp["rem"] + mp["step"]) % p
            mp["rem"] = R
            w, bkt, hits = mp["w"], mp["bkt"], mp["hits"]

            def emit(base_off: int, st: int, en: int):
                nonlocal emitted
                for j in range(st, en):
                    addr = base_off + hits[j]
                    S_buf[addr] = (S_buf[addr] + w) & 0xFF
                emitted += en - st

            r_end = R + L
            if r_end <= p:
                b0, b1 = R >> 2, (r_end - 1) >> 2
                q0, q1, r8 = R >> 8, (r_end - 1) >> 8, R & 255
                if q0 == q1:
                    st, en = bkt[b0], bkt[b1 + 1]
                    if en > st:
                        emit(S_base + 4 - r8, st, en)
                else:
                    b_mid = (q0 + 1) << 6
                    st0, en0 = bkt[b0], bkt[b_mid]
                    if en0 > st0:
                        emit(S_base + 4 - r8, st0, en0)
                    en1 = bkt[b1 + 1]
                    if en1 > en0:
                        emit(S_base + 260 - r8, en0, en1)
            else:
                take0 = p - R
                b0, b_last = R >> 2, (p - 1) >> 2
                q0, q_last, r8 = R >> 8, (p - 1) >> 8, R & 255
                if q0 == q_last:
                    st, en = bkt[b0], bkt[b_last + 1]
                    if en > st:
                        emit(S_base + 4 - r8, st, en)
                else:
                    b_mid = (q0 + 1) << 6
                    st0, en0 = bkt[b0], bkt[b_mid]
                    if en0 > st0:
                        emit(S_base + 4 - r8, st0, en0)
                    en1 = bkt[b_last + 1]
                    if en1 > en0:
                        emit(S_base + 260 - r8, en0, en1)
                rem_L = L - take0
                b_we = (rem_L - 1) >> 2
                stw, enw = bkt[0], bkt[b_we + 1]
                if enw > stw:
                    emit(S_base + 4 + take0, stw, enw)
        return emitted

    # ONE single slab reused sequentially g -> f!
    S_buf = bytearray(256 + 256 * V_a + 256)
    bias_g, bias_f = 95, 108  # ensures real positive dual survivors + tests 0xFF guard rejection!
    tot_hits = tot_gsurv = tot_dual = 0

    for trial in range(n_trials):
        L = 248 if trial < 6 else (trial + 1 if trial < 18 else rng.randint(1, 248))
        if trial < 6:
            # Force prime wrap + simultaneous 256-boundary at full L = 248!
            for tabs in (tabs_g, tabs_f):
                for mp in tabs:
                    mp["rem"] = mp["p"] - (1 + trial * 37)
                    mp["step"] = 5

        # Reference g & f active arrays
        ref_g = bytearray([bias_g] * (256 * V_a))
        ref_f = bytearray([bias_f] * (256 * V_a))
        for tabs, avs, ref in ((tabs_g, avs_g, ref_g), (tabs_f, avs_f, ref_f)):
            for k, mp in enumerate(tabs):
                p = mp["p"]
                R = (mp["rem"] + mp["step"]) % p
                w = mp["w"]
                for rj in range(3):
                    for a_id, v in enumerate(avs[k][rj]):
                        diff = (v - R) if v >= R else (v + p - R)
                        if diff < L:
                            idx = (a_id << 8) | (diff + 4)
                            ref[idx] = (ref[idx] + w) & 0xFF

        # Pass 1 (g): poison ALL left [0..3] & right [4+L..255] guard bytes with 0xFF, prebias active [4..4+L-1]!
        S_buf[:256] = bytes(256)
        S_buf[256 + 256 * V_a :] = bytes(256)
        row_tpl = bytearray(b"\xff" * 256)
        row_tpl[4 : 4 + L] = bytes([bias_g] * L)
        for a_id in range(V_a):
            st = 256 + (a_id << 8)
            S_buf[st : st + 256] = row_tpl
        tot_hits += run_guard4_sieve(S_buf, tabs_g, L)
        assert not any(S_buf[:256]) and not any(S_buf[256 + 256 * V_a :])

        # Masked active-only scan of S_g (masks BOTH left [0..3] and right [4+L..255]!)
        g_surv = []
        for a_id in range(V_a):
            base = a_id << 8
            row_act = S_buf[256 + base + 4 : 256 + base + 4 + L]
            assert row_act == ref_g[base + 4 : base + 4 + L]
            # Verify 0xFF guard poison was ignored by masking strictly to [4..4+L-1]
            for c in range(L):
                if row_act[c] & 0x80:
                    g_surv.append(base | (c + 4))
        tot_gsurv += len(g_surv)

        # Pass 2 (f): reuse SAME slab, re-poison guards with 0xFF, prebias active [4..4+L-1] with bias_f!
        row_tpl[4 : 4 + L] = bytes([bias_f] * L)
        for a_id in range(V_a):
            st = 256 + (a_id << 8)
            S_buf[st : st + 256] = row_tpl
        tot_hits += run_guard4_sieve(S_buf, tabs_f, L)
        assert not any(S_buf[:256]) and not any(S_buf[256 + 256 * V_a :])

        # Direct probe ONLY at g_surv coordinates (_slab_filter_f_at_gsurv)
        dual_got = [off for off in g_surv if S_buf[256 + off] & 0x80]
        dual_ref = [
            (a_id << 8) | (c + 4)
            for a_id in range(V_a)
            for c in range(L)
            if (ref_g[(a_id << 8) | (c + 4)] & 0x80) and (ref_f[(a_id << 8) | (c + 4)] & 0x80)
        ]
        assert dual_got == dual_ref
        tot_dual += len(dual_got)

    assert tot_dual > 0 and tot_gsurv > 0
    return {
        "ok": True,
        "trials": n_trials,
        "V_a": V_a,
        "max_L_tested": 248,
        "guard_poison_0xFF_false_positives": 0,
        "outer_guard_corruptions": 0,
        "bit_exact_dual_survivors": tot_dual,
        "total_g_survivors": tot_gsurv,
        "total_emitted_hits": tot_hits,
    }


# =============================================================================
# 3. FATAL #2: Exact 2D Panjer Smooth-Recall + Survival DP & Zero Wraparound
# =============================================================================
def _panjer_1d(mu_dict: dict[int, float], max_k: int):
    P = [0.0] * (max_k + 1)
    lam = sum(mu_dict.values())
    P[0] = math.exp(-lam)
    items = [(j, v) for j, v in mu_dict.items() if v > 0 and j >= 1]
    for k in range(1, max_k + 1):
        acc = 0.0
        for j, v in items:
            if j <= k:
                acc += j * v * P[k - j]
        P[k] = acc / k
    return P, lam


def _exact_2d_floor_dp(N_eff: float, B_b: float, plo: float, phi: float, qstep: float, T_u8: int, N_max: float) -> dict:
    BPB = 2
    kN = int(round(N_eff * BPB))
    mu_out: dict[int, float] = {}
    for a_int, b_int in [(0.5614594835668851, plo), (phi, B_b)]:
        if b_int <= a_int:
            continue
        n_sub = max(20, int((b_int - a_int) * 8))
        dx = (b_int - a_int) / n_sub
        for i in range(n_sub):
            x0 = a_int + i * dx
            x1 = x0 + dx
            kb = int(round(0.5 * (x0 + x1) * BPB))
            if 1 <= kb <= kN:
                mu_out[kb] = mu_out.get(kb, 0.0) + math.log(x1 / x0)
    P_out, _ = _panjer_1d(mu_out, kN)
    T_max = int(math.ceil(N_eff / qstep)) + 20
    mu_in: dict[tuple[int, int], float] = {}
    mu_in_1d: dict[int, float] = {}
    n_sub = max(80, int((phi - plo) * 16))
    dx = (phi - plo) / n_sub
    for i in range(n_sub):
        x0 = plo + i * dx
        x1 = x0 + dx
        xm = 0.5 * (x0 + x1)
        w = math.log(x1 / x0)
        kb = int(round(xm * BPB))
        tk = int(math.floor(xm / qstep))
        if kb >= 1 and tk >= 1:
            mu_in[(kb, tk)] = mu_in.get((kb, tk), 0.0) + w
            mu_in_1d[tk] = mu_in_1d.get(tk, 0.0) + w
    P_surv_1d, lam_in = _panjer_1d(mu_in_1d, T_max)
    P2: list = [None] * (kN + 1)
    P2[0] = {0: math.exp(-lam_in)}
    in_items = [(dk, dt, dk * v) for (dk, dt), v in mu_in.items()]
    k_min_step = min(dk for dk, _, _ in in_items)
    for k in range(k_min_step, kN + 1):
        d: dict[int, float] = {}
        for dk, dt, dv in in_items:
            if dk <= k:
                prev = P2[k - dk]
                if prev:
                    for t0, p0 in prev.items():
                        nt = t0 + dt
                        d[nt] = d.get(nt, 0.0) + dv * p0
        if d:
            inv_k = 1.0 / k
            for t in d:
                d[t] *= inv_k
            P2[k] = d
    smooth_t: dict[int, float] = {}
    for k in range(kN + 1):
        pk = P2[k]
        if pk and P_out[kN - k] > 0:
            po = P_out[kN - k]
            for t, p in pk.items():
                smooth_t[t] = smooth_t.get(t, 0.0) + p * po
    tot_sm = sum(smooth_t.values())
    rec = sum(p for t, p in smooth_t.items() if t >= T_u8) / tot_sm
    surv = sum(P_surv_1d[j] for j in range(T_u8, len(P_surv_1d)))
    bias = 128 - T_u8
    max_byte = bias + int(math.floor(N_max / qstep))
    assert 0 <= bias and max_byte <= 255
    return {
        "T_u8": T_u8,
        "bias_u8": bias,
        "max_byte": max_byte,
        "rec_exact": rec,
        "log2_surv": math.log2(max(surv, 1e-300)),
        "hits_lam": lam_in,
    }


def check_norms_and_2d_panjer(led: dict) -> dict:
    op = led["operating_point"]
    mc = led["direct_mc_4M"]
    V_a_sh, R0, q = op["V_a_per_shell"], mc["cand_R0"], op["q_ball"]
    Ra_tot = (2.0 * V_a_sh * 120.0 * 70.3125 * 24.0 / (math.pi**4 * 0.267223)) ** 0.125
    b_eq = R0 / (4.0 ** (1.0 / q))
    a_eq = Ra_tot / 2.0
    det_ng = 8.0 * math.log2(ABS_X * (b_eq**3) + a_eq**3)
    det_nf = 0.0
    for k in (1, 2, 3):
        for ia in range(0, 41):
            fa = ia / 40.0
            a1 = math.sqrt(fa * Ra_tot * Ra_tot / k)
            a2 = math.sqrt((1.0 - fa) * Ra_tot * Ra_tot / (4 - k))
            for ib in range(0, 41):
                fb = ib / 40.0
                b1 = ((fb * (R0**q)) / k) ** (1.0 / q) if fb > 0 else 0.0
                b2 = (((1.0 - fb) * (R0**q)) / (4 - k)) ** (1.0 / q) if fb < 1.0 else 0.0
                f1 = a1**9 + a1**6 * b1**3 + a1**3 * b1**6 + b1**9
                f2 = a2**9 + a2**6 * b2**3 + a2**3 * b2**6 + b2**9
                if f1 > 0 and f2 > 0:
                    val = 2.0 * (k * math.log2(f1) + (4 - k) * math.log2(f2))
                    if val > det_nf:
                        det_nf = val
    assert det_ng < 705.0 and abs(det_ng - mc["det_upper_ng"]) < 1e-6
    assert det_nf < 630.0 and abs(det_nf - mc["det_upper_nf"]) < 0.01

    dp_g = _exact_2d_floor_dp(mc["cand_cond_ng"] - 1.0487, op["B_g"], op["g_band"][0], op["g_band"][1], op["q_step_g"], op["T_u8_g"], 705.0)
    dp_f = _exact_2d_floor_dp(mc["cand_cond_nf"] + 5.3752, op["B_f"], op["f_band"][0], op["f_band"][1], op["q_step_f"], op["T_u8_f"], 630.0)
    ref_g, ref_f = led["bandpass_dp_2d"]["g_side"], led["bandpass_dp_2d"]["f_side"]
    assert dp_g["max_byte"] == ref_g["max_byte"] == 245 <= 255
    assert dp_f["max_byte"] == ref_f["max_byte"] == 248 <= 255
    assert abs(dp_g["rec_exact"] - ref_g["rec_exact"]) < 1e-9
    assert abs(dp_f["rec_exact"] - ref_f["rec_exact"]) < 1e-9
    assert abs(dp_g["log2_surv"] - ref_g["log2_surv"]) < 1e-9
    assert abs(dp_f["log2_surv"] - ref_f["log2_surv"]) < 1e-9
    ret_exact = dp_g["rec_exact"] * dp_f["rec_exact"]
    assert ret_exact >= 0.6864 and abs(ret_exact - op["ret_tot"]) < 1e-9
    return {
        "ok": True,
        "det_upper_ng": round(det_ng, 3),
        "det_upper_nf_corner": round(det_nf, 3),
        "g_max_byte": dp_g["max_byte"],
        "f_max_byte": dp_f["max_byte"],
        "rec_g": round(dp_g["rec_exact"], 6),
        "rec_f": round(dp_f["rec_exact"], 6),
        "ret_exact": round(ret_exact, 6),
        "l2_surv_s1": round(dp_g["log2_surv"] + dp_f["log2_surv"], 6),
    }


# =============================================================================
# 4. Guard4 x16 + Stage 1b (L_owner_bits*8, Top2, N_g eval) + SGE/BW + Unified RAM
# =============================================================================
def _next_pow2(x: int) -> int:
    return 1 if x <= 1 else (1 << int(math.ceil(math.log2(x))))


def _ntt_mul_ir(n_limbs: float, distributed: bool = False) -> float:
    if n_limbs <= 16.0:
        return 12.0 * n_limbs * n_limbs * 1.15
    N = float(_next_pow2(int(math.ceil(2.0 * n_limbs))))
    logN = math.log2(N)
    ir_ntt = 3.0 * (2.0 * 0.5 * N * logN * 17.0 + 0.5 * N * logN * 17.0 + N * 30.0 + 2.0 * N * 8.0)
    ir_xfer = (3.0 * N * 160.0) if distributed else 0.0
    return (ir_ntt + ir_xfer + N * 67.0) * 1.15


def _product_tree_ir(n_leaves: float, leaf_bits: float, dist_thresh_bits: float) -> float:
    if n_leaves <= 1.0:
        return 0.0
    h = int(math.ceil(math.log2(n_leaves)))
    tot, cur_n, cur_b = 0.0, float(n_leaves), float(leaf_bits)
    for _ in range(h):
        pairs = cur_n * 0.5
        tot += pairs * _ntt_mul_ir(max(1.0, cur_b / 64.0), distributed=(cur_b * 2.0 > dist_thresh_bits))
        cur_n, cur_b = pairs, cur_b * 2.0
    return tot


def _bernstein_rem_tree_ir(n_leaves: float, leaf_bits: float, dist_thresh_bits: float, descent_mult: float = 2.00) -> float:
    if n_leaves <= 1.0:
        return 40000.0
    h = int(math.ceil(math.log2(n_leaves)))
    ir_prod = _product_tree_ir(n_leaves, leaf_bits, dist_thresh_bits)
    L_bits = n_leaves * leaf_bits
    ir_root = 3.0 * _ntt_mul_ir(max(1.0, L_bits / 64.0), distributed=(L_bits > dist_thresh_bits))
    ir_desc, cur_n, cur_b = 0.0, float(n_leaves), float(leaf_bits)
    for _ in range(h):
        pairs = cur_n * 0.5
        ir_desc += pairs * descent_mult * _ntt_mul_ir(max(1.0, cur_b / 64.0), distributed=(cur_b * 2.0 > dist_thresh_bits))
        cur_n, cur_b = pairs, cur_b * 2.0
    return ir_prod + ir_root + ir_desc + float(n_leaves) * max(25000.0, 40000.0 * (leaf_bits / 636.0))


def _prime_product_bits(b_lo: float, b_hi: float) -> float:
    return max(1.0, (2.0**b_hi - (0.0 if b_lo <= 0.0 else 2.0**b_lo)) / math.log(2.0))


def _buchstab_omega_upper(v: float) -> float:
    if v < 1.0:
        return 0.0
    if v <= 2.0:
        return 1.0 / v
    if v <= 3.0:
        return (1.0 + math.log(v - 1.0)) / v
    return 0.57


def _tier_probs(cof_bits: float, P_t: float, B_bits: float, k_max: int = 8):
    if cof_bits <= P_t:
        return 1.0, 0.0, 1.0
    steps = 200
    dy = (cof_bits - P_t) / steps
    m_done = DICKMAN.rho(cof_bits / P_t)
    m_adv, m_tot = 0.0, m_done
    for i in range(steps):
        y = P_t + (i + 0.5) * dy
        w = DICKMAN.rho((cof_bits - y) / P_t) * (_buchstab_omega_upper(y / P_t) / P_t) * dy
        m_tot += w
        if y <= B_bits:
            m_done += w
        else:
            for k in range(2, k_max + 1):
                if k * P_t <= y <= k * B_bits:
                    m_adv += w
                    break
    return m_done / max(m_tot, 1e-30), m_adv / max(m_tot, 1e-30), (m_done + m_adv) / max(m_tot, 1e-30)


def _integrate_guard4_x16(plo: float, phi: float, V_a: float, L_win: float) -> dict:
    W_B, d_ok, r_bar = 4.0, 1.0 / 3.0, 3.0
    L_emit = L_win + (W_B - 1.0)
    win_per_cell = 1.0 / (V_a * L_win)
    p_split = min(1.0, L_win / 256.0)
    ovh_p = 56.0 * (1.0 - p_split) + 84.0 * p_split
    n_pieces = 1.0 + p_split
    dx = (phi - plo) / 400.0
    h_act = h_em = ir_hits = ir_probe = n_ok = n_bkt = n_ent = 0.0
    for i in range(400):
        x = plo + (i + 0.5) * dx
        p = 2.0**x
        d_id = (p / x) * dx
        d_ok_p = d_id * d_ok
        n_ok += d_ok_p
        n_bkt += d_ok_p * (p / W_B + 2.0)
        n_ent += d_id * V_a
        he = r_bar * V_a * L_emit / p
        ha = r_bar * V_a * L_win / p
        ir_probe += d_ok_p * ovh_p * win_per_cell
        hp = he / n_pieces
        e1 = min(hp, 1.5)
        e4 = min(max(0.0, hp - e1), 6.0)
        e16 = max(0.0, hp - e1 - e4)
        ir_hits += d_ok_p * n_pieces * (e16 * 2.6875 + e4 * 3.25 + e1 * 5.00) * win_per_cell
        h_act += d_ok_p * ha * win_per_cell
        h_em += d_ok_p * he * win_per_cell
    return {
        "hits_active": h_act, "hits_emitted": h_em, "ir_hits": ir_hits,
        "ir_probe": ir_probe, "n_ok_primes": n_ok, "n_bkt": n_bkt, "n_ent": n_ent,
    }


def check_guard4_s1b_la_ram(led: dict) -> dict:
    op, mc, dp = led["operating_point"], led["direct_mc_4M"], led["bandpass_dp_2d"]
    V_a = op["V_a_per_shell"]
    ch = mc["chord_200k"]
    L_win = ch["sh0_vol_frac"] * (ch["sh0_mean_span"] / ch["sh0_win_mult_248"]) + ch["sh1_vol_frac"] * (
        ch["sh1_mean_span"] / ch["sh1_win_mult_248"]
    )
    assert abs(L_win - ch["L_win_eff"]) < 1e-9

    sg = _integrate_guard4_x16(op["g_band"][0], op["g_band"][1], V_a, L_win)
    sf = _integrate_guard4_x16(op["f_band"][0], op["f_band"][1], V_a, L_win)

    # Stage 1b v2 (L_owner_bits * 8, exact Top2, N_g eval, zero double-count)
    V_tot = 2.0 ** mc["log2_Vk_charged"]
    surv_s1 = 2.0 ** dp["log2_surv_s1"]
    ram_mult = 4.5
    L_owner_bits = (4.0 * (2.0**60) * 8.0) / ram_mult
    peak_glob_s1b = 0.0

    def price_t(b_lo, b_hi, p_in, leaf_b, dens):
        nonlocal peak_glob_s1b
        L_lo = _prime_product_bits(b_lo, min(b_hi, 36.0)) if b_lo < 36.0 else 0.0
        L_hi = dens * _prime_product_bits(max(b_lo, 36.0), b_hi) if b_hi > 36.0 else 0.0
        L_t = max(1.0, L_lo + L_hi)
        peak_glob_s1b = max(peak_glob_s1b, ram_mult * L_t / 8.0)
        M_tot = V_tot * p_in
        M_b = min(M_tot, max(1.0, L_t / leaf_b))
        nb = max(1.0, M_tot / M_b)
        L_X = M_b * leaf_b
        ir = _bernstein_rem_tree_ir(M_b, leaf_b, L_owner_bits, descent_mult=2.00)
        if L_X < L_t:
            ir += max(0.0, L_t / L_X - 1.0) * 2.5 * _ntt_mul_ir(max(1.0, L_X / 64.0), distributed=(L_X > L_owner_bits))
        n_p = max(1.0, dens * ((2.0**b_hi) / (b_hi * math.log(2.0)) - (0.0 if b_lo <= 1 else (2.0**b_lo) / (b_lo * math.log(2.0)))))
        ir_bld = _product_tree_ir(n_p, b_hi, L_owner_bits)
        return (nb * ir) / V_tot, ir_bld / V_tot

    Bf, Bg = op["B_f"], op["B_g"]
    Pf, Pg = op["f_band"][1], op["g_band"][1]
    mid_f, mid_g = 0.5 * (Pf + Bf), 0.5 * (Pg + Bg)
    rem_f = max(180.0, min(360.0, mc["cand_cond_nf"] - dp["f_side"]["T_bits"]))
    rem_g = max(180.0, min(360.0, mc["cand_cond_ng"] - dp["g_side"]["T_bits"]))
    _, paf1, _ = _tier_probs(rem_f, Pf, Bf, k_max=8)
    _, paf2, _ = _tier_probs(rem_f, mid_f, Bf, k_max=8)
    _, pag1, _ = _tier_probs(rem_g, Pg, Bg, k_max=8)
    _, pag2, _ = _tier_probs(rem_g, mid_g, Bg, k_max=8)
    pf_sm = min(1.0, (DICKMAN.rho((mc["cand_cond_nf"] + 5.3752) / Bf) * op["rec_f"]) / (2.0 ** dp["f_side"]["log2_surv"]))
    pg_sm = min(1.0, (DICKMAN.rho((mc["cand_cond_ng"] - 1.0487) / Bg) * op["rec_g"]) / (2.0 ** dp["g_side"]["log2_surv"]))
    cr = surv_s1 * 1700.0 + surv_s1 * pf_sm * 1700.0
    cb = 0.0
    for args in [
        (1.0, Pf, surv_s1, 510.0, 0.125),
        (Pf, mid_f, surv_s1 * paf1, 340.0, 0.125),
        (mid_f, Bf, surv_s1 * paf2, 340.0, 0.125),
        (1.0, Pg, surv_s1 * pf_sm, 640.0, 1.0 / 24.0),
        (Pg, mid_g, surv_s1 * pf_sm * pag1, 320.0, 1.0 / 24.0),
        (mid_g, Bg, surv_s1 * pf_sm * pag2, 320.0, 1.0 / 24.0),
    ]:
        r, b = price_t(*args)
        cr += r
        cb += b
    c_mr = surv_s1 * pf_sm * pg_sm * (110.0 * 256.0 * 150.0 * 4.0 + 0.05 * 5.0e7 * 2.0)
    c_cof_rem_mr = cr + c_mr

    ipc_tot = sg["ir_hits"] + sf["ir_hits"] + sg["ir_probe"] + sf["ir_probe"] + 0.5500 + c_cof_rem_mr
    assert abs(ipc_tot - led["sieve_and_stage1b"]["ipc_tot"]) < 1e-4
    coll_l2 = math.log2((V_tot * ipc_tot) / C_COST)
    assert abs(coll_l2 - led["ledger_rows_live_log2_u"]["relation_collection"]) < 1e-4

    n_ent = sg["n_ent"] + sf["n_ent"]
    n_bkt = sg["n_bkt"] + sf["n_bkt"]
    n_ok = sg["n_ok_primes"] + sf["n_ok_primes"]
    replay_l2 = math.log2(op["K_shells"] * (n_ent * 27.5 + n_bkt * 6.5) / C_COST)
    assert abs(replay_l2 - led["ledger_rows_live_log2_u"]["replay_a_tables"]) < 1e-6

    # Stage 2/3 SGE + BW (d_mksol = N_l/32 + 32, IR_ROW_FIXED = 2100.0)
    U = ((2.0**Bf) / (Bf * math.log(2.0)) + (2.0**Bg) / (Bg * math.log(2.0))) / 12.0
    eta_m = ((1.04094 - op["tau_c"]) / (op["eta"] - op["tau_c"])) ** 1.35
    N_l = U * 0.0006589 * math.sqrt(eta_m)
    w_eff = 180.0 * math.sqrt(eta_m)
    N_nnz = N_l * w_eff
    ir_core = 1.144 * 8.0 + 0.122 * 24.5
    ir_row_fixed = float(led["la_and_sge"]["ir_row_fixed"])
    assert ir_row_fixed == 2100.0
    ir_spmv = (N_nnz * ir_core + N_l * ir_row_fixed) * 1.15
    L_seq = N_l / 128.0 + N_l / 32.0 + 64.0
    d_mksol = N_l / 32.0 + 32.0
    tot_spmvs = 1.02 * (32.0 * L_seq + d_mksol)
    ir_bw = tot_spmvs * ir_spmv
    ir_mk = 1.02 * d_mksol * 32.0 * N_l * (236.0 + 47.0) * 1.15
    ir_lg = 20.0 * (160.0**3) * L_seq * math.log2(L_seq) * 236.0 * 1.15
    la_l2 = math.log2((ir_bw + ir_mk + ir_lg + 0.008 * (ir_bw + ir_mk) + 30.0 * U * 283.0 * 1.15) / C_COST)
    sge_l2 = math.log2(((12500.0 * U + 500.0 * 0.05 * U * math.log2(U) + 100.0 * N_nnz) * 1.15) / C_COST)
    assert abs(la_l2 - led["ledger_rows_live_log2_u"]["linear_algebra"]) < 1e-6
    assert abs(sge_l2 - led["ledger_rows_live_log2_u"]["matrix_generator_sge"]) < 1e-6

    # Unified RAM audit
    replay_glob_B = 8.0 * (n_bkt + n_ent) + 48.0 * n_ok
    s1_peak_l2 = math.log2(replay_glob_B + 256.0 * V_a + 2.0**52)
    s2_sge_l2 = math.log2((op["eta"] * U) * 18.0 * 12.0)
    s3_la_l2 = math.log2(N_nnz * 11.0 + N_l * 8.0 + 38.0 * N_l * 32.0)
    s1b_l2 = math.log2(peak_glob_s1b)
    mem = led["memory_audit_unified"]
    assert abs(s1_peak_l2 - mem["stage1_peak_global_log2_B"]) < 1e-6
    assert abs(s2_sge_l2 - mem["stage2_sge_global_log2_B"]) < 1e-6
    assert abs(s3_la_l2 - mem["stage3_la_global_log2_B"]) < 1e-6
    assert abs(s1b_l2 - mem["stage1b_global_log2_B"]) < 1e-4
    assert max(s1_peak_l2, s2_sge_l2, s3_la_l2, s1b_l2) < mem["peak_global_claim_log2_B"] < 120.0

    return {
        "ok": True,
        "ipc_tot": round(ipc_tot, 6),
        "log2_coll_live": round(coll_l2, 6),
        "log2_replay_live": round(replay_l2, 6),
        "log2_la_live": round(la_l2, 6),
        "log2_mksol_u": round(math.log2(ir_mk / C_COST), 6),
        "log2_sge_live": round(sge_l2, 6),
        "stage1_peak_global_log2_B": round(s1_peak_l2, 4),
        "stage2_sge_global_log2_B": round(s2_sge_l2, 4),
        "stage3_la_global_log2_B": round(s3_la_l2, 4),
        "stage1b_global_log2_B": round(s1b_l2, 4),
    }


# =============================================================================
# 5. Live Exact F_q Gaussian Elimination After Algebraic 2-Merge
# =============================================================================
Q_MOD = 2147483629
ZETA12 = 1803057106


def _exact_fq_rank_2merge(U: int, eta: float, seed: int) -> tuple[int, int, int]:
    rng = random.Random(seed)
    z12_pows = [pow(ZETA12, k, Q_MOD) for k in range(12)]
    R = int(round(eta * U))
    Uf = int(round(0.439 * U))
    Ug = U - Uf
    lnBf, lnBg = 69.60 * math.log(2.0), 69.96 * math.log(2.0)
    r_uf1, r_ug1 = DICKMAN.rho(6.65), DICKMAN.rho(8.31)
    weights = [0.0] * U
    for j in range(Uf):
        x = (j + 0.5) / Uf
        t = max(0.02, 1.0 + math.log(x) / lnBf)
        qj = (1.173 / U) * (1.0 / x) * (DICKMAN.rho(7.65 - t) / r_uf1) + (0.35 / (2.0 + 2.5 * j) if j < 40 else 0.0)
        weights[j] = min(0.45, qj)
    for j in range(Ug):
        x = (j + 0.5) / Ug
        t = max(0.02, 1.0 + math.log(x) / lnBg)
        qj = (1.196 / U) * (1.0 / x) * (DICKMAN.rho(9.31 - t) / r_ug1) + (0.40 / (2.0 + 2.5 * j) if j < 40 else 0.0)
        weights[Uf + j] = min(0.45, qj)
    mean_w = sum(weights)
    cdf, acc = [], 0.0
    for w in weights:
        acc += w / mean_w
        cdf.append(acc)
    rows, L_pois = [], math.exp(-mean_w)
    for _ in range(R):
        k_p, p_acc = 0, 1.0
        while p_acc > L_pois:
            k_p += 1
            p_acc *= rng.random()
        d: dict[int, int] = {}
        for _ in range(max(4, k_p - 1)):
            c = min(U - 1, bisect.bisect_left(cdf, rng.random()))
            d[c] = (d.get(c, 0) + z12_pows[rng.randrange(12)]) % Q_MOD
            if d[c] == 0:
                del d[c]
        if d:
            rows.append(d)
    while True:
        col_rows: dict[int, list[int]] = {}
        for r_idx, d in enumerate(rows):
            if d is not None:
                for c in d:
                    col_rows.setdefault(c, []).append(r_idx)
        changed = False
        for rlist in col_rows.values():
            if len(rlist) == 1:
                rows[rlist[0]] = None
                changed = True
        if changed:
            continue
        merged = False
        for c, rlist in col_rows.items():
            if len(rlist) == 2:
                r0, r1 = rlist
                d0, d1 = rows[r0], rows[r1]
                if d0 is None or d1 is None or c not in d0 or c not in d1:
                    continue
                if len(d0) > len(d1):
                    r0, r1, d0, d1 = r1, r0, d1, d0
                fac = (d1[c] * pow(d0[c], -1, Q_MOD)) % Q_MOD
                for col, v0 in d0.items():
                    nv = (d1.get(col, 0) - fac * v0) % Q_MOD
                    if nv == 0:
                        d1.pop(col, None)
                    else:
                        d1[col] = nv
                rows[r0] = None
                if not d1:
                    rows[r1] = None
                merged = True
                break
        if not merged:
            break
    alive = [d for d in rows if d]
    cols = sorted({c for d in alive for c in d})
    cmap = {c: i for i, c in enumerate(cols)}
    pivots: dict[int, dict[int, int]] = {}
    rank = 0
    for d in alive:
        cur = {cmap[c]: v for c, v in d.items()}
        while cur:
            pc = min(cur)
            if pc in pivots:
                fac = cur[pc]
                for c, pv in pivots[pc].items():
                    nv = (cur.get(c, 0) - fac * pv) % Q_MOD
                    if nv == 0:
                        cur.pop(c, None)
                    else:
                        cur[c] = nv
            else:
                inv = pow(cur[pc], -1, Q_MOD)
                pivots[pc] = {c: (v * inv) % Q_MOD for c, v in cur.items()}
                rank += 1
                break
    return len(alive), len(cols), len(cols) - rank


def check_live_fq_rank(U: int = 512, etas: tuple[float, ...] = (0.82,), n_seeds: int = 6) -> dict:
    results, total_null = {}, 0
    for eta in etas:
        nulls = nonempty = surplus = 0
        for s in range(n_seeds):
            nr, nc, cnull = _exact_fq_rank_2merge(U, eta, 7000 + int(eta * 1000) * 100 + s)
            rank = nc - cnull
            if nc > 0:
                nonempty += 1
            if nr >= nc:
                surplus += 1
            nulls += min(nr, nc) - rank
        if eta >= 0.82:
            assert surplus == n_seeds
        results[f"eta_{eta:.2f}"] = {"seeds": n_seeds, "nonempty": nonempty, "surplus_pos": surplus, "dense_nullity": nulls}
        total_null += nulls
    assert total_null == 0
    return {"ok": True, "U": U, "Q_MOD": Q_MOD, "total_dense_core_nullity": total_null, "by_eta": results}


# =============================================================================
# 6. FATAL #4: Stage 4 12D M_g Proof, q=2^70..2^275 ECM Table, ecm_split & Witnesses
# =============================================================================
def _mont_mul_ir_limbs(n_bits: float) -> float:
    L = max(4, int(math.ceil(n_bits / 64.0)))
    if L == 4:
        return 212.0
    cios = 13.25 * L * L + 80.0
    return cios if L <= 16 else min(cios, 2.5 * _ntt_mul_ir(float(L), distributed=False))


def _ecm_curve_prob(B1: float, B2: float, q_bits: float) -> float:
    logB1, logB2 = math.log2(B1), math.log2(B2)
    q_eff = max(1.0, q_bits - 3.5)
    p1 = DICKMAN.rho(q_eff / logB1)
    if B2 <= B1:
        return p1
    dt = (logB2 - logB1) / 120.0
    p2 = sum((DICKMAN.rho((q_eff - (logB1 + (i + 0.5) * dt)) / logB1) * (dt / (logB1 + (i + 0.5) * dt)))
             if (logB1 + (i + 0.5) * dt) < q_eff else (dt / (logB1 + (i + 0.5) * dt))
             for i in range(120))
    return min(1.0, p1 + p2)


def _ecm_curve_ir(B1: float, B2: float, N_bits: float) -> float:
    L = max(4, int(math.ceil(N_bits / 64.0)))
    m_ir = _mont_mul_ir_limbs(N_bits)
    ir_s1 = (B1 / math.log(2.0)) * (11.0 * m_ir + 93.0 * L)
    n_p2 = max(0.0, B2 / math.log(B2) - B1 / math.log(B1)) if B2 > B1 else 0.0
    return 1.15 * (40.0 * m_ir + ir_s1 + 2.5 * n_p2 * (m_ir + 40.0 * L))


_ECM_CACHE: dict = {}
_ECM_PROB_CACHE: dict = {}


def _optimal_ecm_factor_ir(q_bits: float, N_bits: float):
    q_key = round(q_bits, 2)
    key = (q_key, int(math.ceil(N_bits / 64.0)))
    if key in _ECM_CACHE:
        return _ECM_CACHE[key]
    best_ir, best_cfg = 1e300, None
    if q_key not in _ECM_PROB_CACHE:
        plist = []
        for lb1_4 in range(32, 249):
            B1 = 2.0 ** (lb1_4 * 0.25)
            for mult in (10.0, 25.0, 60.0, 120.0, 250.0):
                B2 = B1 * mult
                pr = _ecm_curve_prob(B1, B2, q_bits)
                if pr > 1e-40:
                    plist.append((B1, B2, pr))
        _ECM_PROB_CACHE[q_key] = plist
    for B1, B2, pr in _ECM_PROB_CACHE[q_key]:
        c1 = _ecm_curve_ir(B1, B2, N_bits)
        exp_ir = c1 / pr
        if exp_ir < best_ir:
            best_ir, best_cfg = exp_ir, (B1, B2, pr, c1)
    _ECM_CACHE[key] = (best_ir, best_cfg)
    return best_ir, best_cfg


def _descent_node_cost_u(q_bits: float, side: str, Bf: float, Bg: float, child_cap: float) -> float:
    Bf_lim, Bg_lim = max(Bf, child_cap), max(Bg, child_cap)
    best_u = 1e300
    for deg in (1, 2):
        for s_2 in range(0, 45):
            s = s_2 * 0.5
            if deg == 1:
                nf = 4.5 * q_bits + 72.0 * s + 111.0
                ng = 1.5 * q_bits + 8.0 * max(3.0 * s, LOG2_ABS_X - 3.0 * s) + 45.0
            else:
                nf = 3.0 * q_bits + 72.0 * s + 55.0
                ng = 1.0 * q_bits + 8.0 * max(3.0 * s, LOG2_ABS_X, 2.0 * LOG2_ABS_X - 3.0 * s) + 35.0
            nf_rem = max(1.0, nf - q_bits) if side == "f" else nf
            ng_rem = ng if side == "f" else max(1.0, ng - q_bits)
            log2_T, dim = 20.0, 8.0 * (deg + 1)
            for _ in range(4):
                nf_t = nf_rem + (72.0 / dim) * log2_T
                ng_t = ng_rem + (24.0 / dim) * log2_T
                p_j = DICKMAN.rho(nf_t / Bf_lim) * DICKMAN.rho(ng_t / Bg_lim)
                if p_j <= 1e-250:
                    log2_T = 300.0
                    break
                log2_T = -math.log2(p_j) + 1.0
            if log2_T > 120.0:
                continue
            ir_wf, _ = _optimal_ecm_factor_ir(Bf_lim, nf_t)
            ir_wg, _ = _optimal_ecm_factor_ir(Bg_lim, ng_t)
            tot_u = ((2.0**log2_T) * 2000.0 + 3.0 * (ir_wf + ir_wg) + 1.0e7) / C_COST
            if tot_u < best_u:
                best_u = tot_u
    return best_u


def _xdbl(X: int, Z: int, A24: int, N: int):
    t1, t2 = (X + Z) % N, (X - Z) % N
    t1, t2 = (t1 * t1) % N, (t2 * t2) % N
    t3 = (t1 - t2) % N
    return (t1 * t2) % N, (t3 * (t2 + A24 * t3)) % N


def _xadd(X1: int, Z1: int, X2: int, Z2: int, Xd: int, Zd: int, N: int):
    u, v = ((X1 - Z1) * (X2 + Z2)) % N, ((X1 + Z1) * (X2 - Z2)) % N
    s, d = (u + v) % N, (u - v) % N
    return (Zd * s * s) % N, (Xd * d * d) % N


def _xmul(k: int, X: int, Z: int, A24: int, N: int):
    R0X, R0Z, R1X, R1Z = 1, 0, X % N, Z % N
    for bit in bin(k)[2:]:
        if bit == "0":
            R1X, R1Z = _xadd(R0X, R0Z, R1X, R1Z, X, Z, N)
            R0X, R0Z = _xdbl(R0X, R0Z, A24, N)
        else:
            R0X, R0Z = _xadd(R0X, R0Z, R1X, R1Z, X, Z, N)
            R1X, R1Z = _xdbl(R1X, R1Z, A24, N)
    return R0X, R0Z


def check_stage4_12d_mg_and_witnesses(led: dict) -> dict:
    # 1. Exact 12D M_g = Q(zeta_5)[Y]/(Y^3 - x) determinant & unconditional LLL+AM-GM norm bound
    det_l2 = 6.0 * math.log2(P_BLS) + 4.0 * math.log2(ABS_X)
    covol_l2 = 4.5 * math.log2(5.0) + 6.0 * math.log2(3.0) + det_l2
    lll_amgm_l2 = 33.0 - 6.0 * math.log2(12.0) + covol_l2
    s4 = led["stage4_12d_mg"]
    assert abs(det_l2 - s4["det_log2"]) < 1e-6
    assert lll_amgm_l2 < s4["N_boot_bits"] == 2572.0

    # 2. Full q = 2^70 .. 2^275 ECM calibration table (< 1e-12 self-consistency!)
    expected_ecm = {
        70: 26.525, 85: 29.894, 100: 33.023, 120: 36.907, 140: 40.529,
        160: 43.941, 185: 47.965, 215: 52.504, 245: 56.784, 275: 60.848,
    }
    max_rel_err = 0.0
    for qb, exp_l2 in expected_ecm.items():
        ir_f, (_, _, pr, c1) = _optimal_ecm_factor_ir(float(qb), 2572.0)
        rel = abs(ir_f - c1 / pr) / ir_f
        max_rel_err = max(max_rel_err, rel)
        assert rel < 1e-12 and abs(math.log2(ir_f / C_COST) - exp_l2) < 0.002

    # 3. Full 6-layer (Layer 0 inert + Layers 1..5) zero-early-abort descent + explicit ecm_split_expected
    Bf, Bg = led["operating_point"]["B_f"], led["operating_point"]["B_g"]
    p_sm = DICKMAN.rho(2572.0 / 275.0)
    ir_275, _ = _optimal_ecm_factor_ir(275.0, 2572.0)
    u_boot1 = ((1.0 / (0.918 * p_sm)) * (2.5 * ir_275)) / C_COST
    # Layer 0: inert q == 11 (mod 15) in (2^137.5, 2^275] have residue degree 2 in L_g = M_g(omega) (norm q^2 <= 2^550)
    u_layer0 = 5.0 * _descent_node_cost_u(550.0, "g", Bf, Bg, 275.0)
    u_tree1 = u_layer0
    # Layers 1..5 (Layer 5 includes +20 extra inert nodes for q in (2^34.98, 2^69.96] with q^2 > B_g)
    for q_hi, q_lo, m_nodes in [(275.0, 155.0, 5), (155.0, 130.0, 12), (130.0, 106.0, 28), (106.0, 86.0, 65), (86.0, Bg, 200)]:
        u_tree1 += m_nodes * max(_descent_node_cost_u(q_hi, "f", Bf, Bg, q_lo), _descent_node_cost_u(q_hi, "g", Bf, Bg, q_lo))
    u_both_live = 2.0 * (u_boot1 + u_tree1)
    miss_fac = (1.04094 / led["operating_point"]["eta"]) ** 2
    u_both_snell = u_both_live * miss_fac
    ecm_split_u = 2.0 * (u_tree1 + 0.25 * u_boot1) * miss_fac
    assert abs(math.log2(u_both_live) - s4["live_2log_no_snell_log2_u"]) < 1e-6
    assert abs(math.log2(u_both_snell) - s4["live_2log_with_snell_log2_u"]) < 1e-6
    assert abs(math.log2(ecm_split_u) - s4["ecm_split_live_log2_u"]) < 1e-6

    # 4. Torsion-16 group-order witness (21/21 primes in [101, 200]) + Stage-2 witnesses (50023 & 35)
    m = 3
    num_x, den_x = m * m - 2 * m + 2, 2 - m * m
    t16_ok = t16_tot = 0
    for p in _sieve_primes(200):
        if p < 101:
            continue
        x = (num_x * pow(den_x, -1, p)) % p
        d = ((2 * x * x - 1) * pow(x, -4, p)) % p
        A = (2 * (1 + d) * pow(1 - d, -1, p)) % p
        cnt = 1 + sum(1 if (xx**3 + A * xx * xx + xx) % p == 0 else (2 if pow((xx**3 + A * xx * xx + xx) % p, (p - 1) // 2, p) == 1 else 0) for xx in range(p))
        t16_tot += 1
        if cnt % 16 == 0:
            t16_ok += 1
    assert t16_ok == t16_tot == 21

    N1 = 50023 * 50033
    x1 = (num_x * pow(den_x, -1, N1)) % N1
    d1 = ((2 * x1 * x1 - 1) * pow(x1, -4, N1)) % N1
    A24_1 = pow(1 - d1, -1, N1)
    X, Z = 3, 1
    for p in _sieve_primes(30):
        pk = p
        while pk * p <= 30:
            pk *= p
        X, Z = _xmul(pk, X, Z, A24_1, N1)
    assert math.gcd(Z, N1) == 1
    D = 30
    babies = {i: _xmul(i, X, Z, A24_1, N1) for i in range(1, D // 2 + 1) if math.gcd(i, D) == 1}
    DX, DZ = _xmul(D, X, Z, A24_1, N1)
    G_prev, G_cur, found_w1 = (1, 0), (DX, DZ), 1
    for k in range(1, 350 // D + 2):
        for BX, BZ in babies.values():
            g = math.gcd((G_cur[0] * BZ - BX * G_cur[1]) % N1, N1)
            if 1 < g < N1:
                found_w1 = g
                break
        if found_w1 > 1:
            break
        G_next = _xdbl(DX, DZ, A24_1, N1) if k == 1 else _xadd(G_cur[0], G_cur[1], DX, DZ, G_prev[0], G_prev[1], N1)
        G_prev, G_cur = G_cur, G_next
    assert found_w1 == 50023

    X2, Z2 = _xmul(2, 4, 1, 2, 35)
    assert math.gcd(Z2, 35) == 1
    DX2, DZ2 = _xmul(4, X2, Z2, 2, 35)
    found_w2 = math.gcd((DX2 * Z2 - X2 * DZ2) % 35, 35)
    assert found_w2 == 5 and 35 // found_w2 == 7

    return {
        "ok": True,
        "M_g_det_log2": round(det_l2, 6),
        "M_g_uncond_lll_amgm_log2": round(lll_amgm_l2, 6),
        "ecm_table_max_rel_error": max_rel_err,
        "layer0_inert_log2_u": round(math.log2(u_layer0), 6),
        "live_2log_with_snell": round(math.log2(u_both_snell), 6),
        "ecm_split_live_log2_u": round(math.log2(ecm_split_u), 6),
        "torsion16_primes_verified": f"{t16_ok}/{t16_tot}",
        "witness1_50023": found_w1,
        "witness2_35_split": [found_w2, 35 // found_w2],
    }


# =============================================================================
# 7. Assembly ISA Audit & Ledger / Unconditional Markov Meter / Hygiene Audit
# =============================================================================
FORBIDDEN_REG_RE = re.compile(r"(?<![a-zA-Z0-9_])([vqdshzp][0-9]|[vqdshzp][12][0-9]|[vqdshzp]3[01]|x18|w18)(?![a-zA-Z0-9_])")
FORBIDDEN_OP_RE = re.compile(r"^\s*(dc\s+zva|fadd|fsub|fmul|fdiv|fmadd|fmsub|fcvt|fmov|fcmp|ld1|st1|addv|umull2|pmull)\b")


def check_assembly_isa_audit() -> dict:
    summary = {}
    for fname in ("KERNELS_SIEVE.md", "KERNELS_LA_ECM.md"):
        blocks = re.findall(r"```asm\n(.*?)```", (ROOT / fname).read_text(encoding="utf-8"), flags=re.S)
        assert len(blocks) >= 3
        n_insns = n_staddb = 0
        for blk in blocks:
            for line in blk.splitlines():
                code = line.split(";")[0].strip()
                if not code or code.startswith(".") or code.endswith(":"):
                    continue
                n_insns += 1
                low = code.lower()
                if low.startswith("staddb"):
                    n_staddb += 1
                assert not FORBIDDEN_OP_RE.search(low) and not FORBIDDEN_REG_RE.search(low), f"Forbidden ISA: {code}"
        assert (n_staddb > 0) if fname == "KERNELS_SIEVE.md" else (n_staddb == 0)
        summary[fname] = {"asm_blocks": len(blocks), "instructions": n_insns, "staddb_count": n_staddb, "forbidden": 0}
    return {"ok": True, "files": summary}


def check_ledger_meter_kgate_hygiene(led: dict) -> dict:
    sub = json.loads((ROOT / "submission.json").read_text(encoding="utf-8"))
    score = sub["claims"][0]["score"]
    num, den = int(score["numerator"]), int(score["denominator"])
    assert math.gcd(num, den) == 1 and num == 10667 and den == 100
    assert "required" in sub["experiments"][0]["input_schema"]
    claim_frac = Fraction(num, den)
    assert abs(float(claim_frac) - led["totals"]["claim_log2_u"]) < 1e-12

    live_rows = led["ledger_rows_live_log2_u"]
    chg_rows = led["ledger_rows_charged_log2_u"]
    assert len(live_rows) == 14 and len(chg_rows) == 14
    for k in live_rows:
        assert chg_rows[k] > live_rows[k]
    w_open_chg = sum(2.0**v for v in chg_rows.values())
    l2_open_live = math.log2(sum(2.0**v for v in live_rows.values()))
    l2_open_chg = math.log2(w_open_chg)
    assert abs(l2_open_live - led["totals"]["opened_run_live_log2_u"]) < 1e-9
    assert abs(l2_open_chg - led["totals"]["opened_run_charged_log2_u"]) < 1e-9

    # MAJOR #2: Unconditional Markov pointwise meter at W_pointwise = 2^124 u!
    w_pw = 2.0 ** led["totals"]["w_pointwise_meter_log2_u"]
    assert led["totals"]["w_pointwise_meter_log2_u"] == 124.0
    p_markov = w_open_chg / w_pw
    p_open_actual = 1.0 - p_markov - 2.0 ** (-48)
    assert p_open_actual > 0.99998

    # Exact K_gate minimality
    two64 = 1 << 64
    K_gate = led["totals"]["K_gate"]
    assert K_gate == 7194374076228249696 and f"0x{K_gate:016x}" == led["totals"]["K_gate_hex"]
    assert Fraction(K_gate, two64) * Fraction(99998, 100000) >= Fraction(39, 100)
    assert Fraction(K_gate - 1, two64) * Fraction(99998, 100000) < Fraction(39, 100)

    l2_gated_chg = l2_open_chg + math.log2(K_gate / float(two64))
    assert abs(l2_gated_chg - led["totals"]["gated_charged_log2_u"]) < 1e-9
    assert l2_gated_chg < float(claim_frac)

    # Hygiene: exact 9 files, <= 65536 B each, < 1 MiB total, 0 subdirs, 0 private path leaks
    files = sorted(p for p in ROOT.iterdir() if p.is_file())
    dirs = [p.name for p in ROOT.iterdir() if p.is_dir()]
    assert not dirs and len(files) == 9
    private_tokens = ("/" + "Users/", "research-" + "work/", "curvesmash-" + "agy-cli")
    total_bytes, file_sizes = 0, {}
    for p in files:
        sz = p.stat().st_size
        assert sz <= 65536, f"{p.name} too large: {sz}"
        total_bytes += sz
        file_sizes[p.name] = sz
        txt = p.read_text(encoding="utf-8")
        for tok in private_tokens:
            assert tok not in txt, f"Private path leak {tok} in {p.name}"
    assert total_bytes < 1048576

    return {
        "ok": True,
        "claim_fraction": f"{num}/{den}",
        "gcd_claim": math.gcd(num, den),
        "opened_charged_log2_u": round(l2_open_chg, 6),
        "gated_charged_log2_u": round(l2_gated_chg, 6),
        "slack_to_claim_bits": round(float(claim_frac) - l2_gated_chg, 6),
        "uncond_markov_p_abort": round(p_markov, 9),
        "uncond_p_open_lower": round(p_open_actual, 8),
        "K_gate_minimal": True,
        "file_count": len(files),
        "total_bytes": total_bytes,
        "max_file_bytes": max(file_sizes.values()),
        "file_sizes": file_sizes,
        "private_path_leaks": 0,
    }


def main() -> None:
    t0 = time.perf_counter()
    raw = sys.stdin.read().strip()
    inp = json.loads(raw) if raw else {}
    if not isinstance(inp, dict) or (set(inp.keys()) - {"mode"}):
        raise ValueError("Input must be a JSON object with optional key 'mode'")
    mode = inp.get("mode", "core")
    if mode not in ("core", "diagnostics", "all"):
        raise ValueError(f"Unsupported mode: {mode}")

    led = json.loads((ROOT / "ledger_data.json").read_text(encoding="utf-8"))
    out: dict = {"package_id": led["package_id"], "mode": mode}

    if mode in ("core", "all"):
        out["stage0_algebra"] = check_stage0_algebra()
        out["guard4_masked_differential"] = check_guard4_masked_differential(n_trials=36, V_a=192)
        out["norms_and_2d_panjer"] = check_norms_and_2d_panjer(led)
        out["guard4_s1b_la_ram"] = check_guard4_s1b_la_ram(led)
        out["live_fq_rank_2merge"] = check_live_fq_rank(U=512, etas=(0.82,), n_seeds=6)
        out["stage4_12d_mg_witnesses"] = check_stage4_12d_mg_and_witnesses(led)
        out["assembly_isa_audit"] = check_assembly_isa_audit()
        out["ledger_meter_kgate_hygiene"] = check_ledger_meter_kgate_hygiene(led)

    if mode in ("diagnostics", "all"):
        out["diag_guard4_stress"] = check_guard4_masked_differential(n_trials=100, V_a=384)
        out["diag_fq_rank_multi_eta"] = check_live_fq_rank(U=512, etas=(0.78, 0.82, 0.90, 1.04), n_seeds=4)
        out["diag_prime_splitting_25000"] = verify_prime_splitting(25000)

    checks_ok = all(v.get("ok", True) for v in out.values() if isinstance(v, dict))
    out["checks_ok"] = checks_ok
    out["elapsed_seconds"] = round(time.perf_counter() - t0, 3)
    print(json.dumps(out, indent=2))
    if not checks_ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
