"""Toy-scale run of the exact attack in proof.md on small j=0 curves.

Finds a prime p = 1 mod 3 and b with E: y^2 = x^3 + b of prime order r,
6d | r-1, then for uniform secrets s runs: Cheon p-1 step on A_d=[s^d]P0
as BSGS in H/mu6 (H = d-th powers in F_r^*, mu6 = <-1, lambda>), using
x-coordinate class keys min{x, beta*x, beta^2*x} and 3-window fixed-base
tables, then resolves s among 6d candidates against A_1. Every group
addition is counted. Evidence at toy scope only; not a production cost proof.
"""
import json
import math
import random
import sys
import time

inp = json.load(sys.stdin)
BITS = inp["bits"]
DLOG = inp["d_log2"]
TRIALS = inp["trials"]
SEED = inp["seed"]
P_TARGET = 0.39
rng = random.Random(SEED)
t0 = time.time()


def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    dd, s = n - 1, 0
    while dd % 2 == 0:
        dd //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, dd, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def factor(n):
    fs, q = {}, 2
    while q * q <= n:
        while n % q == 0:
            fs[q] = fs.get(q, 0) + 1
            n //= q
        q += 1
    if n > 1:
        fs[n] = fs.get(n, 0) + 1
    return fs


ADDS = [0]
EXC = [0]


def add(P, Q, p):
    ADDS[0] += 1
    if P is None:
        EXC[0] += 1
        return Q
    if Q is None:
        EXC[0] += 1
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        EXC[0] += 1
        if (y1 + y2) % p == 0:
            return None
        lam = 3 * x1 * x1 * pow(2 * y1, -1, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


def smul(k, P, p):
    R = None
    while k:
        if k & 1:
            R = add(R, P, p)
        P = add(P, P, p)
        k >>= 1
    return R


def find_params():
    d = 1 << DLOG
    while True:
        p = rng.randrange(1 << (BITS - 1), 1 << BITS)
        if p % 3 != 1 or not is_prime(p):
            continue
        t = None
        for v in range(1, math.isqrt(4 * p // 3) + 1):
            t2 = 4 * p - 3 * v * v
            tt = math.isqrt(t2)
            if tt * tt == t2:
                t = (tt, v)
                break
        tt, v = t
        traces = {tt, -tt, (tt + 3 * v) // 2, -(tt + 3 * v) // 2,
                  (tt - 3 * v) // 2, -(tt - 3 * v) // 2}
        for T in traces:
            r = p + 1 - T
            if r % (6 * d) != 1 % (6 * d) or not is_prime(r):
                continue
            for b in range(1, 200):
                for x in range(1, 200):
                    rhs = (x ** 3 + b) % p
                    if pow(rhs, (p - 1) // 2, p) == 1:
                        y = pow(rhs, (p + 1) // 4, p) if p % 4 == 3 else None
                        if y is None:
                            y = tonelli(rhs, p)
                        P0 = (x, y)
                        break
                else:
                    continue
                if smul(r, P0, p) is None:
                    return p, b, r, P0


def tonelli(a, p):
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, R = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1:
            tt = tt * tt % p
            i += 1
        bb = pow(c, 1 << (m - i - 1), p)
        m, c, t, R = i, bb * bb % p, t * bb * bb % p, R * bb % p
    return R


p, b, r, P0 = find_params()
d = 1 << DLOG
N = (r - 1) // d
Np = N // 6
fr = factor(r - 1)
zeta = next(g for g in range(2, r) if all(pow(g, (r - 1) // q, r) != 1 for q in fr))
xi = pow(zeta, d, r)
beta = next(pow(a, (p - 1) // 3, p) for a in range(2, p) if pow(a, (p - 1) // 3, p) != 1)
assert pow(-1 % r, N, r) == 1, "-1 must lie in H"
lam = next(l for l in (pow(zeta, (r - 1) // 3, r), pow(zeta, 2 * (r - 1) // 3, r))
           if smul(l, P0, p) == (beta * P0[0] % p, P0[1]))
assert pow(lam, N, r) == 1, "lambda must lie in H"


def dlog_zeta(h):
    """Ground truth only (Pohlig-Hellman at toy scale); not part of the attack."""
    n, X = r - 1, 0
    for q, e in fr.items():
        qe = q ** e
        gq, hq, x = pow(zeta, n // qe, r), pow(h, n // qe, r), 0
        gam = pow(gq, qe // q, r)
        M = math.isqrt(q) + 1
        tab, c = {}, 1
        for j in range(M):
            tab.setdefault(c, j)
            c = c * gam % r
        step = pow(gam, -M, r)
        for k in range(e):
            y = pow(pow(gq, -x, r) * hq % r, qe // q ** (k + 1), r)
            for i in range(M + 1):
                if y in tab:
                    x += (i * M + tab[y]) * q ** k
                    break
                y = y * step % r
        X += x * (n // qe) * pow(n // qe, -1, qe)
    return X % n


W = -(-r.bit_length() // 3)
MASK = (1 << W) - 1


def build_tables(B):
    tabs = []
    for k in range(3):
        base = smul(1 << (W * k), B, p)
        row, acc = [None], None
        for _ in range((1 << W) - 1):
            acc = add(acc, base, p)
            row.append(acc)
        tabs.append(row)
    return tabs


def fixed_base(tabs, c):
    return add(add(tabs[0][c & MASK], tabs[1][(c >> W) & MASK], p), tabs[2][c >> (2 * W)], p)


def key(Pt):
    x = Pt[0]
    bx = beta * x % p
    return min(x, bx, (-x - bx) % p)


need = -(-int(P_TARGET * Np * 10**6) // 10**6)
while need * 10**6 < P_TARGET * 10**6 * Np:
    need += 1
mb = math.isqrt(need)
if mb * mb < need:
    mb += 1
m = g = mb
assert m * g <= Np, "toy group too small for the requested coverage"

a0 = ADDS[0]
TP = build_tables(P0)
table_adds_P0 = ADDS[0] - a0
a0 = ADDS[0]
baby, bj = {}, 1
for j in range(m):
    kk = key(fixed_base(TP, bj))
    assert kk not in baby, "baby classes must be distinct in H/mu6"
    baby[kk] = j
    bj = bj * xi % r
baby_adds = ADDS[0] - a0
gamma = pow(xi, -m, r)
omega = pow(zeta, Np, r)

succ = agree = disagree = wrong = false_match = done = 0
giant_adds = table_adds_Ad = recover_adds = 0
for trial in range(TRIALS):
    if time.time() - t0 > 15:
        break
    s = rng.randrange(1, r)
    A1 = smul(s, P0, p)
    Ad = smul(pow(s, d, r), P0, p)
    a0 = ADDS[0]
    TA = build_tables(Ad)
    table_adds_Ad += ADDS[0] - a0
    a0 = ADDS[0]
    gi, found = 1, None
    for i in range(g):
        kk = key(fixed_base(TA, gi))
        if kk in baby:
            found = (baby[kk] + m * i) % Np
            break
        gi = gi * gamma % r
    giant_adds += ADDS[0] - a0
    a0 = ADDS[0]
    got = None
    if found is not None:
        cand = pow(zeta, found, r)
        for u in range(6 * d):
            if fixed_base(TP, cand) == A1:
                got = cand
                break
            cand = cand * omega % r
        if got is None:
            false_match += 1
    recover_adds += ADDS[0] - a0
    done += 1
    covered = dlog_zeta(s) % Np < m * g
    if covered == (got == s):
        agree += 1
    else:
        disagree += 1
    if got is not None:
        if got == s:
            succ += 1
        else:
            wrong += 1

# Exact predicted success = m*g/Np, since s -> e = log_zeta s is uniform on Z/(r-1) and N' | r-1.
pred = m * g / Np
print(json.dumps({
    "params": {"p": p, "b": b, "r": r, "d": d, "N_prime": Np, "window_bits": W,
               "m": m, "g": g},
    "trials_run": done,
    "recovered_correct": succ,
    "recovered_wrong": wrong,
    "false_matches": false_match,
    "trials_where_success_equals_coverage_rule": agree,
    "trials_disagreeing_with_coverage_rule": disagree,
    "empirical_success": succ / done if done else None,
    "predicted_success_mg_over_Nprime": pred,
    "binomial_sd": math.sqrt(pred * (1 - pred) / done) if done else None,
    "adds": {"tables_P0": table_adds_P0, "baby": baby_adds,
             "tables_Ad_per_trial": table_adds_Ad / max(done, 1),
             "giant_per_trial_worst": 2 * g,
             "giant_per_trial_mean": giant_adds / max(done, 1),
             "recovery_per_trial_mean": recover_adds / max(done, 1)},
    "baby_adds_per_element": baby_adds / m,
    "exceptional_add_cases_incl_input_generation": EXC[0],
    "seconds": round(time.time() - t0, 2),
    "scope": "toy curve only; checks the algorithm, coverage and success law, not production cost",
}))
