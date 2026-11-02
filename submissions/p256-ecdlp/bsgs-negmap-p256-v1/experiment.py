"""Finite-scope checks for the p256-ecdlp plain negation-map BSGS package.

Evidence only at the stated scope. Nothing here executes production attack
work, allocates the exponential table, times anything, or proves the
production cost bound; that bound is analytic (proof.md).

Modes:
  additive - signed-window coverage/sign logic on Z/nZ (class key min(v,n-v)).
  curve    - the full algorithm on small curves y^2=x^3-3x+b (the P-256
             shape) restricted to prime-order subgroups (cofactor 1 like
             production, plus cofactor cases), real x-coordinate keys.
  kernel   - the exact 6x64-bit CIOS Montgomery multiplication of proof.md over
             the real P-256 base field, checked against Python integers,
             with the per-instruction tally; the Fermat chain for p-2; and the
             protected generator P0 on the curve with [n]P0 = infinity.
  ledger   - exact integer checks of full coverage vs the real order n and
             of the work ledger against 1000*2^145.5.
"""
import json
import sys

P256_P = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
P256_N = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551
P256_A = P256_P - 3
P256_B = 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b
P0X = 0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296
P0Y = 0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5
MASK = (1 << 64) - 1
LIMBS = 6


def plan(n, m):
    s = 2 * m - 1
    c = m - 1
    h = n // s
    return s, c, h


def additive(n, m):
    s, c, h = plan(n, m)
    key = lambda v: min(v % n, (-v) % n)
    table = {}
    for j in range(m):
        k = key(j)
        assert k not in table
        table[k] = j
    covered = h * s
    recovered = 0
    for x in range(n):
        out = None
        for i in range(h):
            t = (x - i * s - c) % n
            k = key(t)
            if k in table:
                j = table[k]
                for a in ((i * s + c + j) % n, (i * s + c - j) % n):
                    if a == x:
                        out = a
                break
        if x < covered:
            assert out == x, ("covered scalar not recovered", n, m, x)
            recovered += 1
        else:
            assert out is None, ("uncovered scalar returned a value", n, m, x)
    return {"order": n, "m": m, "stride": s, "h": h, "covered": covered,
            "recovered": recovered, "exact_coverage": recovered == covered}


def make_curve(p, a, b):
    def add(P, Q):
        # Complete affine addition: infinity, opposite points and doubling
        # resolved before the single inversion, as in proof.md.
        if P is None:
            return Q
        if Q is None:
            return P
        if P[0] == Q[0] and (P[1] + Q[1]) % p == 0:
            return None
        if P == Q:
            num = (3 * P[0] * P[0] + a) % p
            den = (2 * P[1]) % p
        else:
            num = (Q[1] - P[1]) % p
            den = (Q[0] - P[0]) % p
        l = num * pow(den, p - 2, p) % p
        x3 = (l * l - P[0] - Q[0]) % p
        return (x3, (l * (P[0] - x3) - P[1]) % p)

    def mul(k, P):
        R = None
        while k:
            if k & 1:
                R = add(R, P)
            P = add(P, P)
            k >>= 1
        return R
    return add, mul


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


def largest_prime_factor(N):
    f, d = 1, 2
    while d * d <= N:
        while N % d == 0:
            f, N = d, N // d
        d += 1
    return max(f, N) if N > 1 else f


def curve(p, a, b, m, need_cofactor):
    if (4 * a ** 3 + 27 * b * b) % p == 0:
        return None
    add, mul = make_curve(p, a, b)
    pts = [(x, y) for x in range(p) for y in range(p)
           if (y * y - x ** 3 - a * x - b) % p == 0]
    N = len(pts) + 1
    n = largest_prime_factor(N)
    cof = N // n
    if n < 2 * m - 1 or n % 2 == 0 or cof % n == 0 or (need_cofactor and cof == 1):
        return None
    G = None
    for Pt in pts:  # generator of the order-n subgroup: cofactor-clear a point
        R = mul(cof, Pt)
        if R is not None:
            G = R
            break
    assert G is not None and mul(n, G) is None
    key = lambda R: ("inf",) if R is None else ("x", R[0])
    s, c, h = plan(n, m)
    table = {}
    R = None
    for j in range(m):
        k = key(R)
        assert k not in table, "two distinct classes shared an x-coordinate"
        table[k] = (j, R)
        R = add(R, G)
    D = mul((n - s) % n, G)
    covered = h * s
    recovered = 0
    for x in range(n):
        Q = mul(x, G)
        assert mul(n, Q) is None          # subgroup check [n]Q = O
        gamma = add(Q, mul((n - c) % n, G))
        out = None
        for i in range(h):
            k = key(gamma)
            if k in table:
                j, Rj = table[k]
                a_ = (i * s + c + j) % n if Rj == gamma else (i * s + c - j) % n
                if mul(a_, G) == Q:
                    out = a_
                break
            gamma = add(gamma, D)
        if x < covered:
            assert out == x, ("curve: covered scalar not recovered", p, a, b, x)
            recovered += 1
        else:
            assert out is None or mul(out, G) == Q
    return {"field": p, "a": a, "b": b, "curve_order": N, "cofactor": cof,
            "group_order": n, "m": m, "stride": s, "h": h, "covered": covered,
            "recovered": recovered, "exact_coverage": recovered == covered}


def curves():
    runs = []
    # P-256 shape y^2 = x^3 - 3x + b, cofactor 1 (like production).
    for p in (103, 109, 127, 139, 151, 157, 163):
        found = 0
        for b in range(1, p):
            out = curve(p, p - 3, b, 4, False)
            if out and out["cofactor"] == 1:
                out["family"] = "a=-3 prime"
                runs.append(out)
                found += 1
                if found == 2:
                    break
    # P-256 shape with cofactor (robustness).
    for p in (101, 107, 113):
        for b in range(1, p):
            out = curve(p, p - 3, b, 4, True)
            if out:
                out["family"] = "a=-3 cofactor"
                runs.append(out)
                break
    return runs


class Tally:
    def __init__(self):
        self.n = 0

    def __call__(self, k=1):
        self.n += k


def mont_mul(a, b, p, n0, t_):
    """6x64-bit CIOS Montgomery multiplication with the proof.md tally."""
    L = LIMBS
    A = [(a >> (64 * i)) & MASK for i in range(L)]
    B = [(b >> (64 * i)) & MASK for i in range(L)]
    P = [(p >> (64 * i)) & MASK for i in range(L)]
    t_(2 * L)  # LDR a, b limbs
    t_(L)      # materialize p limbs
    T = [0] * (L + 2)
    for i in range(L):
        C = 0
        for j in range(L):
            prod = A[j] * B[i]
            t_(2)                     # MUL, UMULH
            s = T[j] + (prod & MASK) + C
            t_(3)                     # ADDS lo, ADC hi, ADDS carry-in
            T[j] = s & MASK
            C = (prod >> 64) + (s >> 64)
            t_(1)                     # ADC carry-out
        s = T[L] + C
        T[L] = s & MASK
        T[L + 1] = s >> 64
        t_(2)
        mq = (T[0] * n0) & MASK
        t_(1)                         # MUL
        C = 0
        for j in range(L):
            prod = mq * P[j]
            t_(2)
            s = T[j] + (prod & MASK) + C
            t_(3)
            T[j] = s & MASK
            C = (prod >> 64) + (s >> 64)
            t_(1)
        s = T[L] + C
        T[L] = s & MASK
        T[L + 1] += s >> 64
        t_(2)
        assert T[0] == 0
        T = T[1:] + [0]
        t_(L + 1)                     # limb shift moves
    val = sum(T[i] << (64 * i) for i in range(L + 1))
    t_(L + 1)                         # SUBS/SBCS chain incl. top
    if val >= p:
        val -= p
    t_(L)                             # CSEL
    t_(L)                             # STR
    t_(24)                            # call/ret, prologue/epilogue
    return val


def kernel(trials=300):
    p = P256_P
    n0 = (-pow(p, -1, 1 << 64)) % (1 << 64)
    Rm = 1 << (64 * LIMBS)
    Rinv = pow(Rm, -1, p)
    seed = 0x5eed
    worst = 0
    ok = True
    for k in range(trials):
        seed = (seed * 6364136223846793005 + 1442695040888963407) % (1 << 64)
        a = pow(seed, 7, p)
        b = pow(seed + 1, 11, p)
        if k == 0:
            a, b = p - 1, p - 1
        t = Tally()
        got = mont_mul(a, b, p, n0, t)
        ok &= (got == a * b * Rinv % p) and got < p
        worst = max(worst, t.n)
    e = p - 2
    t = Tally()
    x = (P0X * Rm) % p
    acc = x
    for bit in bin(e)[3:]:
        acc = mont_mul(acc, acc, p, n0, t)
        if bit == "1":
            acc = mont_mul(acc, x, p, n0, t)
    inv_ok = (acc * Rinv % p) * P0X % p == 1
    on_curve = (P0Y ** 2 - P0X ** 3 - P256_A * P0X - P256_B) % p == 0
    add, mul = make_curve(p, P256_A, P256_B)
    order_ok = mul(P256_N, (P0X, P0Y)) is None
    return {"field": "P-256 Fp", "limbs": LIMBS, "n0_prime": hex(n0),
            "trials": trials, "mont_mul_correct": bool(ok),
            "mont_mul_worst_tally_instructions": worst,
            "mont_mul_budget_charged": 8192,
            "fermat_squarings": e.bit_length() - 1,
            "fermat_multiplies": bin(e).count("1") - 1,
            "fermat_total_mults": e.bit_length() - 1 + bin(e).count("1") - 1,
            "fermat_budget_mults": 1024, "fermat_chain_inverts": bool(inv_ok),
            "generator_on_curve": bool(on_curve),
            "generator_order_r": bool(order_ok)}


def ledger():
    n = P256_N
    m = 1 << 127
    s, c, h = 2 * m - 1, m - 1, 1 << 128
    cov = h * s
    R_pre = (1 << 148) + m * (1 << 27) + (1 << 48)
    R_total = R_pre + h * (1 << 26)
    inv = 1024 * 8192 + (1 << 14)
    step = inv + 19 * 8192 + 32 * 256 + 4096
    rec = 127 * (1 << 19) + (1 << 18) + step
    it = 128 * (1 << 18) + step
    # W < 2^(291/2)  <=>  R^2 < 1000^2 * 2^291 (exact integers)
    bound_sq = (1000 ** 2) * (1 << 291)
    return {"coverage": str(cov), "coverage_above_n": cov > n - 1,
            "full_success_probability_one": True,
            "n_bits": n.bit_length(), "n_below_2_256": n < (1 << 256),
            "R_pre": str(R_pre), "R_total": str(R_total),
            "R_total_sq_below_1000sq_2_291": R_total * R_total < bound_sq,
            "R_pre_sq_below_1000sq_2_291": R_pre * R_pre < bound_sq,
            "inversion_bound": inv, "point_step_bound": step,
            "point_step_below_2_24": step < (1 << 24),
            "per_record": rec, "per_record_below_2_27": rec < (1 << 27),
            "per_iteration": it, "per_iteration_below_2_26": it < (1 << 26)}


def main():
    cfg = json.load(sys.stdin)
    mode = cfg.get("mode", "all")
    res = {"scope": "toy groups, one real-field kernel model and exact integer "
                    "checks; does not validate production cost bounds"}
    if mode in ("all", "additive"):
        runs = []
        for n in (101, 257, 521, 1009):
            for m in (2, 4, 8, 16):
                if 2 * m - 1 <= n:
                    runs.append(additive(n, m))
        res["additive_runs"] = len(runs)
        res["additive_all_exact"] = all(x["exact_coverage"] for x in runs)
    if mode in ("all", "curve"):
        runs = curves()
        res["curve"] = runs
        res["curve_runs_prime"] = sum(r["family"] == "a=-3 prime" for r in runs)
        res["curve_runs_cofactor"] = sum(r["family"] == "a=-3 cofactor" for r in runs)
        res["curve_all_exact"] = all(x["exact_coverage"] for x in runs)
    if mode in ("all", "kernel"):
        res["kernel"] = kernel()
    if mode in ("all", "ledger"):
        res["ledger"] = ledger()
    print(json.dumps(res))


if __name__ == "__main__":
    main()
