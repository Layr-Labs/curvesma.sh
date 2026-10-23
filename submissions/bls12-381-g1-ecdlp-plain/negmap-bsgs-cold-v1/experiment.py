"""Finite-scope checks for the BLS12-381 G1 plain negation-map BSGS package.

Evidence only at the stated scope. Nothing here executes production attack
work, allocates the exponential table, times anything, or proves the
production cost bound; that bound is analytic (proof.md).

Modes:
  additive - signed-window coverage/sign logic on Z/nZ (class key min(v,n-v)).
  curve    - the full algorithm on small curves y^2=x^3+b (j=0, the BLS G1
             shape) restricted to a prime-order subgroup WITH a nontrivial
             cofactor, real x-coordinate keys; plus generic-a curves.
  kernel   - the exact 6x64-bit CIOS Montgomery multiplication of proof.md over
             the real BLS12-381 base field, checked against Python integers,
             with the per-instruction tally; the Fermat chain for p-2; and the
             protected generator P0 on the curve with [r]P0 = infinity.
  ledger   - exact integer checks of coverage vs the real subgroup order r and
             of the work ledger against 2000*2^143.5.
"""
import json
import sys

BLS_P = int("1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f624"
            "1eabfffeb153ffffb9feffffffffaaab", 16)
BLS_R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
BLS_B = 4
P0X = int("17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac58"
          "6c55e83ff97a1aeffb3af00adb22c6bb", 16)
P0Y = int("08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3ed"
          "d03cc744a2888ae40caa232946c5e7e1", 16)
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
    # j = 0 (the BLS G1 shape y^2 = x^3 + b), p = 1 mod 3, nontrivial cofactor.
    for p in (103, 109, 127, 139, 151, 157, 163):
        found = 0
        for b in range(1, p):
            out = curve(p, 0, b, 4, True)
            if out:
                out["family"] = "j=0 cofactor"
                runs.append(out)
                found += 1
                if found == 2:
                    break
    # generic short Weierstrass a != 0.
    for p in (101, 107, 113):
        for a in (2, 5):
            for b in range(1, p):
                out = curve(p, a, b, 4, False)
                if out:
                    out["family"] = "generic"
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
    p = BLS_P
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
    on_curve = (P0Y ** 2 - P0X ** 3 - BLS_B) % p == 0
    add, mul = make_curve(p, 0, BLS_B)
    order_ok = mul(BLS_R, (P0X, P0Y)) is None
    return {"field": "BLS12-381 Fp", "limbs": LIMBS, "n0_prime": hex(n0),
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
    r = BLS_R
    m = 1 << 126
    s, c, h = 2 * m - 1, m - 1, 1 << 127
    cov = h * s
    R_pre = (1 << 147) + m * (1 << 27) + (1 << 48)
    R_total = R_pre + h * (1 << 26)
    inv = 1024 * 8192 + (1 << 14)
    step = inv + 19 * 8192 + 32 * 256 + 4096
    rec = 126 * (1 << 19) + (1 << 18) + step
    it = 128 * (1 << 18) + step
    # W < 2^(287/2)  <=>  R^2 < 2000^2 * 2^287 (exact integers)
    bound_sq = (2000 ** 2) * (1 << 287)
    return {"coverage": str(cov), "coverage_below_r": cov < r,
            "twice_coverage_above_r": 2 * cov > r,
            "m_minus_1_below_r_half": m - 1 < r // 2,
            "R_pre": str(R_pre), "R_total": str(R_total),
            "R_total_sq_below_2000sq_2_287": R_total * R_total < bound_sq,
            "R_pre_sq_below_2000sq_2_287": R_pre * R_pre < bound_sq,
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
        res["curve_runs_j0_cofactor"] = sum(r["family"] == "j=0 cofactor" for r in runs)
        res["curve_runs_generic"] = sum(r["family"] == "generic" for r in runs)
        res["curve_all_exact"] = all(x["exact_coverage"] for x in runs)
    if mode in ("all", "kernel"):
        res["kernel"] = kernel()
    if mode in ("all", "ledger"):
        res["ledger"] = ledger()
    print(json.dumps(res))


if __name__ == "__main__":
    main()
