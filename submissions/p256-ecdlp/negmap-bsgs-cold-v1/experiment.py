"""Finite-scope checks for the P-256 negation-map BSGS package (proof.md).

Evidence only at the stated scope. Nothing here executes production P-256
attack work, allocates the exponential table, times anything, or proves the
production cost bound; that bound is analytic (proof.md).

Modes:
  additive - signed-window coverage/sign logic on Z/nZ (class key min(v,n-v)).
  curve    - same logic on small prime-order curves y^2=x^3+ax+b with real
             x-coordinate keys: a=-3 curves (the P-256 shape) and generic a.
  kernel   - the exact 4x64-bit CIOS Montgomery multiplication of proof.md over
             the real P-256 field, checked against Python big integers, with the
             proof's per-instruction tally; plus the Fermat inversion chain
             length for the real exponent p-2.
  ledger   - exact integer checks of coverage vs the real P-256 order r and of
             the work ledger against 1000*2^145.
"""
import json
import sys

P256_P = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
P256_R = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
P256_A = (P256_P - 3)
P256_B = 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b
P256_GX = 0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296
P256_GY = 0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5
MASK = (1 << 64) - 1


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


def curve(p, a, b, m):
    if (4 * a ** 3 + 27 * b * b) % p == 0:
        return None
    add, mul = make_curve(p, a, b)
    pts = [(x, y) for x in range(p) for y in range(p)
           if (y * y - x ** 3 - a * x - b) % p == 0]
    n = len(pts) + 1
    if not is_prime(n) or 2 * m - 1 > n:
        return None
    G = pts[0]
    # Repeated generic additions must reach infinity exactly at n.
    assert mul(n, G) is None
    key = lambda R: ("inf",) if R is None else ("x", R[0])
    s, c, h = plan(n, m)
    table = {}
    R = None
    for j in range(m):
        k = key(R)
        assert k not in table, "two distinct classes shared an x-coordinate"
        table[k] = (j, R)   # stored record: index and the point (its y decides sign)
        R = add(R, G)
    D = mul((n - s) % n, G)           # D = -[s]P
    covered = h * s
    recovered = 0
    for x in range(n):
        Q = mul(x, G)
        gamma = add(Q, mul((n - c) % n, G))   # Q - [c]P
        out = None
        for i in range(h):
            k = key(gamma)
            if k in table:
                j, Rj = table[k]
                ysign_plus = (Rj == gamma)   # same x: equal y means +jP
                a_ = (i * s + c + j) % n if ysign_plus else (i * s + c - j) % n
                if mul(a_, G) == Q:
                    out = a_
                break
            gamma = add(gamma, D)
        if x < covered:
            assert out == x, ("curve: covered scalar not recovered", p, a, b, x)
            recovered += 1
        else:
            assert out is None or mul(out, G) == Q
    return {"field": p, "a": a, "b": b, "group_order": n, "m": m, "stride": s,
            "h": h, "covered": covered, "recovered": recovered,
            "exact_coverage": recovered == covered}


def curves():
    runs = []
    # a = -3 (the P-256 shape): first few prime-order curves per field.
    for p in (101, 103, 107, 109, 113, 127):
        found = 0
        for b in range(1, p):
            out = curve(p, (p - 3) % p, b, 4)
            if out:
                out["family"] = "a=-3"
                runs.append(out)
                found += 1
                if found == 2:
                    break
    # generic short Weierstrass a not in {0,-3}.
    for p in (101, 107, 113):
        found = 0
        for a in (2, 5, 7):
            for b in range(1, p):
                out = curve(p, a, b, 4)
                if out:
                    out["family"] = "generic"
                    runs.append(out)
                    found += 1
                    break
    return runs


# ---- 4x64-bit CIOS Montgomery multiplication with instruction tally ----
# Tally follows proof.md: every 64x64 product = MUL+UMULH (2), every
# accumulate = ADDS/ADCS pair or ADC (counted individually), every limb load
# (LDR) and store (STR) counted, plus final conditional subtraction (SUBS/SBCS
# chain + CSEL) and a fixed call/prologue/epilogue allowance.

class Tally:
    def __init__(self):
        self.n = 0

    def __call__(self, k=1):
        self.n += k


def mont_mul(a, b, p, n0, t_):
    A = [(a >> (64 * i)) & MASK for i in range(4)]
    B = [(b >> (64 * i)) & MASK for i in range(4)]
    P = [(p >> (64 * i)) & MASK for i in range(4)]
    t_(8)  # LDR a[0..3], b[0..3] (p limbs are immediates/regs: 4 MOV below)
    t_(4)  # materialize p limbs
    T = [0] * 6
    for i in range(4):
        C = 0
        for j in range(4):
            prod = A[j] * B[i]
            t_(2)                     # MUL, UMULH
            s = T[j] + (prod & MASK) + C
            t_(3)                     # ADDS lo, ADC hi, ADDS carry-in
            T[j] = s & MASK
            C = (prod >> 64) + (s >> 64)
            t_(1)                     # ADC carry-out
        s = T[4] + C
        T[4] = s & MASK
        T[5] = s >> 64
        t_(2)                         # ADDS, ADC
        mq = (T[0] * n0) & MASK
        t_(1)                         # MUL
        C = 0
        for j in range(4):
            prod = mq * P[j]
            t_(2)
            s = T[j] + (prod & MASK) + C
            t_(3)
            T[j] = s & MASK
            C = (prod >> 64) + (s >> 64)
            t_(1)
        s = T[4] + C
        T[4] = s & MASK
        T[5] += s >> 64
        t_(2)
        T = T[1:6] + [0]              # limb shift = register renaming
        t_(5)                         # up to 5 MOVs if not unrolled-renamed
    val = sum(T[i] << (64 * i) for i in range(5))
    t_(5)                             # SUBS/SBCS x4 + SBCS top
    if val >= p:
        val -= p
    t_(4)                             # CSEL x4
    t_(4)                             # STR result
    t_(24)                            # call/ret, prologue/epilogue saves
    return val


def kernel(trials=300):
    p = P256_P
    n0 = (-pow(p, -1, 1 << 64)) % (1 << 64)
    Rm = 1 << 256
    Rinv = pow(Rm, -1, p)
    seed = 0x5eed
    worst = 0
    ok = True
    for k in range(trials):
        seed = (seed * 6364136223846793005 + 1442695040888963407) % (1 << 64)
        a = pow(seed, 7, p)
        b = pow(seed + 1, 11, p)
        if k == 0:
            a, b = p - 1, p - 1        # edge case
        t = Tally()
        got = mont_mul(a, b, p, n0, t)
        ok &= (got == a * b * Rinv % p) and got < p
        worst = max(worst, t.n)
    # Fermat inversion chain for exponent p-2 (left-to-right binary).
    e = p - 2
    squarings = e.bit_length() - 1
    multiplies = bin(e).count("1") - 1
    # check the chain computes the inverse on Montgomery representatives
    t = Tally()
    x = (P256_GX * Rm) % p
    acc = x
    for bit in bin(e)[3:]:
        acc = mont_mul(acc, acc, p, n0, t)
        if bit == "1":
            acc = mont_mul(acc, x, p, n0, t)
    inv_ok = (acc * Rinv % p) * P256_GX % p == 1
    # P-256 generator on curve (a=-3) and canonical
    on_curve = (P256_GY ** 2 - P256_GX ** 3 - P256_A * P256_GX - P256_B) % p == 0
    return {"field": "P-256", "n0_prime": n0, "trials": trials,
            "mont_mul_correct": bool(ok),
            "mont_mul_worst_tally_instructions": worst,
            "mont_mul_budget_charged": 8192,
            "fermat_squarings": squarings, "fermat_multiplies": multiplies,
            "fermat_total_mults": squarings + multiplies,
            "fermat_budget_mults": 512, "fermat_chain_inverts": bool(inv_ok),
            "generator_on_curve": bool(on_curve)}


def ledger():
    r = P256_R
    m = 1 << 127
    s, c, h = 2 * m - 1, m - 1, (1 << 127) + 1
    cov = h * s
    R_pre = (1 << 148) + m * (1 << 27) + (1 << 48)
    R_total = R_pre + h * (1 << 26)
    step = 4210688 + 19 * 8192 + 32 * 256 + 4096
    return {"coverage": str(cov), "coverage_below_r": cov < r,
            "twice_coverage_above_r": 2 * cov > r,
            "no_duplicate_keys_condition_m_minus_1_below_r_half": m - 1 < r // 2,
            "R_pre": str(R_pre), "R_total": str(R_total),
            "R_total_below_1000_2_145": R_total < 1000 * (1 << 145),
            "R_pre_below_1000_2_145": R_pre < 1000 * (1 << 145),
            "point_step_bound": step, "point_step_below_2_23": step < (1 << 23),
            "baby_record_below_2_27": (1 << 26) + (1 << 23) < (1 << 27),
            "giant_iter_below_2_26": 128 * (1 << 18) + (1 << 23) < (1 << 26)}


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
        res["curve_runs_a_minus_3"] = sum(r["family"] == "a=-3" for r in runs)
        res["curve_runs_generic"] = sum(r["family"] == "generic" for r in runs)
        res["curve_all_exact"] = all(x["exact_coverage"] for x in runs)
    if mode in ("all", "kernel"):
        res["kernel"] = kernel()
    if mode in ("all", "ledger"):
        res["ledger"] = ledger()
    print(json.dumps(res))


if __name__ == "__main__":
    main()
