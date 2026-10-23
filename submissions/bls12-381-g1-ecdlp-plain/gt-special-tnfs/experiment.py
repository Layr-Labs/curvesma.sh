"""BLS12-381 constant checks, exact ledger arithmetic, and a toy index calculus.

This is bounded evidence at toy scale. It does not execute or validate a
production number field sieve, pairing, or the submitted cost bound.
"""
import json
import random
import sys
from decimal import Decimal, ROUND_CEILING, getcontext
from fractions import Fraction

P = int("1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab", 16)
R = int("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001", 16)
U = -0xD201000000010000
P0 = (
    int("17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb", 16),
    int("08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1", 16),
)
Q0X = (
    int("024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8", 16),
    int("13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e", 16),
)
Q0Y = (
    int("0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801", 16),
    int("0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be", 16),
)


def fp_inv(a):
    return pow(a % P, P - 2, P)


def g1_add(a, b):
    if a is None:
        return b
    if b is None:
        return a
    x1, y1 = a
    x2, y2 = b
    if x1 == x2:
        if (y1 + y2) % P == 0:
            return None
        lam = 3 * x1 * x1 * fp_inv(2 * y1) % P
    else:
        lam = (y2 - y1) * fp_inv(x2 - x1) % P
    x3 = (lam * lam - x1 - x2) % P
    return x3, (lam * (x1 - x3) - y1) % P


def g1_mul(k, a):
    out = None
    for bit in bin(k)[2:]:
        out = g1_add(out, out)
        if bit == "1":
            out = g1_add(out, a)
    return out


def f2_add(a, b):
    return (a[0] + b[0]) % P, (a[1] + b[1]) % P


def f2_sub(a, b):
    return (a[0] - b[0]) % P, (a[1] - b[1]) % P


def f2_mul(a, b):
    return (a[0] * b[0] - a[1] * b[1]) % P, (a[0] * b[1] + a[1] * b[0]) % P


def f2_inv(a):
    n = (a[0] * a[0] + a[1] * a[1]) % P
    ni = fp_inv(n)
    return a[0] * ni % P, (-a[1]) * ni % P


def g2_add(a, b):
    if a is None:
        return b
    if b is None:
        return a
    x1, y1 = a
    x2, y2 = b
    if x1 == x2:
        if f2_add(y1, y2) == (0, 0):
            return None
        three = (3, 0)
        lam = f2_mul(f2_mul(three, f2_mul(x1, x1)), f2_inv(f2_add(y1, y1)))
    else:
        lam = f2_mul(f2_sub(y2, y1), f2_inv(f2_sub(x2, x1)))
    x3 = f2_sub(f2_sub(f2_mul(lam, lam), x1), x2)
    return x3, f2_sub(f2_mul(lam, f2_sub(x1, x3)), y1)


def g2_mul(k, a):
    out = None
    for bit in bin(k)[2:]:
        out = g2_add(out, out)
        if bit == "1":
            out = g2_add(out, a)
    return out


def check_constants():
    p_formula = ((U - 1) ** 2 * (U**4 - U**2 + 1) + 3 * U) // 3
    r_formula = U**4 - U**2 + 1
    checks = {
        "p_matches_profile": p_formula == P,
        "r_matches_profile": r_formula == R,
        "p_mod_3_is_1": P % 3 == 1,
        "r_divides_p4_minus_p2_plus_1": (P**4 - P**2 + 1) % R == 0,
        "P0_on_curve": (P0[1] ** 2 - P0[0] ** 3 - 4) % P == 0,
        "rP0_is_infinity": g1_mul(R, P0) is None,
        "Q0_on_curve": f2_sub(f2_mul(Q0Y, Q0Y), f2_add(f2_mul(f2_mul(Q0X, Q0X), Q0X), (4, 4))) == (0, 0),
        "rQ0_is_infinity": g2_mul(R, (Q0X, Q0Y)) is None,
    }
    if not all(checks.values()):
        raise AssertionError(checks)
    return checks


def check_ledger():
    elementary_per_log = 2**126
    relation_and_linear_algebra = elementary_per_log
    two_descents = elementary_per_log
    polynomial_selection = 2**20
    e = relation_and_linear_algebra + two_descents + polynomial_selection
    instructions_per_elementary = 64
    per_execution_instructions = instructions_per_elementary * e
    executions = 4
    amplified = executions * per_execution_instructions
    pairing_instructions = 2**42
    representation_instructions = 2**36
    io_instructions = 2**20
    other_instructions = pairing_instructions + representation_instructions + io_instructions
    r_total = amplified + other_instructions
    h_reference_operations = 2**10
    c_reference = 2000
    w = Fraction(h_reference_operations * c_reference + r_total, c_reference)
    getcontext().prec = 140
    log2_w = (Decimal(w.numerator).ln() - Decimal(w.denominator).ln()) / Decimal(2).ln()
    scaled = log2_w * Decimal(10**6)
    ceiling_scaled = int(scaled.to_integral_value(rounding=ROUND_CEILING))
    submitted = Fraction(ceiling_scaled, 10**6)
    if submitted != Fraction(15504277, 125000):
        raise AssertionError(submitted)
    expected = {
        "elementary": 170141183460469231731687303715885154304,
        "per_execution_instructions": 10889035741470030830827987437816649875456,
        "amplified_instructions": 43556142965880123323311949751266599501824,
        "other_instructions": 4466767036416,
        "r_total": 43556142965880123323311949755733366538240,
        "w_numerator": 544451787073501541541399371946667107328,
        "w_denominator": 25,
    }
    actual = {
        "elementary": e,
        "per_execution_instructions": per_execution_instructions,
        "amplified_instructions": amplified,
        "other_instructions": other_instructions,
        "r_total": r_total,
        "w_numerator": w.numerator,
        "w_denominator": w.denominator,
    }
    if actual != expected:
        raise AssertionError(actual)
    preprocessing = Fraction(2**126 * 64, 2000)
    prep_log2 = (Decimal(preprocessing.numerator).ln() - Decimal(preprocessing.denominator).ln()) / Decimal(2).ln()
    prep_ceiling = int((prep_log2 * Decimal(10**6)).to_integral_value(rounding=ROUND_CEILING))
    return {
        "elementary_per_log": str(elementary_per_log),
        "elementary_total": str(e),
        "instructions_per_elementary": instructions_per_elementary,
        "per_execution_instructions": str(per_execution_instructions),
        "executions": executions,
        "amplified_instructions": str(amplified),
        "other_instructions": str(other_instructions),
        "r_total": str(r_total),
        "h_reference_operations": h_reference_operations,
        "c_reference": c_reference,
        "w": f"{w.numerator}/{w.denominator}",
        "log2_w": str(log2_w),
        "submitted_time_log2": str(submitted),
        "submitted_fraction": "15504277/125000",
        "preprocessing_log2_ceiling": f"{prep_ceiling / 10**6:.6f}",
    }


def is_probable_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d = n - 1
    s = 0
    while d % 2 == 0:
        s += 1
        d //= 2
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def safe_prime(bits, rng):
    while True:
        q = rng.getrandbits(bits - 1) | (1 << (bits - 2)) | 1
        p = 2 * q + 1
        if p.bit_length() != bits:
            continue
        if is_probable_prime(q) and is_probable_prime(p):
            return p, q


def primes_upto(n):
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = b"\x00" * len(sieve[i * i :: i])
    return [i for i, v in enumerate(sieve) if v]


def toy_index_calculus(bits, seed, trials=2):
    rng = random.Random(seed)
    p, q = safe_prime(bits, rng)
    while True:
        a = rng.randrange(2, p - 1)
        g = a * a % p
        if g != 1:
            break
    bound = {20: 100, 24: 200, 28: 400, 32: 1000}.get(bits, max(50, 2 * bits))
    factor_base = primes_upto(bound)
    nfb = len(factor_base)

    def factor_smooth(m):
        exps = []
        for prime in factor_base:
            e = 0
            while m % prime == 0:
                m //= prime
                e += 1
            exps.append(e)
        return exps, m

    rows = []
    rhs = []
    val = 1
    k = 0
    target_relations = 8 * nfb + 64
    while len(rows) < target_relations and k < 400000:
        k += 1
        val = val * g % p
        exps, rem = factor_smooth(val)
        if rem == 1:
            rows.append(exps)
            rhs.append(k)
    if len(rows) < nfb + 8:
        raise AssertionError("toy relation collection failed")

    matrix = [row[:] + [rhs[i]] for i, row in enumerate(rows)]
    m = len(matrix)
    rank = 0
    pivots = []
    for col in range(nfb):
        pivot = None
        for r in range(rank, m):
            if matrix[r][col] % q:
                pivot = r
                break
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inv_pivot = pow(matrix[rank][col], q - 2, q)
        matrix[rank] = [(x * inv_pivot) % q for x in matrix[rank]]
        for r in range(m):
            if r != rank and matrix[r][col] % q:
                f = matrix[r][col]
                matrix[r] = [(matrix[r][j] - f * matrix[rank][j]) % q for j in range(nfb + 1)]
        pivots.append(col)
        rank += 1
        if rank == m:
            break
    logs = {col: matrix[i][nfb] for i, col in enumerate(pivots)}
    pivot_set = set(pivots)
    if rank != nfb:
        raise AssertionError(f"toy matrix rank {rank} below factor base {nfb}")

    recovered = []
    for _ in range(trials):
        secret = rng.randrange(1, q)
        h = pow(g, secret, p)
        found = None
        for t in range(20000):
            v = h * pow(g, t, p) % p
            exps, rem = factor_smooth(v)
            if rem != 1:
                continue
            if any(exps[col] and col not in pivot_set for col in range(nfb)):
                continue
            cand = (sum(exps[col] * logs[col] for col in pivot_set) - t) % q
            if pow(g, cand, p) == h:
                found = {"secret_matches": cand == secret, "offset": t}
                break
        if found is None or not found["secret_matches"]:
            raise AssertionError("toy individual logarithm failed")
        recovered.append(found)
    return {
        "bits": bits,
        "seed": seed,
        "group_order_bits": q.bit_length(),
        "factor_base_size": nfb,
        "relations": len(rows),
        "rank": rank,
        "trials": trials,
        "recovered": recovered,
    }


def main():
    payload = json.load(sys.stdin)
    bits = int(payload["bits"])
    seed = int(payload.get("seed", 1))
    if not 20 <= bits <= 32:
        raise AssertionError("bits out of range")
    constants = check_constants()
    ledger = check_ledger()
    toy = toy_index_calculus(bits, seed)
    print(json.dumps({
        "constants_ok": all(constants.values()),
        "ledger_ok": True,
        "ledger": ledger,
        "toy": toy,
        "scope": "Toy index calculus and exact ledger arithmetic only; no production NFS, pairing, or cost validation.",
    }, indent=2))


if __name__ == "__main__":
    main()
