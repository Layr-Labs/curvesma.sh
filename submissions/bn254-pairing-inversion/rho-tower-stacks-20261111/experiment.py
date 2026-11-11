"""Supporting checks for rho-tower-stacks (standard library only, bounded runtime).

Scope: correctness of the tower arithmetic against the protected basis, exact
operation counts per Fp12 product, the word-level CIOS Montgomery product with
an explicit scalar-instruction tally, BN254 walk invariants and output, toy
runs of the exact walk/stack control flow on small prime-order groups, the
success-model integral and the score.  Nothing here solves production DLP or
proves a production success probability.
"""
import importlib.util
import json
import math
import random
import sys
import time
from decimal import Decimal, ROUND_CEILING, getcontext
from pathlib import Path

base = Path("/snapshot/files")
if not base.is_dir():
    base = Path(__file__).parent
spec = importlib.util.spec_from_file_location("rho_ref", base / "rho_ref.py")
R = importlib.util.module_from_spec(spec)
spec.loader.exec_module(R)

req = json.load(sys.stdin)
seed = int(req["seed"])
toy_trials = int(req["toy_trials"])
coins = random.Random(seed)
t_start = time.time()
profile = json.loads((base / "target-profile.json").read_text())
p = int(profile["field_modulus"], 16)
r = int(profile["subgroup_order"], 16)
assert p == R.P and r == R.R
out = {"seed": seed}

# ---------------------------------------------------------------- 1. tower == protected basis
pairs = 0
for _ in range(48):
    a = tuple(coins.randrange(p) for _ in range(12))
    b = tuple(coins.randrange(p) for _ in range(12))
    assert R.tower_to_poly(R.poly_to_tower(a)) == a
    assert R.tower_to_poly(R.f12mul(R.poly_to_tower(a), R.poly_to_tower(b))) == R.poly_mul(a, b)
    pairs += 1
edge = [tuple([p - 1] * 12), tuple([0] * 11 + [p - 1]), tuple([1] + [0] * 11)]
for a in edge:
    for b in edge:
        assert R.tower_to_poly(R.f12mul(R.poly_to_tower(a), R.poly_to_tower(b))) == R.poly_mul(a, b)
        pairs += 1
R.CNT.reset()
R.f12mul(R.poly_to_tower(a), R.poly_to_tower(b))
out["tower_vs_protected_basis_products_checked"] = pairs
out["fp12_product_fp_mults"] = R.CNT.mul
out["fp12_product_fp_addsubs"] = R.CNT.add
assert (R.CNT.mul, R.CNT.add) == (54, 280)

# ---------------------------------------------------------------- 2. word-level CIOS with instruction tally
MASK = (1 << 64) - 1


def cios(a, b, m):
    """4-limb CIOS Montgomery product; returns (a*b/2^256 mod m, scalar-instruction tally).

    Tally convention (AArch64): MUL and UMULH are one instruction each; each
    128-bit accumulate t+x*y+c is MUL,UMULH,ADDS,ADC,ADDS,ADC (6); a carry
    absorption is ADDS,ADC (2); the final conditional subtraction is 5 SUBS/SBCS
    + 4 CSEL; loads of a,b,m,n0 = 13; stores = 4; call/return and callee-saved
    register save/restore = 12.
    """
    n0 = (-pow(m, -1, 1 << 64)) & MASK
    A = [(a >> (64 * i)) & MASK for i in range(4)]
    B = [(b >> (64 * i)) & MASK for i in range(4)]
    M = [(m >> (64 * i)) & MASK for i in range(4)]
    t = [0] * 6
    ins = 13 + 4 + 12
    for i in range(4):
        C = 0
        for j in range(4):
            v = t[j] + A[j] * B[i] + C
            t[j], C = v & MASK, v >> 64
            ins += 6
        v = t[4] + C
        t[4], t[5] = v & MASK, v >> 64
        ins += 2
        q = (t[0] * n0) & MASK
        ins += 1
        v = t[0] + q * M[0]
        assert v & MASK == 0
        C = v >> 64
        ins += 4
        for j in range(1, 4):
            v = t[j] + q * M[j] + C
            t[j - 1], C = v & MASK, v >> 64
            ins += 6
        v = t[4] + C
        t[3], C = v & MASK, v >> 64
        t[4] = t[5] + C
        ins += 3
    u = sum(t[i] << (64 * i) for i in range(5))
    assert u < 2 * m
    ins += 9
    return (u - m if u >= m else u), ins


word_cases = 0
tally = 0
for m in (p, r):
    rinv = pow(1 << 256, -1, m)
    cases = [(x, y) for x in (0, 1, m - 2, m - 1) for y in (0, 1, m - 2, m - 1)]
    cases += [(coins.randrange(m), coins.randrange(m)) for _ in range(200)]
    for x, y in cases:
        v, ins = cios(x, y, m)
        assert v == x * y * rinv % m
        tally = max(tally, ins)
        word_cases += 1
out["cios_cases_passed"] = word_cases
out["cios_instruction_tally"] = tally
assert tally <= 384

# ---------------------------------------------------------------- 3. BN254 group checks and walk invariant
g = R.poly_to_tower(tuple(int(c, 16) for c in profile["generator_g"]["coefficients"]))
assert g != R.ONE12 and R.f12pow(g, r) == R.ONE12
x = coins.randrange(1, r)
z = R.f12pow(g, x)
a0, b0 = coins.randrange(r), coins.randrange(r)
y = R.f12mul(R.f12pow(g, a0), R.f12pow(z, b0))
tab = []
for _ in range(4):
    c, d = coins.randrange(r), coins.randrange(r)
    tab.append((R.f12mul(R.f12pow(g, c), R.f12pow(z, d)), c, d))
h1 = coins.getrandbits(64) | 1
for _ in range(8):
    j = ((R.tower_words(y)[0] * h1) & MASK) >> 62   # 4-adding toy index on the real field
    m_, c, d = tab[j]
    y = R.f12mul(y, m_)
    a0, b0 = (a0 + c) % r, (b0 + d) % r
assert y == R.f12pow(g, (a0 + b0 * x) % r)
out["bn254_generator_order_and_walk_invariant"] = True

# ---------------------------------------------------------------- 4. output Q=[x]Q0
Q0 = (tuple(int(v, 16) for v in profile["generator_Q0"]["x"]),
      tuple(int(v, 16) for v in profile["generator_Q0"]["y"]))
assert R.g2_on_curve(Q0) and R.g2_mul(Q0, r) is None
Q = R.g2_mul(Q0, x)
assert Q is not None and R.g2_on_curve(Q)
out["output_Q_on_twist_and_nonidentity"] = True

# ---------------------------------------------------------------- 5. toy walks (same control flow, small groups)


def is_prime(n):
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for q in small:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:
        v = pow(a, d, n)
        if v in (1, n - 1):
            continue
        for _ in range(s - 1):
            v = v * v % n
            if v == n - 1:
                break
        else:
            return False
    return True


S_RATIO = Decimal("1.15")


def cap_for(n):
    # smallest T with T^2 >= 1.15^2 * n * 1024/1023
    num, den = 529 * n * 1024, 400 * 1023
    T = math.isqrt(num // den)
    while T * T * den < num:
        T += 1
    return T


def run_toys(bits, trials, kind):
    q = coins.getrandbits(bits) | (1 << (bits - 1)) | 1
    while not is_prime(q):
        q += 2
    if kind == "additive":
        # Z_q, residues hashed through a 384-byte-equivalent fixed word (the residue itself).
        gp = lambda e: e % q
        mul = lambda u, v: (u + v) % q
        mk_t = lambda zz: (lambda e: zz * e % q)
        mk_z = lambda xx: xx
        chk = lambda zz: zz
    else:
        k = 1
        while not is_prime(2 * k * q + 1):
            k += 1
        Pm = 2 * k * q + 1
        while True:
            gg = pow(coins.randrange(2, Pm - 1), 2 * k, Pm)
            if gg != 1:
                break
        gp = lambda e: pow(gg, e, Pm)
        mul = lambda u, v: u * v % Pm
        mk_t = lambda zz: (lambda e: pow(zz, e, Pm))
        mk_z = lambda xx: pow(gg, xx, Pm)
    cap = cap_for(q)
    outcomes = {}
    steps = 0
    peak = 0
    for _ in range(trials):
        xx = coins.randrange(q)
        zz = mk_z(xx)
        got, n, o, pk = R.walk(q, gp, mk_t(zz), mul, lambda u: u & MASK, lambda u: u, cap, coins)
        if got is not None:
            assert got == xx
        outcomes[o] = outcomes.get(o, 0) + 1
        steps += n
        peak = max(peak, pk)
    succ = outcomes.get("success", 0)
    return {"kind": kind, "prime_order_bits": q.bit_length(), "prime_order": q, "step_cap": cap,
            "cap_over_sqrt_order": round(cap / math.sqrt(q), 4), "trials": trials,
            "outcomes": outcomes, "success_fraction": round(succ / trials, 4),
            "binomial_2se": round(2 * math.sqrt(max(succ, 1) * (trials - succ) / trials) / trials, 4),
            "mean_steps": round(steps / trials, 1), "max_stack_depth": peak}


toys = [run_toys(16, toy_trials, "additive"),
        run_toys(20, toy_trials, "additive"),
        run_toys(22, max(16, toy_trials // 4), "additive"),
        run_toys(16, max(16, toy_trials // 8), "multiplicative")]
tot_s = sum(t["outcomes"].get("success", 0) for t in toys)
tot_n = sum(t["trials"] for t in toys)
out["toy_walks"] = toys
out["toy_pooled_success"] = round(tot_s / tot_n, 4)

# ---------------------------------------------------------------- 6. success model integral


def model(s, K=256, N=3000, M=200):
    def F(c):
        if c >= 1:
            return 1.0
        tot = 0.0
        for i in range(M):
            v = c + (1 - c) * (i + 0.5) / M
            tot += (1 - c / v) ** K
        return 1 - tot * (1 - c) / M
    tot, h = 0.0, s / N
    for i in range(N):
        xr = (i + 0.5) * h
        tot += xr * math.exp(-xr * xr / 2) * F((s - xr) / xr) * h
    return tot


out["model_success_at_1.15"] = round(model(1.15), 4)
out["model_success_if_effective_order_x1.2"] = round(model(1.15 / math.sqrt(1.2)), 4)
out["upper_bound_1_minus_exp(-s^2/2)"] = round(1 - math.exp(-1.15 ** 2 / 2), 4)

# ---------------------------------------------------------------- 7. score
T = cap_for(r)
S = 49152
F = 1 << 40
Rt = T * S + F
getcontext().prec = 90
score = (Decimal(Rt) / 20000).ln() / Decimal(2).ln()
out["production_step_cap_T"] = str(T)
out["per_step_instruction_bound"] = S
out["setup_and_final_instruction_allowance"] = "2^40"
out["total_instruction_bound"] = str(Rt)
out["score_rounded_up"] = str(score.quantize(Decimal("0.000001"), rounding=ROUND_CEILING))
out["peak_bytes_bound"] = 1 << 26
out["seconds"] = round(time.time() - t_start, 2)
out["scope"] = ("Arithmetic identities, exact operation counts, CIOS tally, BN254 invariants, "
                "toy-scale walk statistics and analytic model only. No production DLP run, no "
                "production success proof, no hardware instruction measurement.")
print(json.dumps(out))
