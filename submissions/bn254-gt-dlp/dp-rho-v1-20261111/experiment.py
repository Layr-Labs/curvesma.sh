import importlib.util
import json
import math
import random
import re
import sys
from decimal import Decimal, ROUND_CEILING, getcontext
from pathlib import Path

base = Path("/snapshot/files")
if not base.is_dir():
    base = Path(__file__).parent
spec = importlib.util.spec_from_file_location("rho_ref", base / "rho.py")
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)
request = json.load(sys.stdin)
seed, trials = request["seed"], request["trials"]
coins = random.Random(seed)
profile = json.loads((base / "target-profile.json").read_text())
p = int(profile["field_modulus"], 16)
r = int(profile["subgroup_order"], 16)

# 1. Word-level Montgomery multiplication (the fp_mul algorithm) -------------
word_checks = 0
for m in (p, r):
    Rinv = pow(1 << 256, -1, m)
    cases = [(a, b) for a in (0, 1, m - 2, m - 1) for b in (0, 1, m - 2, m - 1)]
    cases += [(coins.randrange(m), coins.randrange(m)) for _ in range(256)]
    for a, b in cases:
        assert ref.mont_mul_words(a, b, m) == (a * b * Rinv) % m
        word_checks += 1
assert (-pow(p, -1, 1 << 64)) % (1 << 64) == 0x87D20782E4866389

# 2. Tower Karatsuba formulas against the protected polynomial basis ---------
tower_checks = 0
for _ in range(64):
    a = tuple(coins.randrange(p) for _ in range(12))
    b = tuple(coins.randrange(p) for _ in range(12))
    ta, tb = ref.poly_to_tower(a, p), ref.poly_to_tower(b, p)
    assert ref.tower_to_poly(ta, p) == a
    assert ref.tower_to_poly(ref.fp12_mul(ta, tb, p), p) == ref.poly_mul(a, b, p)
    tower_checks += 1

# 3. Protected generator, full-width walk invariant with the real partition --
g = ref.poly_to_tower(tuple(int(v, 16) for v in profile["generator_g"]["coefficients"]), p)
one = ref.tower_one()
assert ref.fp12_pow(g, r, p) == one and g != one
x = coins.randrange(r)
z = ref.fp12_pow(g, x, p)
salt = coins.getrandbits(64)
tab = []
for _ in range(4):
    c, d = coins.randrange(r), coins.randrange(r)
    tab.append((ref.fp12_mul(ref.fp12_pow(g, c, p), ref.fp12_pow(z, d, p), p), c, d))
a, b = coins.randrange(r), coins.randrange(r)
y = ref.fp12_mul(ref.fp12_pow(g, a, p), ref.fp12_pow(z, b, p), p)
indices = []
for _ in range(24):
    wa, wb, wc, wd = ref.bn_words(y, p)
    j = ref.mix(wa, wb, salt) >> 54
    indices.append(j)
    mj, c, d = tab[j & 3]
    y, a, b = ref.fp12_mul(y, mj, p), (a + c) % r, (b + d) % r
assert y == ref.fp12_pow(g, (a + b * x) % r, p)

# 4. Partition / distinguished-word statistics on genuine GT elements --------
buckets = [0] * 16
low_bits = [0] * 8
y = z
samples = 1536
for _ in range(samples):
    y = ref.fp12_mul(y, tab[0][0], p)
    wa, wb, wc, wd = ref.bn_words(y, p)
    buckets[ref.mix(wa, wb, salt) >> 60] += 1
    for k in range(8):
        low_bits[k] += (wc >> k) & 1
chi2 = sum((v - samples / 16) ** 2 / (samples / 16) for v in buckets)

# 5. Toy walks: same step rule, mix, 1024 table, DP rule, cap, recovery ------
def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d, s = d // 2, s + 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
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

walks = []
for q, dp_bits in ((10007, 1), (100003, 2), (1000003, 4), (10000019, 5), (100000007, 7)):
    assert is_prime(q)
    k = (1 << 60) // q
    while not is_prime(k * q + 1):
        k += 1
    P = k * q + 1                      # 60-bit prime field, subgroup of order q
    h = 2
    while pow(h, k, P) == 1:
        h += 1
    gen = pow(h, k, P)
    t1, t2, t3 = (coins.randrange(2, P) for _ in range(3))
    words = lambda Y, P=P, t1=t1, t2=t2, t3=t3: (Y, Y * t1 % P, Y * t2 % P, Y * t3 % P)
    counts = {key: 0 for key in ("success", "identity", "degenerate", "cap", "table_full", "verification_failed")}
    steps = 0
    for _ in range(trials):
        secret = coins.randrange(q)
        target = pow(gen, secret, P)
        mul = lambda u, v, P=P: u * v % P
        power = lambda e, basis=gen, P=P: pow(basis, e, P)
        got, n, outcome = ref.rho(q, target, 1, mul, power, words, dp_bits, coins, 1 << 20)
        counts[outcome] += 1
        steps += n
        if got is not None:
            assert got == secret
    cap = ref.cap_steps(q)
    walks.append({"prime_order": q, "field_bits": P.bit_length(), "dp_bits": dp_bits,
                  "step_cap": cap, "trials": trials, "outcomes": counts,
                  "success_fraction": (counts["success"] + counts["identity"]) / trials,
                  "heuristic_prediction_ignoring_dp_delay": 1 - math.exp(-(1 - 1 / 1024) * cap * cap / (2 * q)),
                  "mean_steps": steps / trials})

# 6. Executed-instruction recount from the supplied disassembly --------------
fn, cur = {}, None
for line in (base / "arithmetic.disasm.txt").read_text().splitlines():
    head = re.match(r"^[0-9a-f]+ <_?(\w+)>:", line)
    if head:
        cur = head.group(1)
        fn[cur] = []
        continue
    ins = re.match(r"^\s*[0-9a-f]+:\s+(\S+)\s*(.*)$", line)
    if ins and cur:
        fn[cur].append((ins.group(1), ins.group(2)))
BRANCHES = ("cbz", "cbnz", "tbz", "tbnz", "br", "blr")
memo = {}
def executed(name):
    if name in memo:
        return memo[name]
    total = 0
    for op, arg in fn[name]:
        total += 1
        callee = re.search(r"<_?(\w+)>", arg)
        if op == "bl" or (op == "b" and callee and callee.group(1) != name):
            total += executed(callee.group(1))      # call or tail call
        else:
            assert not (op == "b" or op.startswith("b.") or op in BRANCHES), (name, op)
            assert not re.search(r"\b[vqdsz]\d+\b", arg), (name, op, arg)   # no SIMD/FP registers
    memo[name] = total
    return total
counts = {name: {"static": len(fn[name]), "executed": executed(name)} for name in fn if name != "main"}
assert counts["walk_step"]["executed"] <= 24576

# 7. Score ------------------------------------------------------------------
N = ref.cap_steps(r)
STEP, FIXED = 24576, 1 << 70
R = N * STEP + FIXED
getcontext().prec = 90
score = (Decimal(R) / 20000).ln() / Decimal(2).ln()
score_up = str(score.quantize(Decimal("0.000001"), rounding=ROUND_CEILING))
print(json.dumps({
    "seed": seed, "montgomery_cases_passed": word_checks, "tower_vs_profile_basis_cases": tower_checks,
    "full_bn254_generator_order_check": True, "full_bn254_walk_invariant_check": True,
    "full_bn254_partition_indices_sample": indices,
    "gt_partition_top4bit_buckets": buckets, "gt_partition_chi2_15dof": round(chi2, 3),
    "gt_dp_word_low_bit_ones_of_%d" % samples: low_bits,
    "toy_walks": walks, "instruction_counts_from_disassembly": counts,
    "step_cap_N": str(N), "instructions_per_step_bound": STEP, "fixed_instruction_allowance": str(FIXED),
    "total_instruction_bound": str(R), "score_rounded_up": score_up,
    "heuristic_success_prediction_production": 1 - math.exp(-(1 - 1 / 1024) * 0.845),
    "peak_bytes_bound": str(1 << 38),
    "scope": "Toy subgroup walks (orders up to 1e8), full-width Montgomery/tower arithmetic and BN254 invariant checks, and a static recount of the supplied straight-line disassembly. No production DLP, no production success-frequency proof, and no verifier acceptance is claimed."
}))
