import importlib.util
import json
import math
import random
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
seed = request["seed"]
trials = request["trials"]
coins = random.Random(seed)
profile = json.loads((base / "target-profile.json").read_text())
p = int(profile["field_modulus"], 16)
r = int(profile["subgroup_order"], 16)

word_checks = 0
for m in (p, r):
    Rinv = pow(1 << 256, -1, m)
    cases = [(a, b) for a in (0, 1, m-2, m-1) for b in (0, 1, m-2, m-1)]
    cases += [(coins.randrange(m), coins.randrange(m)) for _ in range(256)]
    for a, b in cases:
        assert ref.montgomery_word_reference(a, b, m) == (a * b * Rinv) % m
        word_checks += 1

# This independent stream preserves the original toy-walk inputs.
field_coins = random.Random(seed ^ 0x4B415241)
zero = (0,) * 12
one = (1,) + (0,) * 11
monomial = lambda i, v=1: tuple(v if j == i else 0 for j in range(12))
edges = [zero, one, (p-1,)*12, monomial(11), monomial(6, p-1)]
field_cases = [(a,b) for a in edges for b in edges]
field_cases += [(tuple(field_coins.randrange(p) for _ in range(12)),
                 tuple(field_coins.randrange(p) for _ in range(12)))
                for _ in range(96)]
field_checks = 0
for aa, bb in field_cases:
    operations = {"mul": 0, "linear": 0}
    actual = ref.fp12_mul(aa, bb, p, operations)
    assert actual == ref.fp12_mul_schoolbook(aa, bb, p)
    assert operations == {"mul": 81, "linear": 362}
    assert all(0 <= x < p for x in actual)
    field_checks += 1
# X^22 reduces to 242*X^10-1476*X^4, exercising cascading reduction.
assert ref.fp12_mul(monomial(11), monomial(11), p) == tuple(
    242 if i == 10 else (p-1476 if i == 4 else 0) for i in range(12))

walks = []
for q in (1009, 10007, 100003, 1000003):
    counts = {k: 0 for k in ("success", "identity", "degenerate", "cap", "verification_failed")}
    steps = 0
    for _ in range(trials):
        secret = coins.randrange(q)
        mul = lambda a, b, q=q: (a + b) % q
        power = lambda a, target=1, q=q: a * target % q
        # Full 384-byte state length, as on BN254; only small additive groups.
        encode = lambda y: y.to_bytes(384, "big")
        got, n, outcome = ref.rho(q, secret, 0, mul, power, encode, coins)
        counts[outcome] += 1
        steps += n
        if got is not None:
            assert got == secret
    walks.append({"prime_order": q, "trials": trials, "outcomes": counts,
                  "success_fraction": (counts["success"]+counts["identity"])/trials,
                  "mean_floyd_iterations": steps/trials})

g = tuple(int(v, 16) for v in profile["generator_g"]["coefficients"])
one = (1,) + (0,) * 11
assert ref.fp12_pow(g, r, p) == one and g != one
# Full-width coefficients/scalars in invariant checks, not a full-size DLP run.
x = coins.randrange(r)
z = ref.fp12_pow(g, x, p)
a, b, c, d = [coins.randrange(r) for _ in range(4)]
y = ref.fp12_mul(ref.fp12_pow(g, a, p), ref.fp12_pow(z, b, p), p)
m = ref.fp12_mul(ref.fp12_pow(g, c, p), ref.fp12_pow(z, d, p), p)
for _ in range(4):
    y = ref.fp12_mul(y, m, p)
    a, b = (a+c)%r, (b+d)%r
assert y == ref.fp12_pow(g, (a+b*x)%r, p)

N = 2 * (math.isqrt(r) + 1)
loop = 3 * ((1 << 19) + (1 << 16) + (1 << 12)) + (1 << 16)
R = N * loop + (1 << 40)
getcontext().prec = 90
score = (Decimal(R)/20000).ln()/Decimal(2).ln()
score_up = str(score.quantize(Decimal("0.000001"), rounding=ROUND_CEILING))
print(json.dumps({
    "seed": seed, "montgomery_cases_passed": word_checks,
    "karatsuba_full_width_differential_cases_passed": field_checks,
    "karatsuba_per_product_primitives": {"mul": 81, "linear": 362},
    "cascading_reduction_check": True,
    "toy_walks": walks, "full_bn254_generator_order_check": True,
    "full_bn254_walk_invariant_check": True,
    "floyd_cap": str(N), "instructions_per_floyd_iteration_bound": loop,
    "total_instruction_bound": str(R), "score_rounded_up": score_up,
    "peak_bytes_bound": 1048576,
    "scope": "Full-width Karatsuba versus schoolbook differential checks, field-primitive counts, small prime additive-group walk experiments, Montgomery arithmetic and BN254 invariant checks. No production DLP, production success proof, compiled-C execution, instruction measurement, or verifier acceptance is claimed."
}))
