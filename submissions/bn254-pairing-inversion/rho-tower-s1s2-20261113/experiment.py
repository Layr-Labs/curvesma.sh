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

N = 3 * (math.isqrt(r) + 1) // 2
loop = 3 * (30000 + 20000 + 4096) + 1000
R = N * loop + (1 << 40)
getcontext().prec = 90
score = (Decimal(R)/20000).ln()/Decimal(2).ln()
score_up = str(score.quantize(Decimal("0.000001"), rounding=ROUND_CEILING))
print(json.dumps({
    "seed": seed, "montgomery_cases_passed": word_checks,
    "toy_walks": walks, "full_bn254_generator_order_check": True,
    "full_bn254_walk_invariant_check": True,
    "floyd_cap": str(N), "instructions_per_floyd_iteration_bound": loop,
    "total_instruction_bound": str(R), "score_rounded_up": score_up,
    "peak_bytes_bound": 1048576,
    "scope": "Small prime additive-group walk experiments, full-width Montgomery arithmetic and BN254 invariant checks. No production DLP, success-frequency proof, instruction measurement, or verifier acceptance is claimed."
}))
