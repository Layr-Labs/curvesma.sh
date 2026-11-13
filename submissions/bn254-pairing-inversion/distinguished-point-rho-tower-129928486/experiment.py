# Structure mirrors the track baseline's experiment.py (Montgomery word
# checks, toy walks on small groups, exact score recomputation from the
# protected order); the tower checks and distinguished-point toy walks
# are new to this package.
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

# 1. Full-width Montgomery word checks for both moduli (baseline's cases).
word_checks = 0
for m in (p, r):
    Rinv = pow(1 << 256, -1, m)
    cases = [(a, b) for a in (0, 1, m-2, m-1) for b in (0, 1, m-2, m-1)]
    cases += [(coins.randrange(m), coins.randrange(m)) for _ in range(256)]
    for a, b in cases:
        assert ref.montgomery_word_reference(a, b, m) == (a * b * Rinv) % m
        word_checks += 1

# 2. Tower embedding identities in the protected polynomial basis.
w = (0, 1) + (0,) * 10
one = (1,) + (0,) * 11
u = tuple((ref.fp12_pow(w, 6, p)[i] - (9 if i == 0 else 0)) % p
          for i in range(12))
neg_one = tuple((-1 if i == 0 else 0) % p for i in range(12))
assert ref.fp12_mul(u, u, p) == neg_one
v = ref.fp12_pow(w, 2, p)
nine_plus_u = tuple((u[i] + (9 if i == 0 else 0)) % p for i in range(12))
assert ref.fp12_mul(ref.fp12_mul(v, v, p), v, p) == nine_plus_u
assert ref.fp12_mul(w, w, p) == v

# 3. Tower/polynomial roundtrip + multiply equality on random full-width pairs,
#    and the exact instrumented Karatsuba operation count for one multiply.
tower_trials = 200
for _ in range(tower_trials):
    a = tuple(coins.randrange(p) for _ in range(12))
    b = tuple(coins.randrange(p) for _ in range(12))
    assert ref.to_poly(ref.to_tower(a, p), p) == a
    assert ref.to_poly(ref.to_tower(b, p), p) == b
    assert ref.to_poly(ref.fp12_tower_mul(ref.to_tower(a, p),
                                          ref.to_tower(b, p), p), p) \
        == ref.fp12_mul(a, b, p)
count = {"muls": 0, "adds": 0}
ta = ref.to_tower(tuple(coins.randrange(p) for _ in range(12)), p)
tb = ref.to_tower(tuple(coins.randrange(p) for _ in range(12)), p)
ref.fp12_tower_mul(ta, tb, p, count)
assert count == {"muls": 54, "adds": 280}

# 4. Protected generator order + walk invariant, in tower arithmetic.
g = tuple(int(vv, 16) for vv in profile["generator_g"]["coefficients"])
g_tower = ref.to_tower(g, p)
assert ref.to_poly(g_tower, p) == g
assert ref.fp12_tower_pow(g_tower, r, p) == ref.fp12_tower_pow(g_tower, 0, p)
assert g_tower != ref.fp12_tower_pow(g_tower, 0, p)
assert ref.fp12_pow(g, r, p) == one and g != one
x = coins.randrange(r)
z = ref.fp12_tower_pow(g_tower, x, p)
a, b, c, d = [coins.randrange(r) for _ in range(4)]
y = ref.fp12_tower_mul(ref.fp12_tower_pow(g_tower, a, p),
                       ref.fp12_tower_pow(z, b, p), p)
m = ref.fp12_tower_mul(ref.fp12_tower_pow(g_tower, c, p),
                       ref.fp12_tower_pow(z, d, p), p)
for _ in range(4):
    y = ref.fp12_tower_mul(y, m, p)
    a, b = (a + c) % r, (b + d) % r
assert y == ref.fp12_tower_pow(g_tower, (a + b * x) % r, p)

# 5. Distinguished-point toy walks on small prime additive groups, using the
#    actual rho() with the production step structure and a small theta so
#    distinguished points occur at toy scale. All outcomes are reported,
#    including failures. Toy scope only; no production claim is drawn.
walks = []
for q in (1009, 10007, 100003, 1000003):
    counts = {k: 0 for k in ("success", "identity", "degenerate", "cap",
                             "verification_failed", "table_full")}
    steps = 0
    for _ in range(trials):
        secret = coins.randrange(q)
        mul = lambda a, b, q=q: (a + b) % q
        power = lambda a, target=1, q=q: a * target % q
        # Full 384-byte state length, as on BN254; only small additive groups.
        encode = lambda y: y.to_bytes(384, "big")
        got, n, outcome = ref.rho(q, secret, 0, mul, power, encode, coins,
                                  dp_bits=4, record_cap=1 << 20)
        counts[outcome] += 1
        steps += n
        if got is not None:
            assert got == secret
    walks.append({"prime_order": q, "trials": trials, "outcomes": counts,
                  "success_fraction": (counts["success"] + counts["identity"]) / trials,
                  "mean_steps": steps / trials,
                  "theta": "2^-4",
                  "step_cap": 2 * (math.isqrt(q) + (math.isqrt(q) ** 2 < q))})

# 6. Exact capped-cost recomputation from the protected order.
K = 2 * (math.isqrt(r) + 1)
assert K == 295893513763578638011461384341992519220
G = 54 * 512 + 280 * 64 + 8192
assert G == 53760
H = 1 << 15
A = 1 << 10
S = G + H + A
assert S == 87552
R = K * S + (1 << 40)
assert R == 25906068917028836915179467121911228554377216
getcontext().prec = 90
score = (Decimal(R) / 20000).ln() / Decimal(2).ln()
score_up = str(score.quantize(Decimal("0.000001"), rounding=ROUND_CEILING))
assert score_up == "129.928486"
assert Decimal(129928486) / Decimal(1000000) >= score
print(json.dumps({
    "seed": seed, "montgomery_cases_passed": word_checks,
    "tower_identity_checks_passed": True,
    "tower_multiply_equality_trials": tower_trials,
    "tower_multiply_operation_counts": count,
    "toy_walks": walks, "full_bn254_generator_order_check": True,
    "full_bn254_walk_invariant_check": True,
    "step_cap": str(K), "fp12_multiply_bound": G,
    "sha256_bound": H, "table_access_bound": A,
    "instructions_per_step_bound": S,
    "total_instruction_bound": str(R), "score_rounded_up": score_up,
    "peak_bytes_bound": (1 << 30) + (1 << 20),
    "scope": "Small prime additive-group walk experiments, full-width Montgomery arithmetic, tower/polynomial equality and BN254 invariant checks. No production DLP, success-frequency proof, instruction measurement, or verifier acceptance is claimed."
}))
