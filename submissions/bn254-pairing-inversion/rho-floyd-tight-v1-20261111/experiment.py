"""Declared experiment: finite walk-statistics and BN254 invariant evidence.

Runs in the shared verifier's bounded, credential-free container (stdlib only,
20 s / 256 MiB / 1 CPU). It does NOT solve a production discrete logarithm,
prove a probability bound at the production order, or measure instructions. It:

  * checks four-word Montgomery arithmetic on full-width edge/random inputs;
  * runs the actual 32-adding Floyd walk on four small cyclic groups, reporting
    every outcome including failures;
  * checks, on the full-width protected BN254 constants, that the public
    generator g has order r (g^r = 1, g != 1) and that the exponent-label
    invariant holds for real Fp12 products;
  * recomputes the production cap, the itemized per-iteration instruction bound
    and the rounded-up log2 score from the protected order.
"""

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

# ---- Full-width four-word Montgomery arithmetic identity ---------------------
word_checks = 0
for m in (p, r):
    Rinv = pow(1 << 256, -1, m)
    cases = [(a, b) for a in (0, 1, m - 2, m - 1) for b in (0, 1, m - 2, m - 1)]
    cases += [(coins.randrange(m), coins.randrange(m)) for _ in range(256)]
    for a, b in cases:
        assert ref.montgomery_word_reference(a, b, m) == (a * b * Rinv) % m
        word_checks += 1

# ---- Small cyclic-group walks, including failures ----------------------------
walks = []
for q in (1009, 10007, 100003, 1000003):
    counts = {k: 0 for k in
              ("success", "identity", "degenerate", "cap", "verification_failed")}
    steps = 0
    for _ in range(trials):
        secret = coins.randrange(q)
        mul = lambda a, b, q=q: (a + b) % q
        scalar_mul = lambda e, k, q=q: (e * k) % q
        # Full 384-byte state length, as on BN254; only small cyclic groups here.
        encode = lambda e: e.to_bytes(384, "big")
        got, n, outcome = ref.rho(q, secret, mul, scalar_mul, 1, secret % q, 0,
                                  encode, coins)
        counts[outcome] += 1
        steps += n
        if got is not None:
            assert got == secret
    walks.append({
        "prime_order": q, "trials": trials, "outcomes": counts,
        "success_fraction": (counts["success"] + counts["identity"]) / trials,
        "mean_floyd_iterations": steps / trials,
    })

# ---- Full-width BN254 invariant checks (not a production DLP run) -------------
g = tuple(int(v, 16) for v in profile["generator_g"]["coefficients"])
one = (1,) + (0,) * 11
generator_order_ok = (ref.fp12_pow(g, r, p) == one) and (g != one)
assert generator_order_ok

x = coins.randrange(r)
z = ref.fp12_pow(g, x, p)
a, b, c, d = [coins.randrange(r) for _ in range(4)]
y = ref.fp12_mul(ref.fp12_pow(g, a, p), ref.fp12_pow(z, b, p), p)
M = ref.fp12_mul(ref.fp12_pow(g, c, p), ref.fp12_pow(z, d, p), p)
for _ in range(4):
    y = ref.fp12_mul(y, M, p)
    a, b = (a + c) % r, (b + d) % r
walk_invariant_ok = (y == ref.fp12_pow(g, (a + b * x) % r, p))
assert walk_invariant_ok

# ---- Production ledger and score, recomputed from the protected order --------
N, per_iter, R = ref.ledger(r)
getcontext().prec = 90
score = (Decimal(R) / 20000).ln() / Decimal(2).ln()
score_up = str(score.quantize(Decimal("0.000001"), rounding=ROUND_CEILING))

print(json.dumps({
    "seed": seed,
    "montgomery_cases_passed": word_checks,
    "toy_walks": walks,
    "full_bn254_generator_order_check": generator_order_ok,
    "full_bn254_walk_invariant_check": walk_invariant_ok,
    "floyd_cap": str(N),
    "per_group_step_bound": ref.STEP,
    "instructions_per_floyd_iteration_bound": per_iter,
    "step_breakdown": {"fp12_mul": ref.FP12_MUL, "sha256_select": ref.SHA_SELECT,
                       "exponent_update": ref.UPDATE, "iter_compare": ref.COMPARE,
                       "setup_allowance": ref.SETUP},
    "total_instruction_bound": str(R),
    "score_rounded_up": score_up,
    "peak_bytes_bound": 1 << 20,
    "scope": ("Small cyclic-group walk experiments, full-width Montgomery "
              "arithmetic and BN254 invariant checks. No production DLP, "
              "success-frequency proof, instruction measurement, or verifier "
              "acceptance is claimed."),
}))
