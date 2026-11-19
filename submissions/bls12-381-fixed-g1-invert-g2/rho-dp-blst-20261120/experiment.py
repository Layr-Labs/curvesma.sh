import hashlib
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
walk_steps = request["full_size_walk_steps"]
coins = random.Random(seed)
profile = json.loads((base / "target-profile.json").read_text())
p = int(profile["field_modulus"], 16)
r = int(profile["subgroup_order"], 16)
assert p < 1 << 381 and p % 4 == 3 and r < 1 << 255
assert (1 << 384) < 10 * p  # <=9 conditional subtractions canonicalize 6 words

# 1. Toy walks: the exact rho_dp loop on small prime additive groups. The
# coefficient bits of production GT elements are modelled by SHA-256 of the
# toy element (index word, DP word, fingerprint); this models, not proves,
# the production partition. DP rate 2^-2, gap 2^7 = 32/theta as in production.
walks = []
for q in (10007, 100003, 1000003, 10000019):
    counts = {k: 0 for k in ("success", "identity", "degenerate", "cap",
                             "table_full", "verification_failed")}
    within = 0
    steps = 0
    for _ in range(trials):
        secret = coins.randrange(q)
        mul = lambda a, b, q=q: (a + b) % q
        power = lambda a, target=1, q=q: a * target % q

        def key(y):
            h = hashlib.sha256(y.to_bytes(8, "big")).digest()
            return (int.from_bytes(h[0:8], "little"), int.from_bytes(h[8:24], "little"),
                    int.from_bytes(h[24:32], "little"))
        got, n, outcome = ref.rho_dp(q, secret, 0, mul, power, key, coins, 2, 1 << 7)
        counts[outcome] += 1
        steps += n
        if got is not None:
            assert got == secret
            if n <= (3 * (math.isqrt(q) + 1)) // 2:
                within += 1
    walks.append({"prime_order": q, "trials": trials, "outcomes": counts,
                  "success_fraction": (counts["success"] + counts["identity"]) / trials,
                  "success_within_1.5sqrt_q_without_gap": within / trials,
                  "mean_steps": steps / trials})

# 2. Tower <-> protected basis consistency on random full-size elements.
one_poly = (1,) + (0,) * 11
for _ in range(8):
    x = tuple(coins.randrange(p) for _ in range(12))
    y = tuple(coins.randrange(p) for _ in range(12))
    assert ref.poly_to_tower(ref.tower_to_poly(x, p), p) == x
    assert ref.tower_to_poly(ref.tower_mul(x, y, p), p) == \
        ref.fp12_mul(ref.tower_to_poly(x, p), ref.tower_to_poly(y, p), p)
assert ref.poly_to_tower(one_poly, p) == one_poly

# 3. Full-size GT: order-r element g = h^((p^12-1)/r), z = g^x, and a real
# production walk of walk_steps steps in the tower representation with the
# production key. Reports partition / low-bit statistics and checks the label
# invariant Y = g^a z^b at the end.
tmul = lambda u, v: ref.tower_mul(u, v, p)
while True:
    h = tuple(coins.randrange(p) for _ in range(12))
    g = ref.fp12_pow(ref.poly_to_tower(h, p), (p ** 12 - 1) // r, p, tmul)
    if g != one_poly:
        break
assert ref.fp12_pow(g, r, p, tmul) == one_poly


def gpow(base_el, n):
    n %= r
    return one_poly if n == 0 else ref.fp12_pow(base_el, n, p, tmul)


xs = coins.randrange(1, r)
z = gpow(g, xs)
tbl = []
for _ in range(32):
    c, d = coins.randrange(r), coins.randrange(r)
    tbl.append((tmul(gpow(g, c), gpow(z, d)), c, d))
a, b = coins.randrange(r), coins.randrange(r)
Y = tmul(gpow(g, a), gpow(z, b))
hist = [0] * 32
low4_zero = 0
for _ in range(walk_steps):
    k = ref.production_key(Y, p)
    j = k[0] & 31
    hist[j] += 1
    m, c, d = tbl[j]
    Y = tmul(Y, m)
    a, b = (a + c) % r, (b + d) % r
    if ref.production_key(Y, p)[1] & 15 == 0:
        low4_zero += 1
assert Y == gpow(g, (a + b * xs) % r)
expected = walk_steps / 32
chi2 = sum((v - expected) ** 2 / expected for v in hist)

# 4. Score from the protected r.
D = 110
GAP = 1 << (D + 5)
S = ref.walk_cap(r, GAP)
E = 8192                      # non-reference instructions per step (machine.md)
F = 1 << 40                   # all setup / DP handling / output instructions
getcontext().prec = 100
W = Decimal(S) + Decimal(S * E + F) / 20000
score = W.ln() / Decimal(2).ln()
score_up = str(score.quantize(Decimal("0.000001"), rounding=ROUND_CEILING))
print(json.dumps({
    "seed": seed, "toy_dp_walks": walks,
    "tower_basis_consistency_checks": 8,
    "full_size_gt_walk": {"steps": walk_steps, "partition_histogram": hist,
                          "partition_chi2_31dof": round(chi2, 2),
                          "dp_word_low4_zero_fraction": low4_zero / walk_steps,
                          "label_invariant_holds": True},
    "walk_cap_S": str(S), "dp_bits": D, "gap": str(GAP),
    "non_reference_instructions_per_step": E, "fixed_instructions": str(F),
    "total_units_W": str(W.quantize(Decimal(1), rounding=ROUND_CEILING)),
    "score_rounded_up": score_up,
    "scope": "Toy additive-group DP walks with a SHA model of coefficient bits, tower/protected-basis consistency, and a short real full-size GT walk with partition statistics. No production DLP, success-frequency proof, pairing evaluation, instruction measurement or verifier acceptance is claimed."
}))
