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
spec = importlib.util.spec_from_file_location("fdp_ref", base / "reference.py")
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)
request = json.load(sys.stdin)
seed = request["seed"]
trials = request["trials"]
coins = random.Random(seed)
profile = json.loads((base / "target-profile.json").read_text())
p = int(profile["field_modulus"], 16)
r = int(profile["subgroup_order"], 16)
g = tuple(int(v, 16) for v in profile["generator_g"]["coefficients"])

# ---- 1. Full-width structural facts used by the argument -------------------
assert pow(p, 12, r) == 1 and pow(p, 6, r) == r - 1
assert pow(p, 4, r) != 1 and pow(p, 6, r) != 1          # exact order 12
assert all((pow(p, k, r) - 1) % r != 0 for k in range(1, 12))  # all orbits size 12
assert (r - 1) % 12 == 0 and p % 6 == 1
G = ref.Fp12Group(p, r, g)
assert all(G.cj[j] != 0 for j in ref.SELECT) and G.cj[6] == 1
assert G.power(g, r) == G.identity and g != G.identity

frobenius_checks = 0
for _ in range(2):
    k = coins.randrange(1, r)
    Y = G.power(g, k)
    assert G.frob(Y, 1) == G.power(Y, p)                 # linear map == p-th power
    conj = Y
    for j in range(1, 13):
        conj = G.frob(conj, 1)
        assert conj == G.frob(Y, j % 12)
        assert conj == G.power(g, k * pow(p, j, r) % r)  # acts as *p^j on logs
        assert G.trace_twelfth(conj) == G.trace_twelfth(Y)   # class invariant
        frobenius_checks += 1
    assert conj == Y

# ---- 2. Full-width planted conjugate walks: invariants and recovery --------
x = coins.randrange(r)
z = G.power(g, x)
salt = coins.getrandbits(128).to_bytes(16, "big")
a1, b1, b2 = (coins.randrange(r) for _ in range(3))
k0 = coins.randrange(1, 12)
a2 = (pow(p, k0, r) * (a1 + b1 * x) - b2 * x) % r         # test-harness planting only
Y1 = G.mul(G.power(g, a1), G.power(z, b1))
Y2 = G.mul(G.power(g, a2), G.power(z, b2))
assert Y2 == G.frob(Y1, k0)
C1 = C2 = 1
selected = []
for _ in range(6):
    selected.append(ref.selector(ref.digest(salt, G.invariant(Y1))))
    assert G.invariant(Y1) == G.invariant(Y2)
    Y1, C1 = ref.step(G, salt, Y1, C1)
    Y2, C2 = ref.step(G, salt, Y2, C2)
    assert Y1 == G.power(g, C1 * (a1 + b1 * x) % r)       # label invariant
    assert Y2 == G.frob(Y1, k0)                           # walk is defined on classes
kk = [k for k in range(12) if G.frob(Y1, k) == Y2]
assert kk == [k0]
pk = pow(p, k0, r)
den = (C2 * b2 - pk * C1 * b1) % r
recovered = (pk * C1 * a1 - C2 * a2) * pow(den, -1, r) % r
assert recovered == x and G.power(g, recovered) == z

# ---- 3. Toy runs of the complete algorithm (order-12 automorphism) ---------
def order12_element(q, rng):
    while True:
        rho = pow(rng.randrange(2, q), (q - 1) // 12, q)
        if pow(rho, 6, q) == q - 1 and pow(rho, 4, q) != 1:
            return rho

toy = []
for q, theta_bits in ((120097, 2), (1200061, 3), (12000097, 4)):
    assert q % 12 == 1 and all(q % d for d in range(2, math.isqrt(q) + 1))
    T = ref.ToyGroup(q, order12_element(q, coins))
    N = (q - 1) // 12
    cap = math.isqrt(3 * (q - 1) // 16) + 1               # ceil(1.5*sqrt(N)), as in production
    counts = {k: 0 for k in ("success", "identity", "degenerate", "cap",
                             "walk_cap", "verification_failed")}
    uncapped = 0
    n_trials = trials if q < 10**7 else max(16, trials // 4)
    for _ in range(n_trials):
        secret = coins.randrange(q)
        got, n, outcome = ref.solve(T, secret, theta_bits, cap, 1 << 20,
                                    64 << theta_bits, coins)
        counts[outcome] += 1
        if got is not None:
            assert got == secret
        got, n, outcome = ref.solve(T, secret, theta_bits, 1 << 40, 1 << 20,
                                    64 << theta_bits, coins)
        assert got == secret
        uncapped += n
    toy.append({
        "prime_order": q, "classes_N": N, "theta": "2^-%d" % theta_bits,
        "step_cap": cap, "trials": n_trials, "outcomes_at_cap": counts,
        "success_fraction_at_cap": (counts["success"] + counts["identity"]) / n_trials,
        "uncapped_mean_steps_over_sqrt_piN_over_2":
            (uncapped / n_trials) / math.sqrt(math.pi * N / 2),
    })

# ---- 4. Score reproduction ---------------------------------------------------
S = math.isqrt(3 * (r - 1) // 16) + 1                     # >= 1.5*sqrt((r-1)/12)
M, A = 256, 32
fp12_mul = 54 * M + 280 * A
fp12_sqr = 36 * M + 210 * A
frobenius = 15 * M + 31 * A
sha_block = 3200
scalar_mul = 256
control = 992
per_step = fp12_mul + fp12_sqr + frobenius + sha_block + scalar_mul + control
assert per_step == 48000
R = S * per_step + (1 << 64)
getcontext().prec = 90
score = (Decimal(R) / 20000).ln() / Decimal(2).ln()
score_up = str(score.quantize(Decimal("0.000001"), rounding=ROUND_CEILING))
print(json.dumps({
    "seed": seed,
    "frobenius_order_12_and_trace_invariance_checks": frobenius_checks,
    "full_width_planted_conjugate_walk": {
        "steps": 6, "selected_frobenius_powers": selected,
        "planted_k": k0, "recovered_scalar_ok": True},
    "toy_runs": toy,
    "heuristic_prediction_at_cap": 1 - math.exp(-(10 / 11) * 1.125),
    "step_cap_S": str(S), "instructions_per_step_bound": per_step,
    "total_instruction_bound": str(R), "score_rounded_up": score_up,
    "scope": "Full-width BN254 checks of the Frobenius/trace/label identities and of "
             "the recovery formula on a planted conjugate pair; complete-algorithm runs "
             "only on toy prime-order additive groups with an order-12 automorphism. "
             "No production DLP, instruction measurement or acceptance is claimed.",
}))
