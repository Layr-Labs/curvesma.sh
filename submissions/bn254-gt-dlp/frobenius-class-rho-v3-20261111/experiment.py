import importlib.util
import json
import math
import random
import sys
from decimal import Decimal, ROUND_CEILING, getcontext
from fractions import Fraction
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
g = tuple(int(v, 16) for v in profile["generator_g"]["coefficients"])
one = (1,) + (0,) * 11
out = {"seed": seed, "trials": trials}

# 1. Word-level Montgomery multiplication (same check as the baseline).
word_checks = 0
for m in (p, r):
    Rinv = pow(1 << 256, -1, m)
    cases = [(a, b) for a in (0, 1, m - 2, m - 1) for b in (0, 1, m - 2, m - 1)]
    cases += [(coins.randrange(m), coins.randrange(m)) for _ in range(256)]
    for a, b in cases:
        assert ref.montgomery_word_reference(a, b, m) == (a * b * Rinv) % m
        word_checks += 1
out["montgomery_cases_passed"] = word_checks

# 2. Frobenius structure on the protected Fp12 basis and GT.
assert ref.fp12_pow(g, r, p) == one and g != one
mats = ref.frobenius_matrices(p)
lam = [pow(p, k, r) for k in range(12)]
assert len(set(lam)) == 12 and lam[6] == r - 1 and pow(p, 12, r) == 1
assert all(pow(p, k, r) != 1 for k in (1, 2, 3, 4, 6))
assert (r - 1) % 12 == 0
ident = [[int(i == j) for j in range(12)] for i in range(12)]
assert [[sum(mats[1][i][t] * mats[11][t][j] for t in range(12)) % p for j in range(12)]
        for i in range(12)] == ident
for k in range(12):
    for i in range(12):
        for j in range(12):
            if mats[k][i][j]:
                assert i % 6 == j % 6          # 2x2 blocks on {i, i+6}
assert all(mats[6][i][j] == ((-1) ** i % p if i == j else 0) for i in range(12) for j in range(12))
for k in range(12):                            # Frob^k(g) = g^(p^k mod r) on GT
    assert ref.apply(mats[k], g, p) == ref.fp12_pow(g, lam[k], p)
nnz = [sum(1 for i in range(12) for j in range(12) if mats[k][i][j]) for k in range(12)]
key_mm = sum(1 for k in range(1, 6) for i in (1, 7) for j in range(12) if mats[k][i][j])
key_add = sum(1 for k in range(1, 6) for i in (1, 7)
              if sum(1 for j in range(12) if mats[k][i][j]) == 2)
full_mm = max(nnz)
full_add = max(sum(max(0, sum(1 for j in range(12) if mats[k][i][j]) - 1) for i in range(12))
               for k in range(12))
out["frobenius"] = {"nonzeros_per_k": nnz, "key_mults_k1_to_5": key_mm,
                    "key_adds": key_add, "full_apply_mults_max": full_mm,
                    "full_apply_adds_max": full_add, "lambda6_is_minus_one": True,
                    "frob_k_of_g_equals_g_pow_p_k_mod_r_all_k": True}

# 3. Canonical representative: class invariance and label consistency on GT.
for _ in range(3):
    y = ref.fp12_pow(g, coins.randrange(1, r), p)
    rep, k, tie = ref.canon_fp12(y, mats, p)
    assert rep == ref.apply(mats[k], y, p) and not tie
    for j in range(12):
        assert ref.canon_fp12(ref.apply(mats[j], y, p), mats, p)[0] == rep
assert rep == ref.fp12_pow(y, lam[k], p)
assert ref.canon_fp12(one, mats, p)[:2] == (one, 0)       # identity: all 12 conjugates equal
out["canonical_class_invariance_checked_elements"] = 3

# 4. Multiplicative order of 2 mod r from the verified factorisation of r-1.
FACT = {2: 28, 3: 2, 13: 1, 29: 1, 983: 1, 11003: 1, 237073: 1, 405928799: 1,
        1670836401704629: 1, 13818364434197438864469338081: 1}


def probable_prime(n):
    if n < 4096:
        return n > 1 and all(n % d for d in range(2, math.isqrt(n) + 1))
    d, s = n - 1, 0
    while d % 2 == 0:
        d, s = d // 2, s + 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61):
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


assert math.prod(q ** e for q, e in FACT.items()) == r - 1
assert all(probable_prime(q) for q in FACT)
ord2 = r - 1
for q, e in FACT.items():
    for _ in range(e):
        if ord2 % q == 0 and pow(2, ord2 // q, r) == 1:
            ord2 //= q
assert pow(2, ord2, r) == 1

# 5. Production-field class walk: invariants only (tiny table, 24 steps).
ops = ref.bn254_ops(p, r, mats)
x = coins.randrange(1, r)
z = ref.fp12_pow(g, x, p)
st = ref.class_rho(r, g, z, ops, coins, cap=24, table_size=4)
yf, af, bf = st["final"]
assert st["outcome"] == "cap" and st["steps"] == 24
assert yf == ref.fp12_mul(ref.fp12_pow(g, af, p), ref.fp12_pow(z, bf, p), p)
assert ref.canon_fp12(yf, mats, p)[0] == yf
out["full_bn254_class_walk_invariant_check"] = {"steps": 24, "ties": st["ties"]}


# 6. Random-map prediction for the stack detector: T = rho*(1+v*w).
def q_tail(y):
    if y <= 1:
        return 0.0
    if y >= 2:
        return 1.0
    t = y - 1
    return t * (1 - math.log(t))


def predicted_success(c, steps=20000, top=12.0):
    h = top / steps
    return sum(((i + .5) * h) * math.exp(-((i + .5) * h) ** 2 / 2) * q_tail(c / ((i + .5) * h)) * h
               for i in range(steps))


MEAN_T = 1.25 * math.sqrt(math.pi / 2)         # E[rho*(1+vw)]/sqrt(n_eff)
C0 = Fraction(9, 5)
DELTA = Fraction(1, 64)


# 7. Toy walks on prime additive groups with an order-12 automorphism group.
def small_prime(n):
    return n > 1 and all(n % d for d in range(2, math.isqrt(n) + 1))


def prime_1mod12(n):
    x = n + (1 - n) % 12
    while not small_prime(x):
        x += 12
    return x


def run_toy(q, m, table_size, n_trials):
    if m == 12:
        h = 2
        while True:
            l = pow(h, (q - 1) // 12, q)
            if pow(l, 4, q) != 1 and pow(l, 6, q) != 1:
                break
            h += 1
        lam_t = [pow(l, k, q) for k in range(12)]
        n = (q - 1) // 12 + 1
    else:
        lam_t, n = [1], q
    ops_t = ref.toy_ops(q, lam_t)
    n_eff = n * table_size / (table_size - 1)
    cap_c0 = math.ceil(float(C0) * math.sqrt(n_eff))
    big = math.ceil(20 * math.sqrt(n_eff)) + 100
    counts, solved, within, esc, two, waste, depth = {}, [], 0, 0, 0, 0, 0
    total_steps = 0
    for _ in range(n_trials):
        secret = coins.randrange(q)
        s = ref.class_rho(q, 1, secret, ops_t, coins, cap=big, table_size=table_size)
        counts[s["outcome"]] = counts.get(s["outcome"], 0) + 1
        if s["outcome"] in ("success", "identity"):
            assert s["x"] == secret
            solved.append(s["steps"])
            within += s["steps"] <= cap_c0
        esc += s["escapes"]
        two += s["two_cycles"]
        waste += s["wasted_steps"]
        depth = max(depth, s["max_depth"])
        total_steps += s["steps"]
    mean = sum(solved) / max(1, len(solved))
    return {"q": q, "automorphism_order": m, "classes": n, "table_size": table_size,
            "trials": n_trials, "outcomes": counts,
            "cap_steps_at_c0": cap_c0, "success_within_c0_cap": within / n_trials,
            "predicted_success_at_c0": round(predicted_success(float(C0)), 4),
            "mean_steps_to_solve": round(mean, 1),
            "predicted_mean_1.2533*1.25*sqrt(n_eff)": round(MEAN_T * math.sqrt(n_eff), 1),
            "measured_over_predicted_mean": round(mean / (MEAN_T * math.sqrt(n_eff)), 4),
            "escapes_per_step": esc / max(1, total_steps),
            "two_cycles_per_step": two / max(1, total_steps),
            "predicted_two_cycle_rate_1/(m*R)": (1 / (m * table_size)) if m > 1 else 0,
            "wasted_step_fraction": waste / max(1, total_steps), "max_stack_depth": depth}


toys = [prime_1mod12(10 ** e) for e in (4, 5, 6, 7)]
class_walks = [run_toy(q, 12, ref.R_T, trials) for q in toys]
plain_walks = [run_toy(q, 1, ref.R_T, trials) for q in toys[:3]]
fruitless_small_table = run_toy(toys[2], 12, 32, trials)
out["toy_class_walks"] = class_walks
out["toy_plain_walks_same_detector"] = plain_walks
out["toy_class_walk_32_table_fruitless_stress"] = fruitless_small_table
out["plain_over_class_mean_steps_ratio"] = [
    round(pw["mean_steps_to_solve"] / cw["mean_steps_to_solve"], 3)
    for pw, cw in zip(plain_walks, class_walks)]

# 8. Production cost chain (scalar AArch64 instruction-equivalent bounds).
MM, AS, WORD = 4096, 256, 8                     # baseline primitive bounds
G, H, A = 1 << 20, 1 << 16, 1 << 12             # Fp12 product/square, SHA+encode, table record
neg, key_cmp = 12, 22
K_tally = ((key_mm + full_mm + 2) * MM + (key_add + full_add + neg + key_cmp) * AS
           + 48 * 4 * WORD + (1 << 14))
K = 1 << 18
assert K_tally <= K
S_tally = 2 * 12 * AS + 56 * WORD + AS + AS     # 2 amortised order compares, push, degeneracy, depth
S = 1 << 13
assert S_tally <= S
CTRL = 1 << 12
I_step = G + H + A + K + S + CTRL
F = 1 << 42
n_cls = (r - 1) // 12 + 1
X = Fraction(C0.numerator ** 2 * n_cls * ref.R_T, C0.denominator ** 2 * (ref.R_T - 1))
E = math.isqrt(X.numerator // X.denominator) + 1
assert E * E >= X and E < 1 << 128
assert ord2 > 12 * E
# Lower bound not relying on primality of the 94-bit factor: the true order is a
# multiple of ord2 / 13818364434197438864469338081.
assert ord2 // 13818364434197438864469338081 > 12 * E
# Setup/final/output F, itemised (same formulas as machine.md), and storage.
EXP = 508 * G + (1 << 20)                        # one GT exponentiation (<=508 products) incl. bit handling
F_items = {
    "parse_validate_z^r": EXP + (1 << 32),
    "randomness_bits_x64": (2 * 1025 * 1024 * 256 + 256) * 64,
    "table_1024_plus_start": 1025 * (2 * EXP + G),
    "start_canonicalise_push_extra_compare": K + S,
    "frobenius_constants": EXP + 11 * 1728 * (MM + AS) + 11 * MM + 12 * 23 * MM,
    "montgomery_constants": 2 * (6 * 64 + 512 * AS),
    "tie_budget_1024_events": 1024 * (12 * (23 * MM + 11 * AS) + 11 * 12 * AS + (1 << 16)),
    "final_inverse_division_verification": (1 << 22) + MM + EXP + (1 << 16),
    "pairing_output_G2_double_and_add_and_check": 1016 * (1 << 24) + (1 << 30),
    "allocation_init_fixed_overhead": 4097 * 64 * 8 + (1 << 32),
}
F_tally = sum(F_items.values())
assert F_tally <= F
storage = {"table": 1024 * 448, "stack": 4097 * 512, "frobenius_constants": 12 * 23 * 32 + 12 * 32,
           "baseline_allocation_state_scratch_code": 1 << 20}
assert sum(storage.values()) < 1 << 22
R_total = E * I_step + F
getcontext().prec = 90
log2W = (Decimal(R_total) / 20000).ln() / Decimal(2).ln()
computed_up = log2W.quantize(Decimal("0.000001"), rounding=ROUND_CEILING)
MARGIN = Decimal("0.020000")
claim = computed_up + MARGIN
p_c0 = predicted_success(float(C0))
p_eff = predicted_success(float(C0 * (1 - DELTA)))
out["cost"] = {"F_items": F_items, "F_tally": F_tally, "log2_F_tally": round(math.log2(F_tally), 3),
               "storage_bytes": storage, "storage_total": sum(storage.values()),
               "K_canonicalisation_tally": K_tally, "K_bound": K, "S_stack_tally": S_tally,
               "S_bound": S, "I_step": I_step, "F_setup_and_output": F,
               "classes_n": str(n_cls), "step_cap_E": str(E), "log2_E": round(math.log2(E), 6),
               "log2_ord_r_2": round(math.log2(ord2), 3), "ord_r_2_exceeds_12E": True,
               "total_instruction_bound": str(R_total), "computed_log2W_rounded_up": str(computed_up),
               "stated_margin": str(MARGIN), "claimed_score": str(claim),
               "baseline_score": "135.216266",
               "improvement_bits": str(Decimal("135.216266") - claim)}
out["success_model"] = {"c0": "9/5", "waste_allowance_delta": "1/64",
                        "predicted_P_T_le_c0": round(p_c0, 4),
                        "predicted_P_T_le_c0_minus_waste": round(p_eff, 4),
                        "log2_expected_tie_events": round(math.log2(E * 66) - math.log2(p), 1),
                        "stack_depth_cap": ref.D_MAX, "claimed_success_floor": 0.5}
# 9. Production-field escape path (own coin stream; 1-entry table forces fruitless 2-cycles).
ecoins = random.Random(seed ^ 0x5EED)
xe = ecoins.randrange(1, r)
ze = ref.fp12_pow(g, xe, p)
se = ref.class_rho(r, g, ze, ops, ecoins, cap=48, table_size=1)
ye, ae, be = se["final"]
assert se["outcome"] in ("cap", "success")
assert se["outcome"] != "success" or se["x"] == xe
assert ye == ref.fp12_mul(ref.fp12_pow(g, ae, p), ref.fp12_pow(ze, be, p), p)
out["full_bn254_escape_path_check"] = {"steps": se["steps"], "outcome": se["outcome"],
                                       "escapes": se["escapes"], "two_cycles": se["two_cycles"]}
out["scope"] = ("Small additive prime-group walks with an order-12 automorphism group, word-level "
                "Montgomery checks, full BN254 Frobenius/canonicalisation/walk invariants and the "
                "complete capped cost. No production DLP, success-frequency proof, instruction "
                "measurement or verifier acceptance is claimed.")
print(json.dumps(out))
