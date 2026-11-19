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
assert p < 1 << 381 and p % 12 == 7 and r < 1 << 255
lam = p % r
assert pow(lam, 12, r) == 1 and all(pow(lam, k, r) != 1 for k in (1, 2, 3, 4, 6))
assert pow(lam, 6, r) == r - 1


# 1. Toy class walks with the exact frob_rho loop. Toy GT: the additive group
# Z_q with q = 1 mod 12, generator 1, and an automorphism y -> lam_q*y of
# order 12, playing Frobenius. Class key and frame are modelled by SHA-256:
# key = H(min orbit), frame = argmin_i H'(lam^i y). Both are class-equivariant
# like the production rule. DP rate 2^-2, gap 2^7.
def toy_prime(lo):
    while True:
        q = coins.randrange(lo, 2 * lo)
        q -= (q - 1) % 12
        if q > 13 and all(q % d for d in range(2, math.isqrt(q) + 1)):
            break
    while True:
        lq = pow(coins.randrange(2, q), (q - 1) // 12, q)
        if all(pow(lq, k, q) != 1 for k in (1, 2, 3, 4, 6)):
            return q, lq


walks = []
for lo in (10**5, 10**6, 10**7):
    q, lq = toy_prime(lo)
    lpq = [pow(lq, i, q) for i in range(12)]

    def canon(y, q=q, lpq=lpq):
        orbit = [lpq[i] * y % q for i in range(12)]
        i = min(range(12), key=lambda k: hashlib.sha256(b"F" + orbit[k].to_bytes(8, "big")).digest())
        key = int.from_bytes(hashlib.sha256(b"K" + min(orbit).to_bytes(8, "big")).digest(), "little")
        return i, key
    n = math.sqrt(q / 12)
    counts = {k: 0 for k in ("success", "identity", "degenerate", "cap", "table_full",
                             "verification_failed", "dp_event_limit", "escape_limit")}
    within = {"1.0": 0, "1.25": 0, "1.5": 0}
    steps_total = esc_total = 0
    for _ in range(trials):
        x = coins.randrange(q)
        got, st, outcome, esc = ref.frob_rho(
            q, lq, x, 0, lambda u, v, q=q: (u + v) % q, lambda t, q=q: t % q,
            lambda t, x=x, q=q: t * x % q, lambda i, y, q=q, lpq=lpq: lpq[i] * y % q,
            canon, coins, 2, 1 << 7)
        counts[outcome] += 1
        steps_total += st
        esc_total += esc
        if got is not None:
            assert got == x
            for c in within:
                if outcome == "success" and st <= float(c) * n:
                    within[c] += 1
    walks.append({"prime_order": q, "trials": trials, "outcomes": counts,
                  "sqrt_q_over_12": round(n, 1),
                  "mean_steps_over_sqrt_q_over_12": round(steps_total / trials / n, 3),
                  "fruitless_escapes_per_step": round(esc_total / max(1, steps_total), 5),
                  "success_within_c_sqrt_q_over_12_without_gap":
                      {c: within[c] / trials for c in within}})

# 2. Tower <-> protected basis consistency on random full-size elements.
one = (1,) + (0,) * 11
for _ in range(4):
    u = tuple(coins.randrange(p) for _ in range(12))
    v = tuple(coins.randrange(p) for _ in range(12))
    assert ref.tower_to_poly(ref.tower_mul(u, v, p), p) == \
        ref.fp12_mul(ref.tower_to_poly(u, p), ref.tower_to_poly(v, p), p)

# 3. Full-size GT: Frobenius formula, eigenvalue, class key invariance and
# frame equivariance on all 12 images.
tmul = lambda u, v: ref.tower_mul(u, v, p)
gam = ref.frob_constants(p)
while True:
    h = tuple(coins.randrange(p) for _ in range(12))
    g = ref.fp12_pow(ref.poly_to_tower(h, p), (p ** 12 - 1) // r, p, tmul)
    if g != one:
        break
assert ref.fp12_pow(g, r, p, tmul) == one


def gpow(el, n):
    n %= r
    return one if n == 0 else ref.fp12_pow(el, n, p, tmul)


assert ref.tower_frob(1, g, p, gam) == gpow(g, lam)          # pi(g) = g^p = g^(p mod r)
assert ref.tower_frob(1, h, p, gam) == ref.fp12_pow(h, p, p, tmul)  # formula on all of Fp12
equivariance = 0
for _ in range(4):
    y = gpow(g, coins.randrange(1, r))
    images = [ref.tower_frob(k, y, p, gam) for k in range(12)]
    assert len(set(images)) == 12
    canon_members = set()
    keys = set()
    for im in images:
        i, key = ref.production_canon(im, p, gam)
        canon_members.add(ref.tower_frob(i, im, p, gam))
        keys.add(key)
    assert len(canon_members) == 1 and len(keys) == 1
    equivariance += 1

# 4. Short real full-size class walk (production key/frame, Frobenius-twisted
# table). Reports index and frame histograms and the observed rate of 2-cycles
# (escapes expected about 1/3072 per step with 256 increments). Checks
# the label invariant at the end.
xs = coins.randrange(1, r)
z = gpow(g, xs)


def table2(el):          # el^(2^k), k < 255: fixed-base exponentiation
    out = [el]
    for _ in range(254):
        out.append(tmul(out[-1], out[-1]))
    return out


def fpow(tab, n):
    n %= r
    acc = one
    for k in range(255):
        if n >> k & 1:
            acc = tmul(acc, tab[k])
    return acc


gtab, ztab = table2(g), table2(z)
st = {}
_, nsteps, outcome, escs = ref.frob_rho(
    r, lam, z, one, tmul, lambda n: fpow(gtab, n), lambda n: fpow(ztab, n),
    lambda i, y: ref.tower_frob(i, y, p, gam), lambda y: ref.production_canon(y, p, gam),
    coins, 110, 1 << 115, step_limit=walk_steps, stats=st)
assert outcome == "cap"
Yf, af, bf = st["final"]
assert Yf == gpow(g, (af + bf * xs) % r)
ihist, fhist = st["index"], st["frame"]
chi2 = sum((v - sum(ihist) / 256) ** 2 / (sum(ihist) / 256) for v in ihist)

# 4b. Order of 2 modulo r (excludes escape cancellations in proof.md (iii)).
FACT = {2: 32, 3: 1, 11: 1, 19: 1, 10177: 1, 125527: 1, 859267: 1, 906349: 2,
        2508409: 1, 2529403: 1, 52437899: 1, 254760293: 2}


def is_prime(n):        # Miller-Rabin, first 12 prime bases: deterministic < 3.3e24
    if n < 2:
        return False
    bases = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    if n in bases:
        return True
    d, s = n - 1, 0
    while d % 2 == 0:
        d, s = d // 2, s + 1
    for a in bases:
        y = pow(a, d, n)
        if y in (1, n - 1):
            continue
        for _ in range(s - 1):
            y = y * y % n
            if y == n - 1:
                break
        else:
            return False
    return True


prod = 1
for fq, fe in FACT.items():
    assert is_prime(fq) and fq < 3 * 10 ** 24
    prod *= fq ** fe
assert prod == r - 1
ord2 = r - 1
for fq in FACT:
    while ord2 % fq == 0 and pow(2, ord2 // fq, r) == 1:
        ord2 //= fq
assert pow(2, ord2, r) == 1
ord2_ok = ord2 // math.gcd(ord2, 12) > 1 << 119

# 5. Score from the protected r.
D = 110
GAP = 1 << (D + 5)
S = ref.walk_cap(r, GAP)
E = 12288
F = 1 << 42
getcontext().prec = 100
W = Decimal(S) + Decimal(S * E + F) / 20000
score = W.ln() / Decimal(2).ln()
print(json.dumps({
    "seed": seed, "toy_class_walks": walks,
    "tower_basis_consistency_checks": 4,
    "frobenius_formula_and_eigenvalue_checked": True,
    "class_key_invariance_and_frame_equivariance_checks": equivariance,
    "full_size_class_walk": {"steps": nsteps, "index_histogram": ihist,
                             "index_chi2_255dof": round(chi2, 2),
                             "frame_histogram": fhist, "fruitless_escapes": escs,
                             "expected_escapes_about": round(walk_steps / 3072, 2),
                             "label_invariant_holds": True},
    "walk_cap_S": str(S), "dp_bits": D, "gap": str(GAP),
    "non_reference_instructions_per_step": E, "fixed_instructions": str(F),
    "total_units_W": str(W.quantize(Decimal(1), rounding=ROUND_CEILING)),
    "ord_r_2_log2": round(math.log2(ord2), 3), "ord_r_2_excludes_escape_cancellation": ord2_ok,
    "trap_case_Z_log2_per_start": round(ref.trap_bound_log2(256, 64, D), 2),
    "trap_case_Z_log2_total": round(ref.trap_bound_log2(256, 64, D) + math.log2(S), 2),
    "trap_case_N_log2_total": round(math.log2(S) + D + 1 - math.log2(r), 2),
    "trap_total_log2_upper": round(math.log2(2 ** (ref.trap_bound_log2(256, 64, D) + math.log2(S))
                                             + 2 ** (math.log2(S) + D + 1 - math.log2(r))), 2),
    "score_rounded_up": str(score.quantize(Decimal("0.000001"), rounding=ROUND_CEILING)),
    "scope": "Toy class walks with SHA-modelled key/frame, tower and Frobenius formula checks, class-key invariance and frame equivariance at full size, and a short real full-size class walk. No production DLP, success-frequency proof, pairing evaluation, instruction measurement or verifier acceptance is claimed."
}))
