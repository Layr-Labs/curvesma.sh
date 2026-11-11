"""Hosted audit: tower facts, exact production-size norms, stored-sample replay and the complete ledger.

Reads {"seed": int, "fresh_cells": int} on stdin. Python standard library only.
Scope: consistency of the submitted model with exact arithmetic at production parameters.
It does not run a production attack and does not prove the heuristics of proof.md.
"""
import importlib.util
import json
import math
import random
import sys
from pathlib import Path

base = Path("/snapshot/files")
if not base.is_dir():
    base = Path(__file__).parent
spec = importlib.util.spec_from_file_location("tnfs", base / "tnfs.py")
T = importlib.util.module_from_spec(spec)
spec.loader.exec_module(T)

request = json.load(sys.stdin)
seed = int(request.get("seed", 20261111))
fresh = int(request.get("fresh_cells", 2000))
rng = random.Random(seed)
profile = json.loads((base / "target-profile.json").read_text())
data = T.load_data(base)
p, r, u = T.P, T.R_ORDER, T.U
checks = {}

# 1. protected constants and tower
checks["p_matches_profile"] = p == int(profile["field_modulus"], 16)
checks["r_matches_profile"] = r == int(profile["subgroup_order"], 16)
checks["u_matches_profile"] = u == int(profile["parameter_x"], 16)
checks["r_divides_p4_minus_p2_plus_1"] = (p**4 - p**2 + 1) % r == 0
checks["r_coprime_to_p6_minus_1"] = math.gcd(r, p**6 - 1) == 1
y = [0, 1, 0, 0, 0, 0]


def frob_power(k):
    return T.hpow(y, p**k, p)


# Rabin: h irreducible mod p iff y^(p^6) = y and gcd(y^(p^k) - y, h) = 1 for k = 2, 3.
def poly_gcd_is_one(a, hmod, q):
    def trim(v):
        v = [c % q for c in v]
        while v and v[-1] == 0:
            v.pop()
        return v

    a, b = trim(hmod), trim(a)
    while b:
        inv = pow(b[-1], -1, q)
        while len(a) >= len(b):
            c = a[-1] * inv % q
            sh = len(a) - len(b)
            for i, bc in enumerate(b):
                a[sh + i] = (a[sh + i] - c * bc) % q
            a = trim(a)
        a, b = b, a
    return len(a) == 1


f6 = frob_power(6)
ok = f6 == y
for k in (2, 3):
    d = frob_power(k)
    d[1] = (d[1] - 1) % p
    ok = ok and poly_gcd_is_one(d, T.H, p)
checks["h_irreducible_mod_p"] = ok
uy = [0, u % p, 0, 0, 0, 0]
checks["u_y_is_a_nonsquare_in_Fp6"] = T.hpow(uy, (p**6 - 1) // 2, p) == [p - 1, 0, 0, 0, 0, 0]
# f(x) at x^2 = u y equals y^4 (36u^4 + 36u^3 + 24u^2 + 6u + 1) = y^4 p = 0 mod p
checks["f_and_g_share_the_root_mod_p"] = (36 * u**4 + 36 * u**3 + 24 * u**2 + 6 * u + 1) % p == 0
units = {"y": ([0, 1, 0, 0, 0, 0], [-2, 0, 0, 1, 0, 1]),
         "e2": ([-1, 0, 1, 1, 0, 1], [-1, 0, 0, 1, 0, 0]),
         "e3": ([-1, 1, 0, 1, 0, 1], [-1, 0, 1, -1, 1, 0])}
checks["units_and_inverses"] = all(T.hmul(a, b) == [1, 0, 0, 0, 0, 0] for a, b in units.values())
checks["norm_of_y_is_minus_one"] = T.hnorm(y) == -1

# 2. replay of stored production-size cells (exact integer norms)
main = data["main"]
A = data["config"]["A"]
worst = 0.0
for (a, b), row in zip(main["check_vectors"], main["rows"]):
    row = [row[0] / 10, row[1] / 10, row[2], 0, 0, row[5] / 10, row[6] / 10]
    nf, ng = T.norm_f(a, b), T.norm_g(a, b)
    worst = max(worst, abs(T.log2_int(nf) - row[0]), abs(T.log2_int(ng) - row[1]),
                abs(T.smooth_part_bits(nf) - row[5]), abs(T.smooth_part_bits(ng) - row[6]))
checks["stored_cells_replayed"] = len(main["check_vectors"])
checks["stored_cells_max_abs_difference_bits"] = round(worst, 4)

# unit multiplicities of the first stored cells (exponents of y, e2, e3 in [-4, 4])
pw = {}
for name, (e, ei) in units.items():
    tab = {0: [1, 0, 0, 0, 0, 0]}
    for k in range(1, 5):
        tab[k] = T.hmul(tab[k - 1], e)
        tab[-k] = T.hmul(tab[1 - k], ei)
    pw[name] = tab
unit_list = []
for i in range(-4, 5):
    for j in range(-4, 5):
        yi = T.hmul(pw["y"][i], pw["e2"][j])
        for k in range(-4, 5):
            if i or j or k:
                unit_list.append(T.hmul(yi, pw["e3"][k]))
mult_ok = True
for (a, b), row in list(zip(main["check_vectors"], main["rows"]))[:12]:
    m = 1
    for e in unit_list:
        if max(abs(c) for c in T.hmul(e, a)) <= A and max(abs(c) for c in T.hmul(e, b)) <= A:
            m += 1
    mult_ok = mult_ok and m == row[2]
checks["unit_multiplicities_replayed"] = mult_ok
# exponent range [-6, 6] on four cells: no unit outside [-4, 4] keeps a cell in the box
for name, (e_, ei_) in units.items():
    for k in (5, 6):
        pw[name][k] = T.hmul(pw[name][k - 1], e_)
        pw[name][-k] = T.hmul(pw[name][1 - k], ei_)
wide_ok = True
for (a, b), row in list(zip(main["check_vectors"], main["rows"]))[:4]:
    m = 0
    for i in range(-6, 7):
        for j in range(-6, 7):
            yi = T.hmul(pw["y"][i], pw["e2"][j])
            for k in range(-6, 7):
                e_ = T.hmul(yi, pw["e3"][k])
                if max(abs(c) for c in T.hmul(e_, a)) <= A and max(abs(c) for c in T.hmul(e_, b)) <= A:
                    m += 1
    wide_ok = wide_ok and m == row[2]
checks["unit_multiplicities_unchanged_at_exponents_up_to_6"] = wide_ok

# 3. fresh uniform cells of the production box (not filtered for primitivity)
fn, gn, fs, gs = [], [], [], []
for _ in range(fresh):
    a = [rng.randint(-A, A) for _ in range(6)]
    b = [rng.randint(-A, A) for _ in range(6)]
    if not any(a) or not any(b):
        continue
    nf, ng = T.norm_f(a, b), T.norm_g(a, b)
    fn.append(T.log2_int(nf))
    gn.append(T.log2_int(ng))
    fs.append(T.smooth_part_bits(nf))
    gs.append(T.smooth_part_bits(ng))
fresh_out = {}
if len(fn) > 1:
    for name, v in (("log2_norm_f", fn), ("log2_norm_g", gn), ("smooth_part_bits_f", fs), ("smooth_part_bits_g", gs)):
        m, se = T.mean_se(v)
        fresh_out[name] = {"mean": round(m, 3), "se": round(se, 3)}

# 4. stored special-q lattices: every basis vector is in the lattice and the determinant is q
lat_ok = True
for L in data["descent"]["check_lattices"]:
    q, s, t = int(L["q"]), int(L["s"]), int(L["root"])
    cols = [[int(v) for v in c] for c in L["basis_columns"]]
    cvec = [pow(s, i, q) for i in range(6)] + [(-t * pow(s, i, q)) % q for i in range(6)]
    for c in cols:
        lat_ok = lat_ok and sum(ci * vi for ci, vi in zip(cvec, c)) % q == 0
    lat_ok = lat_ok and abs(T.det_int([[cols[j][i] for j in range(12)] for i in range(12)])) == q
    hs = sum(hc * pow(s, i, q) for i, hc in enumerate(T.H)) % q
    lat_ok = lat_ok and hs == 0
checks["special_q_lattices_verified"] = lat_ok and len(data["descent"]["check_lattices"])

# 5. stored boot vectors: a + tau b = 0 in F_p[y]/h and the g-side norm is as stored
boot_ok = True
for L, row in zip(data["boot"]["check_lattices"], data["boot"]["rows"]):
    tau = [int(v) for v in L["tau"]]
    a = [int(v) for v in L["a"]]
    b = [int(v) for v in L["b"]]
    z = T.hmul(tau, [v % p for v in b], p)
    boot_ok = boot_ok and all((ai + zi) % p == 0 for ai, zi in zip(a, z))
    boot_ok = boot_ok and abs(T.log2_int(T.norm_g(a, b)) - row[0] / 10) < 0.06
checks["boot_vectors_verified"] = boot_ok and len(data["boot"]["check_lattices"])

# 6. Dickman table spot values and the complete ledger
D = T.Dickman()
dick = {str(k): D(k) for k in (2, 5, 8, 10, 12)}
result = T.model(data)
try:
    claim = json.loads((base / "submission.json").read_text())["claims"][0]["score"]
    claimed = int(claim["numerator"]) / int(claim["denominator"])
except (OSError, KeyError, ValueError):
    claimed = data["config"]["declared_score"]
checks["declared_score"] = claimed
checks["declared_score_covers_computed_bound"] = claimed >= result["log2_expected_total_units"]
print(json.dumps({
    "seed": seed, "checks": checks, "fresh_cells": fresh_out, "dickman": dick, "model": result,
    "scope": "Exact arithmetic at the production parameters, replay of stored samples and ledger arithmetic. "
             "No production attack, no proof of heuristics H1-H7, no verifier acceptance is claimed.",
}))
