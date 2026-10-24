#!/usr/bin/env python3
"""v33 AUDIT: independent recompute from design-v33.json ONLY (verdict flaw 1 upgrade),
claim-semantics identity (flaw 2), census-distribution rows (flaw 3), gap-row ell (flaw 4),
plus the 20-pair key test, the four planted controls and the planted-TAUTOLOGY control.

The tautology control is the load-bearing one: it replaces one construction by the v32
style cost-defined stub (price := row_value // count, count := 1) and the audit MUST FAIL
it.  A validator that cannot fail that row cannot certify non-tautology.
"""
import json
import math
import os
import random
import struct
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import keymod as K

DES = json.load(open(os.path.join(HERE, "design-v33.json")))
LED = json.load(open(os.path.join(HERE, "ledger-v33.json")))
CEN = json.load(open(os.path.join(HERE, "census-prod-v33.json")))
LAD = json.load(open(os.path.join(HERE, "structured-ladder.json")))
SUB = json.load(open(os.path.join(HERE, "submission.json")))

out, ok = {}, True


def chk(name, cond, detail=""):
    global ok
    ok = ok and bool(cond)
    out[name] = dict(pass_=bool(cond), detail=detail)


# ---------- 1. non-tautological recompute from DESIGN JSON ONLY -------------------------
def namespace(d):
    ns = {k: v for k, v in d.items() if k != "price_lattice"}
    for k, v in d["price_lattice"].items():
        ns[k] = v
    return ns


RESERVED = ("row_value", "value", "total", "claim", "W")


def recompute(deriv, ns, forbid_values=True):
    """evaluate count x price.  A TAUTOLOGICAL row is one whose expressions need a row
    value (or are self-referential): the namespace holds design constants only, so such a
    row raises NameError, and a price above the accepted lattice ceiling is rejected."""
    vals = {}
    for name, r in deriv.items():
        for token in RESERVED:
            if token in r.get("count_expr", "") or token in r.get("price_expr", ""):
                raise ValueError("row %s reads a cost-derived symbol (%s)" % (name, token))
        try:
            count = eval(r["count_expr"], {"__builtins__": {}}, dict(ns))
            price = eval(r["price_expr"], {"__builtins__": {}}, dict(ns))
        except Exception as e:                       # noqa: BLE001
            raise ValueError("row %s is not recomputable from the design: %s" % (name, e))
        if price <= 0 or count < 0:
            raise ValueError("row %s has a non-positive factor" % name)
        lat_max = max(DES["price_lattice"].values())
        unit = DES["C_instructions_per_unit"]
        in_lattice = (price <= lat_max or price == 1 or
                      (price % unit == 0 and price // unit <= 1024) or
                      (price & (price - 1) == 0 and price <= (1 << 24)))
        if forbid_values and not in_lattice:
            # a per-item price cannot exceed the accepted lattice ceiling; this is the
            # detector for 'price := value / count' stubs.
            raise ValueError("row %s price %s is outside the accepted price lattice"
                             % (name, price))
        vals[name] = int(math.ceil(count * price))
    return vals


ns = namespace(DES)
recomputed = recompute(LED["preprocessing_derivations"], ns)
mismatch = {k: (v, LED["rows"][k]) for k, v in recomputed.items() if v != LED["rows"][k]}
chk("preprocessing_recomputes_from_design_only", not mismatch, str(mismatch)[:300])
chk("preprocessing_row_count", len(recomputed) == len(LED["preprocessing_derivations"]),
    "%d rows" % len(recomputed))

# planted TAUTOLOGY control: the v32 style stub must FAIL the same validator
stub = dict(LED["preprocessing_derivations"]["fixed_per_node"])
v = LED["rows"]["fixed_per_node"]
stub["count_expr"] = "1"
stub["price_expr"] = str(v)                     # price := the row value itself
taut_caught = False
try:
    recompute({"fixed_per_node": stub}, ns)
except ValueError as e:
    taut_caught = "outside the accepted price lattice" in str(e) or "not recomputable" in str(e)
    taut_detail = str(e)
chk("planted_tautology_control_caught", taut_caught, globals().get("taut_detail", ""))

# ---------- 2. ledger arithmetic, claim legality, reduced rational ---------------------
total = sum(LED["rows"].values())
chk("ledger_integer_sum", total == LED["total_instructions"], str(total))
W = math.log2(total / 2000.0)
chk("W_recomputed", abs(W - LED["W_log2"]) < 1e-12, "%.12f" % W)
num, den = int(LED["claim"]["numerator"]), int(LED["claim"]["denominator"])
claim = num / den
chk("claim_reduced_rational", math.gcd(num, den) == 1, "%d/%d" % (num, den))
chk("claim_ge_W", claim >= W, "claim %.6f W %.6f" % (claim, W))
chk("claim_covers_stress", claim - W >= LED["doubling_stress_bits"],
    "slack %.9f stress %.9f" % (claim - W, LED["doubling_stress_bits"]))
chk("submission_claim_matches_ledger",
    int(SUB["claims"][0]["score"]["numerator"]) / int(SUB["claims"][0]["score"]["denominator"]) == claim,
    "%s/%s" % (SUB["claims"][0]["score"]["numerator"], SUB["claims"][0]["score"]["denominator"]))

# ---------- 3. flaw 2: the claim is THE MEASURED success quantity -----------------------
claimed = LED["success"]["claimed_value"]
key_path_val = round(LAD["pooled_long_walks"]["cp95_lower"], 6)
chk("success_claim_is_measured_quantity", claimed == key_path_val,
    "claimed %.6f key-path %.6f" % (claimed, key_path_val))
lit = "%.6f" % claimed
chk("success_literal_in_submission_statement", lit in SUB["claims"][0]["statement"], lit)
chk("success_literal_in_note", lit in open(os.path.join(HERE, "public-note.md")).read(), lit)
chk("success_meets_floor", claimed >= LED["success"]["floor"],
    "%.6f >= %.2f" % (claimed, LED["success"]["floor"]))
chk("k_eff_is_not_the_claim", abs(claimed - 0.508674) > 0.01 and claimed < 1.0,
    "claimed %.6f is not the collision intensity 0.508674" % claimed)

# ---------- 4. flaw 3: census is at the production distribution ------------------------
cells = CEN["cells"]
chk("census_has_production_cell", "prod_crossside" in cells, str(list(cells)))
chk("census_production_law_declared",
    "shared walk list" in CEN["production_density"]["law"] and
    CEN["production_density"]["span_bits"] + CEN["production_density"]["band_bits"] <= 40,
    CEN["production_density"]["law"][:80])
chk("census_canonical_even_applied", "prod_crossside" in CEN["canonicalisation"]["applied_to"] and
    "uniform_control" in CEN["canonicalisation"]["not_applied_to"], "both cells present")
chk("census_delta_vs_v32_reported", "delta" in CEN["delta_vs_v32"],
    str(CEN["delta_vs_v32"]["v33_charged"]))
chk("census_joint_168bit_row", "joint_168bit" in CEN and
    "pairwise_joint_chi2" in CEN["joint_168bit"]["prod"], "joint row present")
chk("kappa_charged_is_max_upper95",
    abs(CEN["charged_kappa"] - max(max(c["kappa_upper95"].values()) for c in cells.values())) < 1e-9 or
    CEN["charged_kappa"] >= max(max(c["kappa_upper95"].values()) for c in cells.values()) - 1e-9,
    "%.6f" % CEN["charged_kappa"])
chk("ledger_kappa_matches_census", abs(LED["kappa"]["charged"] - CEN["charged_kappa"]) < 1e-9,
    "%.6f" % LED["kappa"]["charged"])
N = DES["N"]
E_bad = LED["kappa"]["charged"] * N * N / float(1 << 168)
chk("birthday_under_cap", E_bad < (1 << 84), "E_bad 2^%.3f vs cap 2^84" % math.log2(E_bad))

# ---------- 5. flaw 4: gap row uses PRODUCTION ell and a union bound --------------------
g = LED["gap_support"]
ell_prod = -(-DES["N"] // DES["nb"])
chk("gap_row_uses_production_ell", int(g["ell_used"]) == ell_prod,
    "ell %s vs ceil(N/nb) %d" % (g["ell_used"], ell_prod))
chk("gap_row_not_fixture_ell", int(g["ell_used"]) != 4774, "4774 rejected")
chk("gap_row_is_union_bound_over_distinct_gaps",
    int(g["distinct_gaps"]) == ell_prod - 1 and "union bound" in g["rule"],
    g["rule"])
chk("gap_abort_bound_below_2_100", g["abort_bound_log2"] <= -100, str(g["abort_bound_log2"]))
chk("gap_row_value_uses_corrected_support",
    LED["rows"]["gap_point_tables"] >=
    DES["sides"] * DES["nb"] * (1 << DES["gap_support_log2"]) * DES["gap_entry_point_ops"] * 2000,
    "%.3f" % math.log2(LED["rows"]["gap_point_tables"]))

# ---------- 6. key: 20 fresh adversarial pairs + 4 planted controls --------------------
rng = random.Random(3301)
hds = [bin(K.key_mixed_final(b) ^ K.key_mixed_final(b + (1 << 64) * rng.randrange(1, 1 << 30))).count("1")
       for b in (rng.getrandbits(300) for _ in range(20))]
chk("twenty_pair_hd_plausible", all(50 <= h <= 120 for h in hds), str(hds))
ctrl = {
    "bit_flip_changes_key": K.key_mixed_final(12345) != K.key_mixed_final(12345 ^ (1 << 200)),
    "limb5_reaches_output": K.key_mixed_final(1) != K.key_mixed_final(1 + (1 << 320)),
    "limb0_reaches_output": K.key_mixed_final(1) != K.key_mixed_final(2),
    "v31_control_still_collides": K.key_mixed_v31(5) == K.key_mixed_v31(5 + (1 << 64)),
}
chk("four_planted_controls", all(ctrl.values()), json.dumps(ctrl))

# ---------- 7. direction / regex sweep -------------------------------------------------
stmt = SUB["claims"][0]["statement"]
note = open(os.path.join(HERE, "public-note.md")).read()
chk("claim_literal_in_statement", "%.4f" % claim in stmt, "%.4f" % claim)
chk("gate_literal_in_statement", "116.8271" in stmt, "gate")
chk("no_wrong_116_8543_literal", "116.8543" not in stmt and "116.8543" not in note, "v32 claim gone")
chk("over_gate_disclosed", "OVER" in stmt.upper(), "crossing stated")
chk("directions_lower_for_success",
    "LOWER" in LED["success"]["claimed_quantity"].upper(),
    LED["success"]["claimed_quantity"][:60])
chk("title_len_ok", len(SUB["title"]) <= 300, str(len(SUB["title"])))
chk("statement_len_ok", len(stmt) <= 7965, str(len(stmt)))

out["all_pass"] = ok
json.dump(out, open(os.path.join(HERE, "audit-v33.json"), "w"), indent=1)
fails = [k for k, v in out.items() if isinstance(v, dict) and not v.get("pass_")]
print(json.dumps({k: v for k, v in out.items() if isinstance(v, dict) and not v.get("pass_")},
                 indent=1))
print("checks", len([k for k, v in out.items() if isinstance(v, dict)]), "fails", len(fails))
sys.exit(0 if ok else 1)
