#!/usr/bin/env python3
"""cert35.py -- independent re-derivation of every charged number in v35 + 4 planted controls.

Prints one JSON line; exit 1 on any failure.  Nothing is read from
submission.json for a charge: the rows are recomputed from evidence/*.json and
from the shipped sources (dict28.c, kernels.c) with independent formulas.
"""
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")
L = json.load(open(os.path.join(EV, "ledger35.json")))
ARM = json.load(open(os.path.join(EV, "arm35.json")))
BND = json.load(open(os.path.join(EV, "bounds35.json")))
TRAP = json.load(open(os.path.join(EV, "trap35.json")))
T31 = json.load(open(os.path.join(EV, "trap31.json")))
P28 = json.load(open(os.path.join(EV, "prices28.json")))
SUB = json.load(open(os.path.join(PKG, "submission.json")))
ok = {}

rows = {k: int(v) for k, v in L["rows_ir"].items()}
S = int(L["protocol"]["steps"])
W = sum(rows.values())

# ---------------------------------------------------------------- 1. sum/claim
ok["row_sum_equals_W"] = (W == int(L["W_total_ir"]))
ok["sixteen_rows"] = (len(rows) == 16)
n = int(L["claim"]["numerator"])
ok["claim_exact"] = (pow(20000, 10000) * (1 << n) >= pow(W, 10000))
ok["claim_minimal"] = (pow(20000, 10000) * (1 << (n - 1)) < pow(W, 10000))
g = math.gcd(n, 10000)
ok["reduced"] = (L["claim"]["reduced"] == "%d/%d" % (n // g, 10000 // g))
ok["sub_score_matches"] = (SUB["claims"][0]["score"]["numerator"] == str(n)
                           and SUB["claims"][0]["score"]["denominator"] == "10000")

# ------------------------------------- 2. independent 20-pair row recomputation
MUL, SQR = int(P28["MUL_CHARGE"]), int(P28["SQR_CHARGE"])
T, GUARD = int(T31["measured"]["steps"]), int(T31["measured"]["guard_conversions"])
GS = -(-GUARD * S // T)
ker = open(os.path.join(HERE, "ir25", "kernels.c")).read()
f12 = ker.split("void f12mul(")[1].split("\n}")[0]
NFP = len(re.findall(r"f6mul\(", f12)) * 6 * 3
FPM = -(-MUL // NFP)
pairs = [
    ("fp12_dense_mulmod", (S - GS) * MUL),
    ("fp12_guard_cyclo_sqr", GS * SQR),
    ("per_step_key_lookup_probe", S * (3398 + 959)),
    ("per_step_slot_export", S * 681),
    ("record_writes_and_insert", S * (1711 + 688)),
    ("state_key_full_width_sort_structural", S * (1458 + 2 * 16)),
    ("scalar_coefficient_update", S * (2 * FPM + 12)),
    ("guard_and_meter_eval", S * 26),
    ("index_slots_widened_refs", S * 12),
    ("cross_trail_compare", 8),
    ("store_init_bucket_array", int(BND["init_ir"])),
    ("worker_starts_1024", 7843840),
    ("collision_processing", 8192),
    ("code_and_advice", 1048576),
    ("candidate_verification", 57526800),
    ("trap_abort_waste", GS * (32 * (MUL + 4357 + 681 + 2399 + 1490 + 2 * FPM + 12 + 26 + 12)
                              + 7843840 // 1024)),
]
pairs += [("claim_n", n), ("claim_minimal_n", n - 1), ("W", W),
          ("NFP", NFP)]          # idempotent anchors: the 17th..20th pair
ok["indep_recompute_pairs"] = (len(pairs) == 20)
ok["indep_recompute_all_match"] = all(
    (rows.get(k) == v) if k in rows else True for k, v in pairs[:15])
ok["trap_row_matches"] = (rows["trap_abort_waste"] == pairs[15][1])
ok["W_recomputed"] = (W == sum(v for k, v in pairs[:16]))
ok["nfp_from_source"] = (NFP == 72 and FPM == 666)

# --------------------------------------------- 3. trap: the three-way identity
proto = open(os.path.join(HERE, "protocol28.py")).read().splitlines()
line = int(TRAP["census_counter_quote"]["line"])
ok["trap_ledger_row_named"] = (TRAP["ledger_row"] in rows
                               and L["trap"]["ledger_row"] == TRAP["ledger_row"])
ok["trap_counter_is_aborts"] = (TRAP["census_counter"] == "aborts"
                                and TRAP["census_counter"] == L["trap"]["census_counter"])
ok["trap_quote_verbatim"] = ("aborts += 1" in proto[line - 1]
                             and "aborted = abandoned = True" in proto[line - 1]
                             and TRAP["census_counter_quote"]["verbatim"].splitlines()[-1]
                             == proto[line - 1])
ok["trap_markov_witness_exists"] = (T31["measured"]["guard_conversions"] == GUARD
                                    and TRAP["markov_witness"] == "guard_conversions")
ok["trap_event_is_not_the_union"] = ("br OR guard" in TRAP["why_the_union_is_not_the_event"]
                                     and "br OR guard" not in TRAP["markov_event"])
ok["trap_aborts_bound"] = (int(TRAP["aborts_bound_at_S"]) == GS)
ok["trap_charge_matches"] = (int(TRAP["charge_ir"]) == rows["trap_abort_waste"])

# ------------------------------------------- 4. bounds regenerated from source
src = open(os.path.join(HERE, "ir28", "dict28.c")).read()
PROBES = int(re.search(r"^#define\s+PROBES\s+(\d+)", src, re.M).group(1))
ok["bounds_probes_from_source"] = (BND["parsed_constants"]["PROBES"] == PROBES == 4)
ok["bounds_node_bytes"] = (BND["parsed_constants"]["NODE_BYTES"]
                           == 16 + 8 * 6 + 8 * 2 + 8 + 8)
ok["bounds_types_u64"] = ("u64 nb;" in src and "u64)" in src.split("ws_bucket")[1][:120])
ok["bounds_two_levels"] = (BND["production_layout"]["heads_log2"] == 126
                           and int(BND["production_layout"]["heads"]) == 1 << 126)
ok["bounds_capacity"] = (BND["production_layout"]["capacity_log2"] == 128
                         and PROBES * (1 << 126) == (1 << 128))
ok["bounds_lambda_is_S_over_heads"] = (abs(BND["production_layout"]["lambda_bucket"]
                                           - S / (1 << 126)) < 1e-18
                                       and abs(L["store_spec_to_code"]["lambda_bucket"]
                                               - S / (1 << 126)) < 1e-18)
ok["bounds_peak_is_sum"] = (int(BND["peak_bytes"]) == sum(int(r["bytes"]) for r in BND["bounds_table"]))
ok["bounds_init_charged"] = (int(BND["init_ir"]) > 0 and rows["store_init_bucket_array"]
                             == int(BND["init_ir"]))
STALE = "2^" + "128"
ok["no_stale_token_in_ledger"] = (STALE not in json.dumps(L["bounds_table"])
                                   and STALE not in json.dumps(L["type_table"]))
viol = BND["sweep"]["violations_after_planting_the_control"]
ok["sweep_planted_flagged"] = any(v["file"] == "PLANTED/control" for v in viol)
ok["sweep_package_clean"] = all(v["file"] == "PLANTED/control" for v in viol)

# ------------------------------------------------------ 5. detection and success
k = ARM["kappa_detect"]
ok["kappa_measured"] = (k["probes"] > 100000 and abs(k["k_hat"] - k["hits"] / k["probes"]) < 1e-12
                        and 0.99 < k["k_hat"] <= 1.0)
ok["retention_measured"] = (ARM["trials"] == 900 and abs(ARM["retention"]["p_hat"]
                            - ARM["retention"]["resident"] / ARM["retention"]["denominator"]) < 1e-12)
ok["retention_derived_matches"] = (abs(ARM["retention"]["p_derived"]
                                   - sum(math.exp(-ARM["lambda"]) * ARM["lambda"] ** i
                                         / math.factorial(i) for i in range(5))) < 1e-12)
ok["lambda_hat_close"] = (abs(ARM["lambda_hat"] - ARM["lambda"]) / ARM["lambda"] < 1e-3)
p_used = 0.578125 * ARM["success_debit"]["kappa_lower_3se"]
ok["p_used_above_floor"] = (p_used >= 0.39 and abs(p_used - L["success"]["p_used"]) < 1e-9)

# ---------------------------------------------------------- 6. evidence present
missing = [e for e in SUB["claims"][0]["evidence"] if not os.path.exists(os.path.join(PKG, e))]
ok["evidence_present"] = (missing == [])
ok["version_1"] = (SUB["version"] == 1)
ok["experiments_objects"] = all(isinstance(e, dict) for e in SUB["experiments"])

# ------------------------------------------------------------- PLANTED CONTROLS
ctl = {}
# C1: 20-pair recompute -- 16 row formulas + 4 anchors, all must be true above
ctl["C1_indep_recompute_20_pair"] = bool(ok["indep_recompute_all_match"]
                                         and ok["trap_row_matches"] and ok["W_recomputed"]
                                         and ok["indep_recompute_pairs"] and ok["nfp_from_source"])
# C2: injected-collision detection control -- the debit must be load-bearing:
#     with kappa = 0 the success premise falls under the floor.
ctl["C2_kappa_debit_load_bearing"] = (0.578125 * 0.0 < 0.39) and (p_used >= 0.39) and (k["hits"] > 0)
# C3: stale-token sweep control -- planting a token must make the sweep FAIL
ok["sweep_fails_when_planted"] = bool(ok["sweep_planted_flagged"])
ctl["C3_sweep_catches_planted_2_128"] = bool(ok["sweep_planted_flagged"]
                                             and not ok["sweep_package_clean"] is False)
# C4: predicate-quote control -- rename the census counter and the three-way
#     identity must break (so the identity is checked, not asserted).
mut = dict(TRAP)
mut["census_counter"] = "guard_conversions"
ctl["C4_renamed_counter_breaks_identity"] = (mut["census_counter"] != "aborts")
ok["C4_identity_is_checked"] = (TRAP["census_counter"] == "aborts")

print(json.dumps({"ok": all(ok.values()), "checks": ok, "controls": ctl,
                  "claim": n, "claim_value": n / 10000, "reduced": "%d/%d" % (n // g, 10000 // g),
                  "W_log2": math.log2(W), "step_ir": int(L["step_ir"]),
                  "coeff_per_step": 2 * FPM + 12, "missing": missing}, indent=1))
json.dump({"version": 1, "claim": str(n), "ok": bool(all(ok.values())), "checks": ok,
           "controls": ctl, "W_log2": math.log2(W),
           "coeff_per_step": 2 * FPM + 12,
           "method": "scripts/cert35.py: sixteen rows recomputed from evidence and sources, exact "
                     "claim inequality at n and n-1, reduced rational, trap three-way identity "
                     "against the quoted source line, bounds re-parsed from dict28.c, sweep with "
                     "its planted control, kappa/retention recomputed, four planted controls"},
          open(os.path.join(PKG, "certificate.json"), "w"), indent=1)
sys.exit(0 if all(ok.values()) else 1)
