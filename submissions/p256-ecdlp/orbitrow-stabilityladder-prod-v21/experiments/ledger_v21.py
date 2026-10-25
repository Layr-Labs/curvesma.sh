"""v20 ledger.  Claim is computed once here; make_submission_v20.py copies it.

WHAT CHANGED AND WHY (verdict d32b4f73, all three flaws):

flaw 1, success transfer.  The carrier is no longer the W=1 production cell (every
collision in it is inside one walker's own trail, which is the regime the review
showed production is NOT in: production per-trail collision load is 2^-40.66 and
the success comes from collisions between distinct walkers).  The carrier is now
`X_R256_W64_p16`, the same production selector at 64 walkers where only 1.6% of
collisions are within-walker (production: ~0%) and success is counted from
CROSS-WALKER relations only, at n = 2000.  The W=1 production cell is retained as
the big-R replication, not as the bound.

flaw 2, return/restart.  Charged from returns_v21.json: per TRAIL, exact-state
tracking, a structural argument that classical rho-cycle deaths cannot occur on
an advancing-counter selector, and the identity-restart rate at its CP-95 upper
bound -- not divided by the walker count.

flaw 3, store.  Charged from store_v21.json: 16 B register-pair transfers, page
init/flush/fill at 256 transfers per 4096 B page, fanout 2^21 x 6 levels, and the
radix directories charged and counted in the peak.
"""
import json, math, os
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "evidence")

LOG2_V = 127.67
WALKERS = 1 << 20
STEPS = 2 ** (LOG2_V - 20)


def cp95_lower(k, n, alpha=0.025):
    if k == 0:
        return 0.0
    lg = math.lgamma

    def tail(p):
        lp, lq = math.log(p), math.log1p(-p)
        return sum(math.exp(lg(n + 1) - lg(i + 1) - lg(n - i + 1) + i * lp + (n - i) * lq)
                   for i in range(k, n + 1))

    lo, hi = 0.0, 1.0
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if tail(mid) < alpha:
            lo = mid
        else:
            hi = mid
    return lo


def main():
    price = json.load(open(os.path.join(OUT, "insn_price.json")))
    fe = price["prices_ir"]["fe_mul"]
    fa = price["prices_ir"]["fe_add"]
    fs = price["prices_ir"]["fe_sub"]
    store = json.load(open(os.path.join(OUT, "store_v21.json")))
    ret = json.load(open(os.path.join(OUT, "returns_v21.json")))
    drift = json.load(open(os.path.join(OUT, "drift_v20.json")))
    carrier = json.load(open(os.path.join(OUT, "ladder20_X_R256_W64_p16.json")))
    rep = json.load(open(os.path.join(OUT, "ladder_R2p20_W1_16p.json")))   # v19 production cell

    sort_ir = store["sort_total_ir_per_record"]
    rows = {
        "affine_add": 3 * fe + 4 * fa + 2 * fs,
        "batch_invert_share": 3 * fe,
        "batch_chain_share": (385 * fe) / (1 << 20),
        "table_read_index": 10,
        "canon_sign_rule": 20,
        "record_append_48B": 12,
        "global_sort_and_witness": sort_ir,
    }
    per_visit = sum(rows.values())
    units = per_visit / 1000.0
    log2_rows = LOG2_V + math.log2(units)

    # ---- success: the CROSS-WALKER carrier cell carries the bound -------------
    assert carrier["walkers"] == 64 and carrier["periods"] == 16
    n = carrier["trials"]
    wins_cross = carrier["wins_cross"]
    wins_all = carrier["wins"]
    p_lo = cp95_lower(wins_cross, n)
    p_hi_all = wins_all / n

    # drift gate: if ANY coverage axis drifts down, the carrier must be replaced
    # by that axis's largest-coverage cell (the honest drift-consistent value).
    drift_down = [a for a in drift["axes"] if a["verdict"] == "DRIFT_DOWN"]
    drift_cells = []
    for a in drift_down:
        c = a["cells"][-1]
        drift_cells.append({"axis": a["axis"], "label": c["label"], "p_hat": c["p_hat"],
                            "cp95_lower": c["cp95_lower"], "trials": c["trials"]})
    p_used = min([p_lo] + [c["cp95_lower"] for c in drift_cells])
    p_point = min([wins_cross / n] + [c["p_hat"] for c in drift_cells])

    # ---- frequency x utility decomposition (the transfer argument) -----------
    # Frequency is ANALYTIC and conceded by the review: the production attack
    # writes N = 2^127.67 records into a store with r = 2^255 sign classes, so
    # E[collisions] = N(N-1)/(2r).  Utility is the only MEASURED factor:
    # P(recovery | collision), which every ladder cell measures as 1.0.
    N_prod = 2 ** LOG2_V
    classes_prod = 2 ** 255
    lam_prod = N_prod * (N_prod - 1) / (2.0 * classes_prod)
    n_coll_carrier = carrier["collisions"]
    assert carrier["useful_relations"] == n_coll_carrier
    q_util_lo = 0.025 ** (1.0 / n_coll_carrier)      # CP-95 lower, k = n = collisions
    success_analytic = 1.0 - math.exp(-lam_prod * q_util_lo)
    p_used = min(p_used, success_analytic)

    # attempts needed so that P(at least one success) >= 0.39
    if p_point > 0:
        attempts = max(1.0, math.log(1 - 0.39) / math.log(1 - p_point))
    else:
        attempts = float("inf")
    assert p_used >= 0.39, ("carrier CP-95 lower %.6f misses the 0.39 floor" % p_used)
    log2_attempts = math.log2(attempts) if attempts > 1 else 0.0

    lam_meas = carrier["collisions"] / n
    useful_frac = carrier["useful_relations"] / max(1, carrier["collisions"])
    assert useful_frac == 1.0, useful_frac

    full_step = 388 * fe + 4 * fa + 2 * fs + 30
    bulk_step = rows["affine_add"] + rows["batch_invert_share"] + 30
    infl = full_step / bulk_step
    trails_per_challenge = 2.0 * (1.0 + lam_meas)
    replay_rel = trails_per_challenge * infl / WALKERS

    rel_restart = ret["production_charge"]["charged_relative_loss"]
    log2_work = log2_rows + math.log2(1.0 + rel_restart) + math.log2(1.0 + replay_rel)
    table_ir = 2 * (1 << 20) * 256 * 3 * fe
    log2_work += math.log2(1.0 + table_ir / 2 ** log2_work)
    log2_work += log2_attempts

    claim = math.ceil(log2_work * 1e4) / 1e4
    assert claim >= log2_work
    num = int(round(claim * 10000))
    g = gcd(num, 10000)
    rational = {"numerator": num // g, "denominator": 10000 // g}

    ev = {
        "rows_ir_per_visit": {k: round(v, 6) for k, v in rows.items()},
        "instructions_per_visit": per_visit,
        "units_per_visit": units,
        "sort": {
            "per_pass_ir_per_record": store["per_pass_total_ir"],
            "passes": store["passes"],
            "store_io_ir_per_record": store["store_io_ir_per_record_per_pass"],
            "readback_ir_per_record": store["readback_ir_per_record_per_pass"],
            "per_pass_components": store["per_pass_ir_per_record"],
            "witness_ir": store["witness_ir_per_record"],
            "directory_construction_ir": store["charges_ir_per_record"]["directory_construction_amortised"],
            "total_sort_ir_per_record": sort_ir,
            "peaks": store["peaks"],
            "basis": "PROTECTED 64-bit GPR accounting (evidence/store_v21.json): 16 B per "
                     "register-pair transfer, 3 transfers per 48 B record each way, 256 per "
                     "4096 B page init/flush/fill, fanout 2^21 x 6 levels over the page number, "
                     "directories charged and included in the peak",
        },
        "success": {
            "floor": 0.39,
            "carrier": {
                "label": carrier["label"], "R": carrier["R"], "walkers": carrier["walkers"],
                "periods": carrier["periods"], "trials": n,
                "wins_cross_walker": wins_cross, "rate_cross_walker": wins_cross / n,
                "cp95_lower_cross_walker": p_lo,
                "wins_any": wins_all, "rate_any": p_hi_all,
                "within_walker_collisions": carrier["within_collisions"],
                "cross_walker_collisions": carrier["cross_collisions"],
                "cross_share": carrier["cross_collisions"] / carrier["collisions"],
                "within_walker_load": carrier["within_walker_load"],
                "why_carrier": "production success comes from collisions between DISTINCT "
                               "walkers (per-trail load 2^-40.66); at W=64 only 1.6% of "
                               "collisions are within-walker, and the bound counts ONLY "
                               "cross-walker relations, so it is the mix-matched statistic",
            },
            "big_r_replication_not_the_bound": {
                "label": rep["label"], "R": rep["R"], "walkers": rep["walkers"],
                "trials": rep["trials"], "wins": rep["verified_successes"],
                "rate": rep["verified_successes"] / rep["trials"],
                "cp95_lower": rep["cp95_lower"],
                "note": "same walk at R=2^20, W=1, 16 periods; retained as the large-R "
                        "replication of the same selector",
            },
            "p_used_cp95_lower": p_used,
            "p_used_point": p_point,
            "drift_down_axes": drift_down and [a["axis"] for a in drift_down] or [],
            "attempts_charged": attempts,
            "log2_attempts_charged": log2_attempts,
            "transfer": {
                "measured_axes": [{"axis": a["axis"], "verdict": a["verdict"],
                                   "slope_per_log2_axis": a["slope_per_log2_axis"]}
                                  for a in drift["axes"]],
                "statement": "success per challenge is measured on the production selector, "
                             "cross-walker only, and is FLAT in coverage (16..1024 periods), "
                             "in walker count (1..64) and in R (256..65536) wherever the exact "
                             "CP intervals overlap; the remaining step to 2^87.67 periods is "
                             "DISCLOSED AS A CONJECTURE (a rho-style heuristic: the per-period "
                             "row hazard has no measured dependence on the period count), not "
                             "as a measurement, and the bound charged is the carrier cell's own "
                             "CP-95 lower bound, not an extrapolated point",
                "not_claimed": "that a 37-53 bit fixture group is the P-256 subgroup; the "
                               "fixture is a stand-in at matched load lambda = 0.634",
            },
        },
        "policy": {
            "record_fields": "32B canonical key || 20B identity pack (20-bit trail || 108-bit step)",
            "rule": "every hit replays both trails from their stored identities; the replay "
                    "decides triviality and, if nontrivial, recovers x=-da/db",
            "measured_collisions_per_challenge": lam_meas,
            "measured_useful_fraction": useful_frac,
            "trail_pairs_charged_per_challenge": trails_per_challenge,
            "replay_inflation_vs_bulk": infl,
            "replay_relative_loss": replay_rel,
        },
        "returns": {
            "relative_loss": rel_restart,
            "correction": ret["correction"],
            "high_load_cell": {k: ret["high_load_cell"][k] for k in
                               ("q", "R", "steps_per_walker", "trials", "expected_repeats_analytic",
                                "measured", "fatal_fraction_cp95_upper", "detects_repeats")},
            "production_charge": ret["production_charge"],
            "source": "returns_v21.json: per-TRAIL exact tracking on the declared walk, "
                      "38,400 trails, 0 exact (S,t) returns, 0 identity restarts, positive "
                      "control detects 1.04e6 S returns at q=1031",
            "structural_argument": ret["correction"]["fixed_here"],
            "measured": ret["high_load_cell"]["measured"],
            "fatal_fraction_cp95_upper": ret["high_load_cell"]["fatal_fraction_cp95_upper"],
            "non_vacuity": "expected repeats computed analytically before the run; the measured "
                           "count is reported against it",
            "v19_error": ret["v20_error"],
        },
        "work": {"log2_visits": LOG2_V, "walkers": WALKERS, "steps_per_walker": STEPS,
                 "log2_rows": log2_rows, "log2_work_final": log2_work},
        "bounds": {"work_log2": claim, "claim": claim,
                   "storage_log2_peak": store["peaks"]["log2_peak_bytes"]},
        "record": {"bytes": 48,
                   "layout": "32B canonical key || 16B identity (20-bit trail || 108-bit step)"},
        "claim_time_log2": claim,
        "claim_reduced_rational": rational,
        "decomposition": {
            "why": "the review's arithmetic separates collision FREQUENCY (analytic) from "
                   "collision UTILITY (measured); production success is a whole-attack event "
                   "over 2^20 walkers' records, not a per-walker hit",
            "frequency": {"records": "2^%.2f" % LOG2_V, "sign_classes": "2^255",
                          "expected_collisions_per_attack": lam_prod,
                          "formula": "E[X] = N(N-1)/(2r)",
                          "conceded_by_review": "the verdict states the same birthday arithmetic "
                                                "((2^107.67)^2/r per trail, and the 0.634 "
                                                "per-challenge load)"},
            "utility": {"definition": "q = P(recovery | collision)",
                        "carrier_collisions": n_coll_carrier,
                        "carrier_recoveries": carrier["useful_relations"],
                        "q_hat": 1.0, "q_cp95_lower": q_util_lo,
                        "all_cells": "every one of the 20 ladder cells has "
                                     "useful_relations == collisions (control (f)), so q = 1.0 "
                                     "is measured at 16/64/256/1024 periods, W 1..64 and "
                                     "R 256..65536"},
            "success_per_attack": {"analytic_from_decomposition": success_analytic,
                                   "measured_cross_walker_cp95_lower": p_lo,
                                   "charged": p_used},
            "assumptions_disclosed": "X ~ Poisson and collisions mutually independent; the "
                                     "stability ladder supports the absence of long-range "
                                     "correlation (flat in 16..1024 periods, W 1..64, "
                                     "R 256..65536) but does not prove independence at 2^87.67 "
                                     "periods. This is a rho-heuristic CONJECTURE, and the bound "
                                     "charged is the smaller of the analytic value and the "
                                     "measured carrier CP-95 lower bound.",
        },
        "delta_vs_v19": {
            "v19_claim": 128.7921,
            "v19_carrier": "R2p20_W1_16p, 95/200 = 0.475, CP-lo 0.404120, W=1 (all collisions "
                           "within one trail)",
            "v20_carrier": "%s cross-walker, %d/%d = %.4f, CP-lo %.6f"
                           % (carrier["label"], wins_cross, n, wins_cross / n, p_lo),
            "price_delta": "store price %.4f -> %.4f Ir/record/pass (16 B rule + page I/O)"
                           % (65.417469, store["per_pass_total_ir"]),
            "return_delta": "5.845e-09 (per-walker/2^20) -> %.3e (per-trail CP-95 upper)"
                            % rel_restart,
        },
    }
    json.dump(ev, open(os.path.join(OUT, "ledger_v21.json"), "w"), indent=1, sort_keys=True)
    json.dump({"claim": claim, "rational": rational, "log2_work": log2_work},
              open(os.path.join(OUT, "claim_v21.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({"per_visit": round(per_visit, 3), "sort_ir": round(sort_ir, 3),
                      "log2_rows": log2_rows, "log2_work": log2_work, "claim": claim,
                      "rational": rational,
                      "carrier": {"wins_cross": wins_cross, "n": n, "cp95_lower": p_lo},
                      "attempts": attempts, "restart_rel": rel_restart,
                      "replay_rel": replay_rel}, indent=1))


if __name__ == "__main__":
    main()
