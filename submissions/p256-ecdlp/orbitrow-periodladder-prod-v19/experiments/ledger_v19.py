"""v19 ledger.

Changes from v18 (all four verdict-v18 flaw classes):
  * the sort row now prices the ANALYTIC production store (store_v19.json):
    36 Ir carried-scatter body + external translation / page I/O / read-back,
    with no bounded fixture in the price basis;
  * the cycle row is MEASURED ON THE DECLARED WALK (cycle_declared_v19.json),
    not on the v16 canonicalised-state walk;
  * the success bound is the PRODUCTION CELL ALONE at 200 trials
    (0.404120 >= 0.39). The ladder is reported as supporting evidence and is no
    longer the object that carries the bound, so no cross-R / cross-order pooling
    is asked to transfer.

Claim is computed once here and copied verbatim by make_submission_v19.py.
"""
import json, math, os
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "evidence")

LOG2_V = 127.67
WALKERS = 1 << 20
STEPS = 2 ** (LOG2_V - 20)


def cp95_lower(k, n, alpha=0.025):
    """Exact Clopper-Pearson lower bound, same alpha convention as the runner."""
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
    store = json.load(open(os.path.join(OUT, "store_v19.json")))
    cyc = json.load(open(os.path.join(OUT, "cycle_declared_v19.json")))
    prod = json.load(open(os.path.join(OUT, "ladder_R2p20_W1_16p.json")))

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

    # ---- success: the PRODUCTION CELL ALONE carries the bound ----
    assert prod["trials"] >= 200 and prod["periods_covered"] >= 16.0
    prod_wins = prod["verified_successes"]
    prod_trials = prod["trials"]
    prod_lo = cp95_lower(prod_wins, prod_trials)
    assert abs(prod_lo - prod["cp95_lower"]) < 1e-9, (prod_lo, prod["cp95_lower"])
    assert prod_lo >= 0.39, prod_lo
    assert prod["useful_per_collision"] == 1.0

    # supporting ladder, reported but NOT pooled into the bound
    ladder = json.load(open(os.path.join(OUT, "period_ladder_v17.json")))["cells"]
    LADDER_LABELS = {"R64_W64_16p", "R256_W64_16p", "R1024_W64_16p", "R4096_W64_16p"}
    sup = [c for c in ladder if c["coverage_meets_16"] and c["label"] in LADDER_LABELS]
    pool_wins = sum(c["verified_successes"] for c in sup) + prod_wins
    pool_trials = sum(c["trials"] for c in sup) + prod_trials

    lam_meas = prod["collisions"] / prod_trials
    useful_frac = prod["useful_relations"] / max(1, prod["collisions"])

    full_step = 388 * fe + 4 * fa + 2 * fs + 30
    bulk_step = rows["affine_add"] + rows["batch_invert_share"] + 30
    infl = full_step / bulk_step
    trails_per_challenge = 2.0 * (1.0 + lam_meas)
    replay_rel = trails_per_challenge * 2 ** (LOG2_V - 20) * infl / (WALKERS * 2 ** (LOG2_V - 20))

    rel_cycle = cyc["pricing"]["relative_work_loss_charge"]
    log2_work = log2_rows + math.log2(1.0 + rel_cycle) + math.log2(1.0 + replay_rel)
    table_ir = 2 * (1 << 20) * 256 * 3 * fe
    log2_work += math.log2(1.0 + table_ir / 2 ** log2_work)

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
            "scatter256x_ir": 36.0,
            "store_io_ir_per_record": store["store_io_ir_per_record_per_pass"],
            "readback_ir_per_record": store["readback_ir_per_record_per_pass"],
            "per_pass_components": store["per_pass_ir_per_record"],
            "witness_ir": store["witness_ir_per_record"],
            "total_sort_ir_per_record": sort_ir,
            "capacity": store["capacity"],
            "basis": "ANALYTIC production store (evidence/store_v19.json): 36 Ir measured "
                     "loop body + external page-table walk, record store, page flush and "
                     "per-pass read-back charged per access. No bounded fixture is priced; "
                     "the shipped C program is a unit test of the addressing arithmetic only.",
            "v18_comparison": store["v18_comparison"],
        },
        "success": {
            "floor": 0.39,
            "claimed_object": "production cell ALONE (R = 2^20, 16 full selector periods, "
                              "fresh secret per trial, lambda 0.634)",
            "production_rung": {
                "label": prod["label"], "R": prod["R"], "trials": prod_trials,
                "wins": prod_wins, "rate": prod_wins / prod_trials, "cp95_lower": prod_lo,
                "collisions": prod["collisions"], "useful_per_collision": prod["useful_per_collision"],
                "periods": prod["periods_covered"],
                "seed_rule": prod["seed_rule"],
                "outcomes": prod.get("outcomes", ""),
                "prefix_check": prod.get("prefix_check"),
            },
            "supporting_ladder_not_the_bound": [
                {"label": c["label"], "R": c["R"], "periods": c["periods_covered"],
                 "trials": c["trials"], "wins": c["verified_successes"],
                 "measured": c["measured_success"], "cp95_lower": c["cp95_lower"]}
                for c in sup],
            "pooled_descriptive_only": {"wins": pool_wins, "trials": pool_trials,
                                        "note": "descriptive; NOT the claimed bound and not "
                                                "used to transfer success across R, group order "
                                                "or walker count"},
            "transfer_chain": {
                "measured": "per-trail success probability at 16 selector periods",
                "direction": "production affords 2^87.67 selector periods per walker, not 16, "
                             "and the within-trail collision probability is monotone "
                             "non-decreasing in the number of periods, so the measured "
                             "16-period rate is a LOWER bound on the production per-walker "
                             "probability, provided the per-period row-hazard is the same in "
                             "both regimes",
                "assumption_stated": "hazard equivalence: the row index is derived from the "
                                     "state representative XOR an advancing counter, and the "
                                     "representative is a uniform element of a PRIME-ORDER "
                                     "cyclic group in both regimes (fixture q prime; P-256 "
                                     "subgroup prime order). This link is an assumption, is "
                                     "the only unmeasured link, and is stated as such; no "
                                     "cross-R pooling is used to support it.",
                "not_claimed": "that a 49-bit fixture is isomorphic to the P-256 subgroup; the "
                               "fixture is a stand-in whose failure modes (trivial collisions, "
                               "cycle deaths, identity breaks) are measured and charged",
            },
            "note": "success = at least one coefficient-nontrivial collision that recovers x, "
                    "measured per trial on the global store; useful_per_collision = 1.0",
        },
        "policy": {
            "record_fields": "32B canonical key || 20-bit trail || 108-bit step",
            "rule": "every hit replays both trails from their stored identities; the replay "
                    "decides triviality and, if nontrivial, recovers x=-da/db; rejected hits "
                    "are charged their own replay",
            "measured_collisions_per_challenge": lam_meas,
            "measured_useful_fraction": useful_frac,
            "trail_pairs_charged_per_challenge": trails_per_challenge,
            "replay_inflation_vs_bulk": infl,
            "replay_relative_loss": replay_rel,
            "wasted_steps_per_return": 0.0,
            "cycle_charge_steps_per_death": 4.0,
        },
        "cycle": {
            "relative_loss": rel_cycle,
            "source": "cycle_declared_v19.json: MEASURED ON THE DECLARED WALK (true-point "
                      "stepping, canonical key only) at R=64/W=64/16 periods, 600 trails, "
                      "tracking on: 0 deaths; charged at the 95% CP upper bound on the "
                      "per-trail death rate divided by the 2^20 walkers.",
            "measured_deaths": cyc["deaths"], "trails": cyc["trails"],
            "death_rate_cp95_upper": cyc["death_rate_cp95_upper"],
            "v18_source": "cycle_channel_v16.json measured the CANONICALISED-STATE walk; that "
                          "number is superseded here and kept only as the control that the "
                          "tracker detects deaths (every v16 trail cycles at gap 2)",
            "gap2_condition": {
                "statement": "from eps2*(eps1*(S+v1)+v2) = S with eps1,eps2 in {+1,-1}: "
                             "v2 = -eps1*v1 (when eps2 = eps1) or v2 = -2*eps1*S - eps1*v1 "
                             "(when eps2 = -eps1), all mod q",
                "derivation": "expand eps2*(eps1*S + eps1*v1 + v2) = S and use 1/eps2 = eps2",
                "sign_correction": "v18 stated v2 = +eps1*v1 or +2*eps1*S + eps1*v1; those "
                                   "signs were WRONG (verdict v18, flaw 4). The corrected "
                                   "congruence above is the one this package uses.",
                "load_bearing": False,
                "why_not": "the shipped walk steps on the true point and canonicalises only "
                           "the stored key, so no state return is asserted; the cycle row is "
                           "the measurement above",
            },
        },
        "work": {"log2_visits": LOG2_V, "walkers": WALKERS, "steps_per_walker": STEPS,
                 "log2_rows": log2_rows, "log2_work_final": log2_work},
        "bounds": {"work_log2": claim, "claim": claim,
                   "storage_log2_peak": math.ceil((LOG2_V + math.log2(96)) * 1e4) / 1e4},
        "record": {"bytes": 48,
                   "layout": "32B canonical key || 16B identity (20-bit trail || 108-bit step)"},
        "claim_time_log2": claim,
        "claim_reduced_rational": rational,
    }
    json.dump(ev, open(os.path.join(OUT, "ledger_v19.json"), "w"), indent=1, sort_keys=True)
    json.dump({"claim": claim, "rational": rational, "log2_work": log2_work},
              open(os.path.join(OUT, "claim_v19.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({"per_visit": round(per_visit, 3), "sort_ir": round(sort_ir, 3),
                      "log2_work": log2_work, "claim": claim, "rational": rational,
                      "prod": {"trials": prod_trials, "wins": prod_wins, "cp95": prod_lo},
                      "pooled_descriptive": {"wins": pool_wins, "trials": pool_trials},
                      "cycle_rel": rel_cycle, "replay_rel": replay_rel}, indent=1))


if __name__ == "__main__":
    main()
