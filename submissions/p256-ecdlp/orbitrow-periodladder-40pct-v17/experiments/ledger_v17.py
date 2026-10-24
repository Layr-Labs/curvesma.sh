"""v17 ledger: global external sort with 256-bit counters, witness extraction,
production collision policy, cycle charge, and a claim rational taken from this
function's own output (never typed by hand).

v16 defect this fixes: the first-claim score in submission.json was typed as
641003/5000 = 128.2006 while the ledger's work was 128.2137; the shipped bound
was below the computed work.  Here `claim` is computed once, asserted >= work,
and written to evidence/claim_v17.json; make_submission_v17.py copies that value
verbatim into submission.json.

Prices come from evidence/radix_v17.json (objdump spans of the shipped arm64 -O2
build of experiments/radix_sort_v17.c) and evidence/insn_price.json (fe_mul
148 Ir, fe_add/fe_sub 27 Ir).
"""
import json, math, os, sys
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "evidence")

LOG2_V = 127.67
WALKERS = 1 << 20
STEPS = 2 ** (LOG2_V - 20)
LOG2_R = 256.0


def main():
    price = json.load(open(os.path.join(OUT, "insn_price.json")))
    fe = price["prices_ir"]["fe_mul"]
    fa = price["prices_ir"]["fe_add"]
    fs = price["prices_ir"]["fe_sub"]
    rad = json.load(open(os.path.join(OUT, "radix_v17.json")))
    cyc = json.load(open(os.path.join(OUT, "cycle_channel_v16.json")))

    # ---- per-pass sort price, all measured on the shipped build ----
    per_pass = (rad["clear256_ir_per_bucket"] * 65536 / (1 << 20)
                + rad["count256_ir_per_record"]
                + rad["prefix256_ir_per_bucket"] * 65536 / (1 << 20)
                + rad["scatter256_ir_per_record"])
    sort_ir = per_pass * rad["passes"] + rad["witness_ir_per_record"]

    rows = {
        "affine_add": 3 * fe + 4 * fa + 2 * fs,        # 2M+1S + chord/coord adds
        "batch_invert_share": 3 * fe,                   # 3 muls per record, one chain per batch
        "batch_chain_share": (385 * fe) / (1 << 20),    # 385-mul chain / batch
        "table_read_index": 10,
        "canon_sign_rule": 20,
        "record_append_48B": 12,
        "global_sort_and_witness": sort_ir,             # 16 x 256-bit LSD passes + witness emit
    }
    per_visit = sum(rows.values())
    units = per_visit / 1000.0
    log2_rows = LOG2_V + math.log2(units)

    # ---- collision-processing policy, from production record fields only ----
    # Records carry (32B canonical key || 20-bit trail || 108-bit step).  A hit is
    # discovered by the global sorted order; the record cannot say whether the
    # match is coefficient-trivial, so EVERY hit replays both trails and the
    # replay decides.  Census measures (a) hits per challenge, (b) the fraction
    # that survive the triviality filter.  Both are charged.
    agg = os.path.join(OUT, "period_ladder_v17.json")
    if os.path.exists(agg):
        cen = json.load(open(agg))
    else:
        import glob
        cells = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(OUT, "ladder_*.json")))]
        cen = {"cells": cells}
    cells = [c for c in cen["cells"] if c["coverage_meets_16"]]
    lam_meas = sum(c["collisions"] for c in cells) / sum(c["trials"] for c in cells)
    useful_frac = (sum(c["useful_relations"] for c in cells)
                   / max(1, sum(c["collisions"] for c in cells)))
    triv_frac = 1.0 - useful_frac
    full_step = 388 * fe + 4 * fa + 2 * fs + 30           # dependent step, own inverse
    bulk_step = rows["affine_add"] + rows["batch_invert_share"] + 30
    infl = full_step / bulk_step
    trails_per_challenge = 2.0 * (1.0 + lam_meas)   # pair for the recovering hit + a pair for every hit that may be a reject
    replay_rel = trails_per_challenge * 2 ** (LOG2_V - 20) * infl / (WALKERS * 2 ** (LOG2_V - 20))
    reject_rel = lam_meas * 2 ** (LOG2_V - 20) * 2 * infl / (WALKERS * 2 ** (LOG2_V - 20))

    # ---- cycle channel: no six-step allowance; it was never derived ----
    # The shipped walk steps on the true point (v17 fix): canon is a storage rule.
    # A return does not stop or restart the walker -- the next step uses the next
    # row because t increments -- so the wasted-step allowance is ZERO extra steps;
    # the only cost is that the 4 steps of a structural return are spent without
    # new progress, which is the cycle-channel relative loss already measured.
    rel_cycle = cyc["pricing"]["relative_work_loss"]
    waste_per_death = 0.0
    cyc_charge = cyc["pricing"]["per_4step_death_rate_prod"] * 4.0   # 4 spent steps, no extra

    rel_cycle_total = rel_cycle
    log2_work = log2_rows
    log2_work += math.log2(1.0 + rel_cycle_total)
    log2_work += math.log2(1.0 + replay_rel)

    # advice / table construction: 2R scalar mults, R = 2^20 rows, negligible but stated
    table_ir = 2 * (1 << 20) * 256 * 3 * fe
    log2_work += math.log2(1.0 + table_ir / 2 ** log2_work)

    claim = math.ceil(log2_work * 1e4) / 1e4
    assert claim >= log2_work

    # storage: live store 48V + one scratch array 48V (external sort) = 96V bytes
    peak_bytes = 96
    store_log2 = LOG2_V + math.log2(peak_bytes)
    queue_log2 = math.log2(WALKERS * 32)
    peak_log2 = math.log2(2 ** store_log2 + 2 ** queue_log2)

    num = int(round(claim * 10000))
    g = gcd(num, 10000)
    rational = {"numerator": num // g, "denominator": 10000 // g}

    ev = {
        "rows_ir_per_visit": {k: round(v, 6) for k, v in rows.items()},
        "instructions_per_visit": per_visit,
        "units_per_visit": units,
        "sort": {
            "per_pass_ir_per_record": per_pass,
            "passes": rad["passes"],
            "count256_ir": rad["count256_ir_per_record"],
            "scatter256_ir": rad["scatter256_ir_per_record"],
            "clear256_ir_per_bucket": rad["clear256_ir_per_bucket"],
            "prefix256_ir_per_bucket": rad["prefix256_ir_per_bucket"],
            "witness_ir": rad["witness_ir_per_record"],
            "total_sort_ir_per_record": sort_ir,
            "scope": "GLOBAL: the sort spans the whole store, so a hit can be between any two "
                     "visits of the challenge, not only inside one 2^20 batch",
            "counters": "4 x u64 (256-bit) per bucket on every pass; at V=2^127.67 a bucket "
                        "holds ~2^111.67 records, so no pass can use a 32-bit count",
            "executed": rad["run"],
        },
        "policy": {
            "record_fields": "32B canonical key || 20-bit trail || 108-bit step; no coefficients",
            "rule": "every hit replays both trails from their stored identities; the replay "
                    "decides triviality (da=db=0) and, if nontrivial, recovers x=-da/db; "
                    "rejected hits are charged their own replay",
            "measured_collisions_per_challenge": lam_meas,
            "measured_useful_fraction": useful_frac,
            "measured_trivial_fraction": triv_frac,
            "trail_pairs_charged_per_challenge": trails_per_challenge,
            "replay_inflation_vs_bulk": infl,
            "replay_relative_loss": replay_rel,
            "reject_relative_loss": reject_rel,
            "scanner": "rows are drawn uniformly: r = (x(S) ^ t*2654435761 ^ trail*7919) mod R, "
                       "so tuple enumeration is uniform over R and carries no weight beyond "
                       "the measured hit rate",
            "wasted_steps_per_return": waste_per_death,
            "cycle_charge_steps_per_death": 4.0,
            "why_no_six_steps": "v16's six-step allowance was not derived from an implemented "
                                "escape rule and the escape was deleted; here a return costs no "
                                "extra step because the walk continues from the true point and "
                                "t advances, so only the 4 spent steps are charged",
        },
        "cycle": {"relative_loss": rel_cycle_total,
                  "source": "cycle_channel_v16.json: structural R^-2 2.283e-13/window + "
                            "additive R^-1 1.510e-7/window, measured at R=32/64/128"},
        "work": {"log2_visits": LOG2_V, "walkers": WALKERS, "steps_per_walker": STEPS,
                 "log2_rows": log2_rows, "log2_work_final": log2_work},
        "bounds": {"work_log2": claim, "storage_log2_peak": math.ceil(peak_log2 * 1e4) / 1e4,
                   "claim": claim},
        "success": {
            "floor": 0.39,
            "ladder": [{"label": c["label"], "R": c["R"], "periods": c["periods_covered"],
                        "trials": c["trials"], "measured": c["measured_success"],
                        "cp95_lower": c["cp95_lower"], "collisions": c["collisions"],
                        "useful_per_collision": c["useful_per_collision"]}
                       for c in cells],
            "min_cp95_lower": min(c["cp95_lower"] for c in cells),
            "note": "success = at least one coefficient-nontrivial collision that recovers x, "
                    "measured per trial on the global store; useful_per_collision = 1.0 in "
                    "every ladder cell, so the success event is witness existence, not a mean "
                    "collision count",
        },
        "record": {"bytes": 48,
                   "layout": "32B canonical key || 16B identity (20-bit trail || 108-bit step); "
                             "no coefficient field, so a hit's triviality is only knowable at replay"},
        "claim_time_log2": claim,
        "claim_reduced_rational": rational,
    }
    json.dump(ev, open(os.path.join(OUT, "ledger_v17.json"), "w"), indent=1, sort_keys=True)
    json.dump({"claim": claim, "rational": rational, "log2_work": log2_work},
              open(os.path.join(OUT, "claim_v17.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({"per_visit": round(per_visit, 3), "units": round(units, 6),
                      "sort_ir": round(sort_ir, 3), "log2_work": log2_work, "claim": claim,
                      "rational": rational, "replay_rel": replay_rel,
                      "lam_meas": lam_meas, "useful_frac": useful_frac}, indent=1))


if __name__ == "__main__":
    main()
