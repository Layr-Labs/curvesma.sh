"""v18 ledger: restored price anchors, carried-scatter sorter price, pooled
success bound that INCLUDES the measured production rung, cycle row unchanged.

Claim is computed here once and copied verbatim into submission.json by
make_submission_v18.py (never typed by hand).
"""
import json, math, os
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "evidence")

LOG2_V = 127.67
WALKERS = 1 << 20
STEPS = 2 ** (LOG2_V - 20)


def cp95_lower(k, n):
    if k == 0:
        return 0.0
    def tail(p):
        lp, lq = math.log(p), math.log1p(-p)
        s = 0.0
        lg = math.lgamma
        for i in range(k, n + 1):
            s += math.exp(lg(n + 1) - lg(i + 1) - lg(n - i + 1) + i * lp + (n - i) * lq)
        return s
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if tail(mid) < 0.025:
            lo = mid
        else:
            hi = mid
    return lo


def main():
    price = json.load(open(os.path.join(OUT, "insn_price.json")))
    fe = price["prices_ir"]["fe_mul"]
    fa = price["prices_ir"]["fe_add"]
    fs = price["prices_ir"]["fe_sub"]
    rad = json.load(open(os.path.join(OUT, "radix_v18.json")))
    cyc = json.load(open(os.path.join(OUT, "cycle_channel_v16.json")))

    per_pass = (rad["clear256_ir_per_bucket"] * 65536 / (1 << 20)
                + rad["count256_ir_per_record"]
                + rad["prefix256_ir_per_bucket"] * 65536 / (1 << 20)
                + rad["scatter256x_ir_per_record"])
    sort_ir = per_pass * rad["passes"] + rad["witness_ir_per_record"]

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

    # ---- success: the whole measured ladder, production rung INCLUDED ----
    ladder = json.load(open(os.path.join(OUT, "period_ladder_v17.json")))["cells"]
    LADDER_LABELS = {"R64_W64_16p", "R256_W64_16p", "R1024_W64_16p", "R4096_W64_16p"}
    cells = [c for c in ladder if c["coverage_meets_16"] and c["label"] in LADDER_LABELS]
    prod = json.load(open(os.path.join(OUT, "ladder_R2p20_W1_16p.json")))
    assert prod["periods_covered"] >= 16.0 and prod["coverage_meets_16"]
    assert prod["trials"] >= 100, prod["trials"]
    cells = cells + [prod]
    wins = sum(c["verified_successes"] for c in cells)
    trials = sum(c["trials"] for c in cells)
    pooled = wins / trials
    pooled_lo = cp95_lower(wins, trials)
    lam_meas = sum(c["collisions"] for c in cells) / trials
    useful_frac = (sum(c["useful_relations"] for c in cells)
                   / max(1, sum(c["collisions"] for c in cells)))

    full_step = 388 * fe + 4 * fa + 2 * fs + 30
    bulk_step = rows["affine_add"] + rows["batch_invert_share"] + 30
    infl = full_step / bulk_step
    trails_per_challenge = 2.0 * (1.0 + lam_meas)
    replay_rel = trails_per_challenge * 2 ** (LOG2_V - 20) * infl / (WALKERS * 2 ** (LOG2_V - 20))

    rel_cycle = cyc["pricing"]["relative_work_loss"]
    log2_work = log2_rows + math.log2(1.0 + rel_cycle) + math.log2(1.0 + replay_rel)
    table_ir = 2 * (1 << 20) * 256 * 3 * fe
    log2_work += math.log2(1.0 + table_ir / 2 ** log2_work)

    claim = math.ceil(log2_work * 1e4) / 1e4
    assert claim >= log2_work and pooled_lo >= 0.39, pooled_lo

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
            "scatter256x_ir": rad["scatter256x_ir_per_record"],
            "scatter256x_basis": rad["scatter256x_basis"],
            "clear256_ir_per_bucket": rad["clear256_ir_per_bucket"],
            "prefix256_ir_per_bucket": rad["prefix256_ir_per_bucket"],
            "witness_ir": rad["witness_ir_per_record"],
            "total_sort_ir_per_record": sort_ir,
            "addressing": rad["addressing"],
            "scope": "GLOBAL over the whole store, through a carried 256-bit cursor and an external chunked store",
            "executed": rad["run"],
        },
        "success": {
            "floor": 0.39,
            "ladder": [{"label": c["label"], "R": c["R"], "periods": c["periods_covered"],
                        "trials": c["trials"], "wins": c["verified_successes"],
                        "measured": c["measured_success"], "cp95_lower": c["cp95_lower"],
                        "collisions": c["collisions"],
                        "useful_per_collision": c["useful_per_collision"]} for c in cells],
            "pooled": {"wins": wins, "trials": trials, "rate": pooled,
                       "cp95_lower": pooled_lo, "includes_production_rung": True},
            "production_rung": {"label": prod["label"], "R": prod["R"],
                                "trials": prod["trials"], "wins": prod["verified_successes"],
                                "rate": prod["measured_success"], "cp95_lower": prod["cp95_lower"],
                                "collisions": prod["collisions"],
                                "useful_per_collision": prod["useful_per_collision"]},
            "min_cell_cp95_lower": min(c["cp95_lower"] for c in cells),
            "note": "success = at least one coefficient-nontrivial collision that recovers x, "
                    "measured per trial on the global store; useful_per_collision = 1.0 in every "
                    "ladder cell, so the success event is witness existence",
        },
        "policy": {
            "record_fields": "32B canonical key || 20-bit trail || 108-bit step; no coefficients",
            "rule": "every hit replays both trails from their stored identities; the replay decides "
                    "triviality (da=db=0) and, if nontrivial, recovers x=-da/db; rejected hits are "
                    "charged their own replay",
            "measured_collisions_per_challenge": lam_meas,
            "measured_useful_fraction": useful_frac,
            "trail_pairs_charged_per_challenge": trails_per_challenge,
            "replay_inflation_vs_bulk": infl,
            "replay_relative_loss": replay_rel,
            "wasted_steps_per_return": 0.0,
            "cycle_charge_steps_per_death": 4.0,
        },
        "cycle": {"relative_loss": rel_cycle,
                  "source": "cycle_channel_v16.json: MEASURED per-window loss at R=32/64/128 "
                            "(structural R^-2 2.283e-13 + additive R^-1 1.510e-7); it does NOT rest "
                            "on any unconditional 2-cycle identity",
                  "precondition_if_stated": "a canonicalised-STATE walk returns to S at gap 2 only "
                                            "if sign(v2) == sign(v1) (row repeat, probability ~1/R "
                                            "with an advancing counter) or v2 == 2S+v1 (mod q) in the "
                                            "negated case / v2 == -v1 or -2S-v1 in the non-negated "
                                            "case; v17/v18 step on the true point and canonicalise "
                                            "only the stored key, so no such return is claimed"},
        "work": {"log2_visits": LOG2_V, "walkers": WALKERS, "steps_per_walker": STEPS,
                 "log2_rows": log2_rows, "log2_work_final": log2_work},
        "bounds": {"work_log2": claim, "claim": claim,
                   "storage_log2_peak": math.ceil((LOG2_V + math.log2(96)) * 1e4) / 1e4},
        "record": {"bytes": 48,
                   "layout": "32B canonical key || 16B identity (20-bit trail || 108-bit step)"},
        "claim_time_log2": claim,
        "claim_reduced_rational": rational,
    }
    json.dump(ev, open(os.path.join(OUT, "ledger_v18.json"), "w"), indent=1, sort_keys=True)
    json.dump({"claim": claim, "rational": rational, "log2_work": log2_work},
              open(os.path.join(OUT, "claim_v18.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({"per_visit": round(per_visit, 3), "sort_ir": round(sort_ir, 3),
                      "log2_work": log2_work, "claim": claim, "rational": rational,
                      "pooled": {"wins": wins, "trials": trials, "rate": pooled,
                                 "cp95_lower": pooled_lo},
                      "prod": {"trials": prod["trials"], "wins": prod["verified_successes"],
                               "rate": prod["measured_success"], "cp95_lower": prod["cp95_lower"]},
                      "lam": lam_meas}, indent=1))


if __name__ == "__main__":
    main()
