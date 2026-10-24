#!/usr/bin/env python3
"""make_reproduction.py -- v18: regeneration.json is GENERATED, never typed.

v17 shipped a reproduction.json whose `expected` block still carried the WITHDRAWN
v17a numbers (success 0.560408, p_set 0.5771875, kernel span 95, step_insn 152,
claim 127.926).  Every one of those contradicted the submitted claim, and a
reviewer reading the reproduction file would have found a second, independent
inconsistency.  The fix is structural: this script reads the evidence files and
ledger.py and WRITES reproduction.json, so the expected block cannot drift from the
shipped evidence.  Run it after every evidence regeneration.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, D)
sys.path.insert(0, os.path.join(D, "scripts"))

import ledger  # noqa: E402


def j(rel):
    with open(os.path.join(D, "evidence", rel)) as f:
        return json.load(f)


def main():
    pair = j("pair-census.json")
    c20 = j("census-v20n.json")
    c11 = j("census-seed11.json")
    lag3 = j("lag3-firing.json")
    rc = j("restart-count.json")
    lad = j("transfer-ladder.json")
    tr = {t["arm"]: t for t in lad["arms"]}
    design = max(tr.values(), key=lambda t: (t["classes_log2"], t["sets"]))
    out = {
        "commands": [
            "python3 scripts/lag3_fixture.py    # planted 3-cycle fixture -> evidence/lag3-firing.json",
            "python3 scripts/absorb_v20n.py     # exhaustive two-size absorption count -> evidence/absorb-direct.json",
            "BT=3200 SEED=20250 python3 scripts/census_v20n.py   # 3 arms x 3200 trials -> evidence/census-v20n.json",
            "BT=3200 SEED=11 python3 scripts/census_v20n.py      # second seed block -> evidence/census-seed11.json",
            "ARM=lag3 BT=3200 python3 scripts/pair_census.py     # 3200 sets x 2 runs, same challenge -> evidence/pair-census.json",
            "ARM=lag1 BT=3200 python3 scripts/pair_census.py     # lag-1 control pair arm -> evidence/pair-census-lag1.json",
            "ARMS=L24,L28,L32,L36 python3 scripts/transfer_arm.py  # ORDER LADDER, same-challenge w=2 -> evidence/transfer-arm-*.json",
            "python3 scripts/restart_count.py   # fresh-start rows at C=1000, 4-limb counts",
            "python3 scripts/make_reproduction.py  # this file, generated from the evidence",
            "python3 ledger.py                  # exact Fractions, log2 W2, claim, floor decision",
        ],
        "environment": {
            "python": ("/opt/homebrew/opt/python@3.14/bin/python3.14 3.14.8 (numpy 2.4.4) -- "
                       "NOT /usr/bin/python3, and NOT the login-shell python3, which has no "
                       "numpy on this host; the background-shell trap is real"),
            "cost_model": "tracks/docs/COST_MODEL.md, C = 1000 on secp256k1-ecdlp (the bls track's C = 2000 is used nowhere)",
            "disasm": ("WITHDRAWN.  No disassembly measurement is claimed.  scripts/kernel/cios4.c "
                       "failed its own full-width reference test twice (verdicts b86f973b, "
                       "cd18d785) and M_insn = 97 is a STRUCTURAL count only."),
        },
        "entrypoint": "ledger.py",
        "expected": {
            "claim": "%.3f" % float(ledger.CLAIM),
            "log2_work": "%.8f" % ledger.LOG_W,
            "log2_work_one_run": "%.8f" % ledger.log2f(ledger.W),
            "w2_factor": ledger.W2_FACTOR,
            "step_insn": ledger.STEP_INSN,
            "step_price": float(ledger.STEP_PRICE),
            "sum_of_step_rows": sum(n for _o, n in ledger.ROWS),
            "kernel_M_insn_structural": ledger.M_insn,
            "kernel_disassembly_claimed": False,
            "window_products_v1": 8,
            "success_design": round(ledger.SUCCESS, 6),
            "success_floor": 0.39,
            "floor_met": bool(ledger.SUCCESS >= 0.39),
            "transfer_design_arm": design["arm"],
            "transfer_design_order_log2": round(design["classes_log2"], 3),
            "transfer_design_sets": design["sets"],
            "transfer_ladder_p_set": {t["arm"]: t["p_set_measured"] for t in tr.values()},
            "transfer_ladder_cp95": {t["arm"]: t["cp95_set_lower_f"] for t in tr.values()},
            "transfer_band_low": round(ledger.T_TOP_MIN, 6),
            "transfer_band_high": round(ledger.T_TOP_MAX, 6),
            "transfer_order_slope_per_doubling": float("%.3e" % ledger.T_SLOPE),
            "équiv_extended_vs_shipped": "evidence/transfer-equiv.json (bit-exact trace)",
            "p_set_surrogate_3200": pair["p_set_measured"],
            "p_set_formula": pair["p_set_formula"],
            "cp95_w2_sets_surrogate": pair["cp95_set_lower_f"],
            "cp_helper_calibration_coverage": pair["calibration"]["coverage"],
            "restart_armed_trials": ledger.RESTART_TRIALS,
            "restart_events": ledger.RESTART_EVENTS,
            "restart_loss_fraction": float(ledger.RESTART_LOSS_FRACTION),
            "fresh_ref": rc["fresh_ref"],
            "lag3_success_block20250": c20["arms"]["escaped_lag3"]["success_rate"],
            "lag1_success_block20250": c20["arms"]["escaped_lag1"]["success_rate"],
            "lag3_success_block11": c11["arms"]["escaped_lag3"]["success_rate"],
            "three_cycle_lag1": c20["cycle_length"]["escaped_lag1"]["len3"],
            "three_cycle_lag3": c20["cycle_length"]["escaped_lag3"]["len3"],
            "lag3_fires_block20250": c20["arms"]["escaped_lag3"]["escape_fires_lag3"],
            "planted_period3_starts": lag3.get("planted_period3_starts"),
            "planted_lag3_channel_fires": lag3["lag3_fires_lag3_channel"],
            "planted_lag1_channel_fires_on_same_starts": lag3["lag1_fires_lag3_channel"],
            "invariant_violations": c20["arms"]["escaped_lag3"]["invariant_violations_sampled"],
        },
        "test": "reproduction",
    }
    with open(os.path.join(D, "reproduction.json"), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print(json.dumps(out["expected"], indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
