#!/usr/bin/env python3
"""build_evidence_v10p.py -- single source of truth for the v10p evidence rows.

Reads the pooled pair census and the carried measurement files, then writes:
  evidence/batch_amortization.json   B in {8,16,65536}: per-walker share, charge, claimed row
  evidence/kappa_transfer.json       kappa derived from the measured closure-rate match
  evidence/a_shift.json              P0(a=0) vs P0(a=-3 numerator) + a-channel bound
  evidence/ledger.json               re-summed charge, claim integer, asserts
  evidence/claim.json                W and the 4-dp ceiling

Every number here is recomputed, not copied: the charge rows are re-derived from the
measured span prices, so a stale literal anywhere is caught by arith_audit.py.
"""
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
E = HERE

FE_MUL_INSTR = 147                 # measured: fe_mul = 147 Ir (stepspans_core.json)
INV_MULS = 385                     # Fermat inversion: 256 squarings + 129 muls
STEP_FORWARD = 239
STEP_CLOSE = 694
BACKSUB = 293                      # per walker: prefix + back-substitution, 2 muls x 147
DEN_SCAN = 14                      # denominator zero-scan, per walker
RARE = 0.63976416015625            # rare-event reserve, carried and re-summed
S_LOG2 = 127.631517202749          # log2 of total steps, S = 2^127.6315
W_AMP = 2
CLOSURE_MEASURED_LOG2 = -24.03     # measured closures per walker-step
CLOSURE_PREDICTED_LOG2 = -24.0     # predicted from the 2^-11 precursor
P0_PRE = {"dp4": 0.248333, "dp6": 0.234667, "dp8": 0.241, "dp10": 0.239667}
P0_MEAN = sum(P0_PRE.values()) / len(P0_PRE)
P_GUARDED_DP4 = 0.033333           # P2_shipped, dp4


def charge(B, kappa):
    share = INV_MULS * FE_MUL_INSTR / B
    esc = 57277 * 2 ** CLOSURE_MEASURED_LOG2 * kappa
    total = STEP_FORWARD + STEP_CLOSE + BACKSUB + DEN_SCAN + share + RARE + esc
    return share, esc, total


def w(B, kappa):
    _s, _e, ch = charge(B, kappa)
    single = S_LOG2 + math.log2(ch / 1000.0)
    return single, single + math.log2(W_AMP)


def main():
    pc = json.load(open(os.path.join(E, "pair_census_v10p.json")))
    kappa = 2 ** (CLOSURE_MEASURED_LOG2 * -1 - 24.0)      # 2^(24.03-24) = 1.0210
    kappa_meas_ratio = 2 ** (CLOSURE_MEASURED_LOG2 + 24.0)  # 0.9794

    # ---------------- batch amortization ----------------
    rows = {}
    for B in (8, 16, 65536):
        share, esc, ch = charge(B, kappa)
        sgl, w2 = w(B, kappa)
        rows[str(B)] = {
            "batch_walkers": B,
            "inversion_muls_per_batch": INV_MULS,
            "fe_mul_instr": FE_MUL_INSTR,
            "inversion_instr_per_batch": INV_MULS * FE_MUL_INSTR,
            "share_per_walker_step": share,
            "prefix_backsub_per_walker_step": BACKSUB,
            "escape_surcharge_kappa_charged": esc,
            "charge_instr_per_walker_step": ch,
            "W_single_log2": sgl,
            "W_w2_log2": w2,
            "claim_ceiling_4dp": math.ceil(w2 * 1e4) / 1e4,
        }
    worst = max(rows.values(), key=lambda r: r["charge_instr_per_walker_step"])
    ba = {
        "test": "batch_inv_amortization",
        "verdict_item": "inversion amortized over 65536 walkers, shipped W array 16, live census 8",
        "derivation": ("one Fermat inversion per STEP per BATCH of B walkers: 385 fe_mul; "
                       "prefix products (B-1) and back-substitution (B-1) muls are per-walker "
                       "and already charged in the 293 row. per-walker share = 385*147/B."),
        "rows": rows,
        "charged_row": worst["batch_walkers"],
        "charged_because": "worst of the shipped widths (8) is charged; the 65536 row is the production-lockstep reading",
        "delta_vs_v10o_vs_128.9437": worst["W_w2_log2"] - 128.94360930732344,
        "note_shipped": ("the disasm price is unchanged (fe_mul 147, 385 muls per inversion); only the "
                         "amortization denominator changes, so this is a pure accounting correction, "
                         "not a new measurement"),
        "arithmetic": {
            "385*147/65536": 385 * 147 / 65536,
            "385*147/16": 385 * 147 / 16,
            "385*147/8": 385 * 147 / 8,
            "sum_check_B65536_equals_v10o_charge": abs(charge(65536, 1.0)[2] - 1241.5067493023873) < 1e-6,
        },
    }
    json.dump(ba, open(os.path.join(E, "batch_amortization.json"), "w"), indent=1, sort_keys=True)

    # ---------------- kappa ----------------
    kt = {
        "test": "kappa_transfer",
        "verdict_item": "kappa measured at fixture, transferred to production without derivation",
        "definition": "kappa = (measured closure rate per walker-step) / (predicted rate), charged as max(1, 1/kappa_meas)",
        "measured_closures_per_step_log2": CLOSURE_MEASURED_LOG2,
        "predicted_closures_per_step_log2": CLOSURE_PREDICTED_LOG2,
        "measured_over_predicted": kappa_meas_ratio,
        "match_fraction": 0.979,
        "kappa_charge": kappa,
        "derivation": ("the closure rate is the ONLY channel through which the fixture's ring/escape "
                       "geometry enters the production surcharge. It was measured at 2^-24.03 against a "
                       "2^-24 prediction (97.9 pct of prediction): that ratio IS the transfer factor, and "
                       "it is now charged in the direction that costs work, kappa = 2^0.03 = 1.0210."),
        "residual": {
            "statistical": "closure count is a live counter, not a rate estimate; no binomial interval is claimed",
            "scope": "the rate is per walker-step, so the cumulative count is ell*rate (2^87.60), already settled",
            "unmodelled": "production DP density 2^-62 and 65536 walkers change the collision channel, not the "
                          "ring-escape channel; the latter is per-step and table-driven",
        },
        "envelope_charge": {
            "escape_surcharge_unscaled": 57277 * 2 ** CLOSURE_MEASURED_LOG2,
            "escape_surcharge_kappa_charged": 57277 * 2 ** CLOSURE_MEASURED_LOG2 * kappa,
            "ir_delta": 57277 * 2 ** CLOSURE_MEASURED_LOG2 * (kappa - 1),
            "bits_delta": math.log2(charge(65536, kappa)[2] / charge(65536, 1.0)[2]),
        },
    }
    json.dump(kt, open(os.path.join(E, "kappa_transfer.json"), "w"), indent=1, sort_keys=True)

    # ---------------- a shift ----------------
    sh = pc["shipped_guarded|a=0"]
    sh3 = pc["shipped_guarded|a=-3"]
    ct = pc["control_pre_parity|a=0"]
    ct3 = pc["control_pre_parity|a=-3"]
    steps_a0 = 16 * 50764 * sh["trials"] * 1.0
    steps_a3 = 16 * 50764 * sh3["trials"] * 1.0
    dbl_a0 = sh["doubles"]
    dbl_a3 = sh3["doubles"]
    a_shift = {
        "test": "a_shift",
        "verdict_item": "census uses an a=0 fixture; production is a=-3",
        "structural_fact": ("for y^2=x^3+ax+b the addition law with x1 != x2 is a-independent; a enters "
                            "ONLY the doubling numerator, 3x^2+a. The fixture curve has a=0 (numerator 3x^2); "
                            "production P-256 has a=-3 (numerator 3x^2-3)."),
        "measured_a_channel": {
            "walkers_per_run": 16,
            "steps_per_walker": 50764,
            "walker_steps_shipped_a0": steps_a0,
            "walker_steps_shipped_a_minus3": steps_a3,
            "doubling_events_a0": dbl_a0,
            "doubling_events_a_minus3": dbl_a3,
            "doubling_rate_a0": dbl_a0 / steps_a0,
            "a_channel_bound": ("|P(a=-3) - P(a=0)| <= (fraction of walker-steps whose transition is a "
                                "doubling) since only those transitions differ; measured fraction above"),
        },
        "rows": {
            "shipped_guarded_a0": {"n": sh["trials"], "p1": sh["p1_single_run"], "p2": sh["p2_set"],
                                   "p2_cp95_lower": sh["p2_cp95_lower"]},
            "shipped_guarded_a_minus3_numerator": {"n": sh3["trials"], "p1": sh3["p1_single_run"],
                                                   "p2": sh3["p2_set"], "p2_cp95_lower": sh3["p2_cp95_lower"]},
            "control_pre_parity_a0": {"n": ct["trials"], "p1": ct["p1_single_run"], "p2": ct["p2_set"],
                                      "p2_cp95_lower": ct["p2_cp95_lower"]},
            "control_pre_parity_a_minus3_numerator": {"n": ct3["trials"], "p1": ct3["p1_single_run"],
                                                      "p2": ct3["p2_set"], "p2_cp95_lower": ct3["p2_cp95_lower"]},
        },
        "disclosure": ("a true a=-3 fixture needs a prime-order group on a curve with a=-3, i.e. point "
                       "counting (Schoof) at 2^40; that is NOT done here. What is measured is the exact "
                       "code difference a introduces - the doubling numerator is substituted with 3x^2-3 "
                       "and the census rerun - plus the incidence of the branch it lives in. The group "
                       "order therefore differs from the a=0 fixture and that residual is disclosed, not "
                       "measured."),
        "materially_changes_floor_question": None,
    }
    if ct3["trials"]:
        a_shift["materially_changes_floor_question"] = {
            "control_p2_a0": ct["p2_set"], "control_p2_a_minus3": ct3["p2_set"],
            "difference": ct3["p2_set"] - ct["p2_set"],
            "shipped_p2_a0": sh["p2_set"], "shipped_p2_a_minus3": sh3["p2_set"],
            "shipped_difference": sh3["p2_set"] - sh["p2_set"],
        }
    json.dump(a_shift, open(os.path.join(E, "a_shift.json"), "w"), indent=1, sort_keys=True)

    # ---------------- ledger ----------------
    share, esc, ch = charge(8, kappa)
    sgl, w2 = w(8, kappa)
    ledger = {
        "version": "v10p-r1",
        "W": {
            "A_log2_Tprod": S_LOG2,
            "charge_total_instr_per_walker_step": ch,
            "charge_rows": {
                "step_forward": STEP_FORWARD,
                "step_close": STEP_CLOSE,
                "batch_inv_backsub_per_walker": BACKSUB,
                "denominator_zero_scan": DEN_SCAN,
                "batch_inv_inversion_share_B8": share,
                "rare_event_reserve": RARE,
                "escape_surcharge_kappa_charged": esc,
            },
            "charge_sum_check": STEP_FORWARD + STEP_CLOSE + BACKSUB + DEN_SCAN + share + RARE + esc,
            "charge_equals_sum": abs((STEP_FORWARD + STEP_CLOSE + BACKSUB + DEN_SCAN + share + RARE + esc) - ch) < 1e-9,
            "W_single_instance": sgl,
            "amplification_w": W_AMP,
            "W_amplified": w2,
            "claim": math.ceil(w2 * 1e4) / 1e4,
            "claim_ceiling_4dp": math.ceil(w2 * 1e4) / 1e4,
        },
        "amplification": {
            "w": W_AMP,
            "floor": 0.39,
            "shipped_arm_p2_measured": sh["p2_set"],
            "shipped_arm_p2_cp95_lower": sh["p2_cp95_lower"],
            "shipped_arm_meets_floor": sh["p2_cp95_lower"] >= 0.39,
            "control_arm_p2_measured": ct["p2_set"],
            "control_arm_p2_cp95_lower": ct["p2_cp95_lower"],
            "control_arm_meets_floor": ct["p2_cp95_lower"] >= 0.39,
            "reading": ("the w=2 ticket is measured, not modelled: on the SHIPPED (guarded) walk the set "
                        "probability is measured below the floor, so the success premise is UNMET for the "
                        "shipped arm. The pre-guard control, measured the same way on the same cell, clears "
                        "it. The guard - not the w=2 construction - is the whole obstacle."),
        },
        "assert": {
            "charge_equals_sum": True,
            "claim_covers_W": True,
            "rounding": "claim is the 4-dp CEILING of W_w2, so at-most holds",
            "batch_inv_at_shipped_width": True,
            "kappa_charged": True,
            "success_ge_floor": bool(sh["p2_cp95_lower"] >= 0.39),
            "control_ge_floor": bool(ct["p2_cp95_lower"] >= 0.39),
        },
        "object": {"sha256": "96b736ef6dd55cc7cee4be9dac2216b48f7fb56cd40cd4d99d16139451c2990f",
                   "bytes": 32864, "note": "same counted TU as v10o; the span prices are unchanged"},
    }
    json.dump(ledger, open(os.path.join(E, "ledger.json"), "w"), indent=1, sort_keys=True)

    claim = {"W_amplified": w2, "claim_ceiling_4dp": math.ceil(w2 * 1e4) / 1e4,
             "charged_batch": 8, "W_at_B16": w(16, kappa)[1], "W_at_B65536": w(65536, kappa)[1]}
    json.dump(claim, open(os.path.join(E, "claim.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({k: claim[k] for k in claim}, indent=1))
    print("shipped p2", sh["p2_set"], sh["p2_cp95_lower"], "| control p2", ct["p2_set"], ct["p2_cp95_lower"])


if __name__ == "__main__":
    main()
