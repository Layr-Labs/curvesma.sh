#!/usr/bin/env python3
"""v24 self-audit: certificate <-> documents consistency + planted-edit sensitivity.

Checks (all must pass):
  * claim digits in argument.md == public-note.md == certificate.json
  * kappa charge 1.231 in both documents == certificate debit.kappa
  * success row in both documents == certificate success.after
  * ledger recomputes: sum(rows)/2000 == W, claim >= W, claim < gate
  * direct census rows are internally consistent (solves/trials == q)
  * gate passed, kappa_168 labelled not-measured, retry priced per regime
  * PLANTED EDIT: perturb the claim by 1e-3 in a copy of argument.md and require
    the check to FAIL (otherwise the audit is vacuous)
"""
import json, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load():
    return (json.load(open(HERE + "/certificate.json")),
            open(HERE + "/argument.md").read(),
            open(HERE + "/public-note.md").read(),
            json.load(open(HERE + "/gate-v24.json")),
            json.load(open(HERE + "/direct-v24.json")))


def check(arg_text, note_text, cert, gate, direct, claim_override=None):
    claim = claim_override or (int(cert["claim"]["numerator"]) /
                               int(cert["claim"]["denominator"]))
    out = {}
    out["claim_in_argument"] = ("%.4f" % claim) in arg_text
    out["claim_in_note"] = ("%.4f" % claim) in note_text
    out["kappa_in_argument"] = ("%.3f" % cert["debit"]["kappa"]) in arg_text
    out["kappa_in_note"] = ("%.3f" % cert["debit"]["kappa"]) in note_text
    succ = cert["success"]["after"]
    out["success_in_argument"] = ("%.6f" % succ) in arg_text
    out["success_in_note"] = ("%.6f" % succ) in note_text
    rows = cert["rows"]
    W = sum(rows.values()) / 2000.0
    out["ledger_recomputes"] = abs(W - 2 ** cert["raw_log2"]) / W < 1e-6  # raw_log2 is rounded to 6 dp
    out["claim_above_W"] = claim > math.log2(W)
    out["claim_below_gate"] = claim < int(cert["gate"]["numerator"]) / int(cert["gate"]["denominator"])
    out["rows_11"] = len(rows) == 11
    out["gate_passed"] = bool(gate["passed"] and gate["ok"] == gate["n"] == 64)
    out["gate_old_path_zero"] = gate["v23_path_solves"] == 0
    out["success_measured"] = direct["summary"]["total_trials"] >= 512
    out["floor_flag_consistent"] = bool(cert["success"]["meets_floor"] ==
                                        (succ >= cert["success"]["claim_floor"]))
    out["below_floor_disclosed"] = bool(cert["success"]["success_below_floor"] ==
                                        (not cert["success"]["meets_floor"]))
    out["kappa_168_not_measured"] = cert["kappa_168"]["measured"] is False
    out["kappa_charged_dominates"] = bool(cert["kappa_168"]["charged_exceeds_fitted"])
    out["retry_two_regimes"] = len(cert["retry"]["per_regime"]) == 2
    out["retry_priced_by_worst"] = abs(
        cert["retry"]["worst_regime_used"]["q_gap_cp95_upper"] -
        max(c["q_gap_cp95_upper"] for c in cert["retry"]["per_regime"])) < 1e-12
    out["route_row_absent"] = "route_retry_worst" not in rows
    for g in direct["groups"]:
        q = g["solves"] / g["trials"]
        out["row_%d_q_consistent" % g["group"]] = abs(q - g["q"]) < 1e-12
        out["row_%d_on_curve" % g["group"]] = g["points_on_curve"] == g["records_per_trial"] * g["trials"]
    out["PASS"] = all(v for k, v in out.items() if k != "PASS" and isinstance(v, bool))
    return out


def main():
    cert, arg, note, gate, direct = load()
    base = check(arg, note, cert, gate, direct)
    # planted edit: +1e-3 on the claim must be caught by the document checks
    bad_claim = int(cert["claim"]["numerator"]) / int(cert["claim"]["denominator"]) + 1e-3
    planted = check(arg, note, cert, gate, direct, claim_override=bad_claim)
    base = {k: v for k, v in base.items() if k not in ("claim_in_argument", "claim_in_note",
                                                       "success_in_argument", "success_in_note")}
    base["planted_edit_detected"] = not (planted["claim_in_argument"] and planted["claim_in_note"])
    base["planted_edit_value"] = round(bad_claim, 6)
    base["evidence_files"] = len([f for f in os.listdir(HERE) if not f.startswith(".")])
    base["PASS"] = base["PASS"] and base["planted_edit_detected"]
    json.dump(base, open(HERE + "/selfaudit-v24.json", "w"), indent=1)
    print(json.dumps(base, indent=1))
    return 0 if base["PASS"] else 1


if __name__ == "__main__":
    sys.exit(main())
