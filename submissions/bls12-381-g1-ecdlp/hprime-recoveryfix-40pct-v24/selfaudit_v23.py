#!/usr/bin/env python3
"""v23 self-audit: ONE canonical kappa, cited files exist, planted-variant control.

v22 defect: 'The manifest and argument specify kappa=1.337 ... build_certificate and
certificate.json instead use 1.231 ... These cannot silently be treated as the same
exact declared calculation.'  Here the certificate is the single source; the note and
argument are checked against it, and a planted edit in either document must FAIL.
"""
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
cert = json.load(open(HERE + "/certificate.json"))
note = open(HERE + "/public-note.md").read()
arg = open(HERE + "/argument.md").read()
kappa = cert["debit"]["kappa"]                       # THE designated field
res = {}

# ---- 1. every document states the certificate's kappa, and no other value
def kappas(text):
    return [float(x) for x in re.findall(r"(?:κ|kappa)\s*=\s*\**\s*([0-9]+\.[0-9]+)", text)]

res["cert_kappa"] = kappa
res["note_kappas"] = kappas(note)
res["argument_kappas"] = kappas(arg)
res["note_matches_cert"] = all(abs(k - kappa) < 5e-4 for k in res["note_kappas"]) and bool(res["note_kappas"])
res["argument_matches_cert"] = all(abs(k - kappa) < 5e-4 for k in res["argument_kappas"]) and bool(res["argument_kappas"])


def audit(note_text, arg_text):
    ok = all(abs(k - kappa) < 5e-4 for k in kappas(note_text)) and kappas(note_text)
    ok = ok and all(abs(k - kappa) < 5e-4 for k in kappas(arg_text)) and kappas(arg_text)
    return bool(ok)


res["audit_pass"] = audit(note, arg)
planted_note = note.replace("κ = %.3f" % kappa, "κ = %.3f" % (kappa + 0.001), 1)
res["planted_note_control_fails"] = (not audit(planted_note, arg)) if planted_note != note else None
planted_arg = arg.replace("%.3f" % kappa, "%.3f" % (kappa + 0.001), 1)
res["planted_argument_control_fails"] = (not audit(note, planted_arg)) if planted_arg != arg else None
res["mutation_control_pass"] = bool(res["planted_note_control_fails"] and
                                    res["planted_argument_control_fails"])

# ---- 2. the withdrawn claims are really gone
res["no_weight_convergence_claim"] = ("extrapolat" not in (note + arg).lower()) and ("monotone approach" not in (note + arg).lower())
res["no_route_zero_event_claim"] = cert["retry"]["route_vacuity"]["max_below_threshold"] is True
res["route_row_absent"] = "route_retry_worst" not in cert["rows"]
res["no_any_budget_wording"] = ("infeasible at any budget" not in (note + arg)) and ("out of reach at the claimed budget" in (note + arg).lower() or "at the claimed budget" in (note + arg).lower())
res["direct_success_null"] = cert["direct"]["success_final"] is None and cert["direct"]["certified"] is False

# ---- 3. every cited evidence file exists
sub = json.load(open(HERE + "/submission.json"))
missing = [f for f in sub["claims"][0]["evidence"] if not os.path.exists(os.path.join(HERE, f))]
res["evidence_files"] = len(sub["claims"][0]["evidence"])
res["missing_evidence"] = missing
res["all_evidence_present"] = not missing

# ---- 4. ledger/claim arithmetic
res["rows"] = len(cert["rows"])
res["W_log2"] = cert["raw_log2"]
res["claim"] = cert["claim_text"]
res["claim_above_W"] = float(cert["claim_text"]) > cert["raw_log2"]
res["claim_below_gate"] = float(cert["claim_text"]) < 116.8271

res["PASS"] = all([res["audit_pass"], res["mutation_control_pass"], res["no_weight_convergence_claim"],
                   res["route_row_absent"], res["no_any_budget_wording"], res["direct_success_null"],
                   res["all_evidence_present"], res["claim_above_W"], res["claim_below_gate"]])
json.dump(res, open(HERE + "/selfaudit-v23.json", "w"), indent=1)
print(json.dumps(res, indent=1))
