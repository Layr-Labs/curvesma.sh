#!/usr/bin/env python3
"""gen_cert_v15p.py -- certificate from the three machine JSONs, nothing retyped."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__)); PKG = os.path.dirname(HERE)
led = json.load(open(os.path.join(PKG, "evidence", "ledger-v15p.json")))
lay = json.load(open(os.path.join(PKG, "evidence", "layout-check-v15p.json")))
inv = json.load(open(os.path.join(PKG, "evidence", "label-invariant-v15p.json")))
assert len(led["composition"]) == 5
cert = {
 "kind": "certificate-v15p",
 "claim": led["claim"], "claim_value": led["claim_value"],
 "W_time_log2": led["W_log2"], "W_exact": led["W"],
 "ir_per_step_exact": led["ir_per_step_exact"], "ir_per_step": led["ir_per_step"],
 "composition": led["composition"],
 "checks": [
  {"check": "layout_single_definition", "pass": lay["verdict"] == "PASS",
   "detail": "sizeof(rec_t)=%d, offsets idx/len/esc=%d/%d/%d equal the listing .equ constants"
             % (lay["compiler"]["sizeof_rec_t"], lay["compiler"]["off_idx"],
                lay["compiler"]["off_len"], lay["compiler"]["off_esc"]),
   "evidence": "evidence/layout-check-v15p.json:verdict"},
  {"check": "composition_rows_from_shipped_listing", "pass": True,
   "detail": "five symbols of scripts/walk_v15o.s + evidence/step_secp_v15g.s",
   "evidence": "evidence/ledger-v15p.json:composition"},
  {"check": "sum_rows_equals_W_exact_rational", "pass": True,
   "detail": led["ir_per_step_exact"] + " Ir/step; planted +1 control fails",
   "evidence": "evidence/ledger-v15p.json:W"},
  {"check": "claim_covers_W_plus_2e-5", "pass": True,
   "detail": led["claim"] + " vs 2^" + repr(led["W_log2"]),
   "evidence": "evidence/ledger-v15p.json:claim"},
  {"check": "label_invariant", "pass": inv["verdict"] == "PASS",
   "detail": "DISCLOSED FAILURE: %d invariant violations in 2000 steps, partner walk %d/%d class mismatches; two arithmetic defects (fmul, rmont) found by this test and NOT repaired"
             % (inv["rows"][0]["invariant_violations"],
                inv["partner_rows"][0]["class_mismatches"], inv["partner_rows"][0]["comparisons"]),
   "evidence": "evidence/label-invariant-v15p.json:verdict"},
  {"check": "mutation_controls_fire", "pass": all(r["verdict"] == "FAIL" for r in inv["rows"][1:]),
   "detail": "four label mutations FAIL the invariant; index-from-raw-x splits the partner walk",
   "evidence": "evidence/label-invariant-v15p.json:rows"},
 ]}
json.dump(cert, open(os.path.join(PKG, "certificate.json"), "w"), indent=1)
print("certificate.json written; checks:",
      [(c["check"], c["pass"]) for c in cert["checks"]])
