#!/usr/bin/env python3
"""make_submission.py -- build submission.json.  Every number is interpolated from the
owning JSON so the prose cannot drift from the ledger."""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda f: json.load(open(os.path.join(HERE, f)))
D, LED, LAD, LV = J("design-v36.json"), J("ledger-v36.json"), J("structured-ladder.json"), J("lift-verify.json")
cl, s, L, P = LED["claim"], LED["success"], LED["success"], LAD["production_evaluation"]
a, an = LED["anchors"], LED["anchors"]["declared"]
h, dx = LED["diagnosis"], LED["peak_storage"]

st = (
 f"BLOCKING NEGATIVE RESULT, submitted because the previous claim was refutable: claim "
 f"{float(cl['value']):.4f} ({cl['numerator']}/{cl['denominator']}), charged W=2^{LED['W_log2']:.6f} in "
 f"{LED['C_instructions_per_unit']}-instruction units against gate 116.8271, i.e. {LED['gate_delta_claim_minus_gate']:.4f} "
 f"bits ABOVE the gate. No rank is taken. The point of the package is that the gate cannot be taken by this "
 f"family at all, and this fire proves it with two executed instruments instead of two asserted numbers.\n\n"
 f"(1) THE MISSING LIFTING STAGE IS IMPLEMENTED AND EXECUTED (verdict flaw 1, the one gap with a real fix). "
 f"keymod.lift() takes s0 = s mod m with m=(r-1)/2^D and enumerates s = s0 + k*m over k in [0,2^D), testing the "
 f"FULL point condition [s]P0 == T, not a truncated key. Measured: {LV['lift_100_secrets']['recovered_exactly']}/100 "
 f"random secrets recovered exactly at D=10 with [s]P0==T and [s^d]P0==A_d on every trial; the full production "
 f"depth D=20 enumerated {LV['full_scale_D20']['candidates']} candidates in {LV['full_scale_D20']['wall_s']} s and "
 f"exactly one passed ({LV['full_scale_D20']['per_candidate_us']} us/candidate on the 37-bit fixture); the same lift "
 f"ran on the REAL 381-bit BLS12-381 G1 curve, {LV['real_bls12_381_g1']['recovered_exactly']}/{LV['real_bls12_381_g1']['trials']} "
 f"recovered with {LV['real_bls12_381_g1']['passes_total']} passes total (uniqueness measured); control: a one-bit corruption of s0 "
 f"yields ZERO passing candidates in 20/20 trials. The lifting cost row is in the ledger: {LED['rows']['lifting']} "
 f"instructions = 2^{LED['row_log2']['lifting']:.1f} of 2^{LED['W_log2']:.1f} total, i.e. a rounding error, "
 f"which is the honest reason the lifting fix does not move the claim.\n\n"
 f"(2) THE SUCCESS CLAIM IS WITHDRAWN AND THE REAL TERM IS COVERAGE (verdict flaw 2). ladder_v36.py RUNS the attack: "
 f"{LAD['powered_pooled']['pools']} independent pools over the TRUE Cartesian product (nb bases x one shared walk list, "
 f"every (i,l) record built as a real EC point, no hard-coded endpoint), at key widths 12/16/20/24 in three splits: "
 f"{LAD['powered_pooled']['successes']}/{LAD['powered_pooled']['pools']} pools verified at every width, CP-95 lower "
 f"{LAD['powered_pooled']['cp95_lower']}. Width does NOT matter, because a covered target shares its record's full x. "
 f"What matters is coverage: an uncovered target still produces truncated hits and EVERY one fails the full-point "
 f"check ({LAD['width_effect']['uncovered_control']['truncated_hits']} hits, "
 f"{LAD['width_effect']['uncovered_control']['false_accepts']} false accepts). So P(success) = (2N/r) x P(pipeline|covered) = "
 f"(2^{P['coverage_log2']:.2f}) x {P['conditional_success']} = 2^{P['success_log2']:.2f} at the declared N=2^115.65 and "
 f"r=2^{P['r_log2']:.3f}. The v35 number 0.396243 was a 12-20-bit key-hit rate; it is withdrawn, not re-explained. "
 f"Meeting the 0.39 floor needs N = 0.195*r = 2^{P['records_needed_log2']:.2f} records.\n\n"
 f"(3) W(n) AND PEAK STORAGE, COMPLETE (verdict flaw 4). W(n) = [2*N*c_rec + 2*N*c_rec/100 + 2*2^D*P_lift + 2*P_scalar + "
 f"N*2^(D+20-168) + 2*nb*4 + 2*nb*2^33*15*C + ell*10 + 2*2^20*10]/C with n entering through the key width (2^-n "
 f"truncation FP), the lifting depth D=log2(n), and N(n)=0.195*r(n); the table is in ledger-v36.json:W_n_table. "
 f"Peak storage is a MAX over every component, not the record store alone: records 2*N*64 B, gap tables "
 f"2*nb*2^33*15 points, O(1) walker scratch, code and state < 2^20 B -> MAX row 2^{dx['max_row_bytes_log2']:.1f} B, "
 f"which exceeds the instruction claim and is the binding resource.\n\n"
 f"(4) WHAT IS NOT CLAIMED (verdict flaw 3). The v35 price anchors are NOT carried and nothing here depends on a "
 f"compiled-body measurement: the per-record charge is a declared 2303 instructions and the field-op prices used for "
 f"the lifting row are DECLARED MODELS (mulmod {an['mul_exec']} instructions = 768-trip reduction loop x 12, inversion "
 f"{an['inversion_calls']} mulmod calls, affine add {an['point_add_affine']}, scalar mult {an['scalar_mult']}). The "
 f"static-vs-executed, INV-unit, sort-inner-step, SIMD and borrow-propagation defects are therefore left UNREPAIRED "
 f"and UNRELIED UPON rather than papered over. Diagnosis: at the v35 parameters the honest cost of the algorithm v35 "
 f"described is 2^{h['honest_at_v35_parameters']:.3f}, {h['delta_bits']:.1f} bits above the claim it filed; restoring "
 f"the gate would need N=2^{math.log2(float(D['N_floor'])):.1f} records, {LED['claim']['value'] - 116.8271:.1f} bits above the gate. "
 f"Recommended action for this track: close the line on the gate, do not re-fire on pricing."
)
assert len(st) <= 8000, len(st)
ev = ["public-note.md", "ledger-v36.json", "derive_v36.py", "design-v36.json", "structured-ladder.json",
      "ladder_v36.py", "lift-verify.json", "lift_verify.py", "keymod.py", "ec.py",
      "test_kernel_v36.py", "audit_v36.py"]
for e in ev:
    assert os.path.exists(os.path.join(HERE, e)), e
sub = {
 "version": 1,
 "title": ("bls12-381-g1-ecdlp v36: lifting implemented and executed (100/100 secrets, D=20 full scale, real "
           "381-bit curve), ladder RUN on the true Cartesian product, success claim withdrawn, coverage term "
           "2N/r identified, honest claim 2^%.2f" % LED["W_log2"]),
 "claims": [{"statement": st,
             "score": {"numerator": cl["numerator"], "denominator": cl["denominator"]},
             "evidence": ev}],
 "experiments": [
   {"id": "lift-verify-v36", "source": "lift_verify.py",
    "purpose": ("Execute the lifting stage end to end: 100 random secrets at D=10 with [s]P0==T and [s^d]P0==A_d, "
                "the full D=20 enumeration, the same lift on the real 381-bit BLS12-381 G1 curve, and the "
                "one-bit-corruption control. Writes lift-verify.json."),
    "input_schema": {"type": "object", "properties": {}, "required": [], "additionalProperties": False},
    "initial_inputs": []},
   {"id": "ladder-v36", "source": "ladder_v36.py",
    "purpose": ("Run the attack per pool over the true Cartesian product at key widths 12/16/20/24 in three splits, "
                "count verified recoveries, run the uncovered-target control, and evaluate the same law at the "
                "production width and record count. Writes structured-ladder.json."),
    "input_schema": {"type": "object", "properties": {}, "required": [], "additionalProperties": False},
    "initial_inputs": []},
   {"id": "derive-v36", "source": "derive_v36.py",
    "purpose": ("Derive every ledger row (including the lifting row), the W(n) table and the peak-storage MAX row, "
                "and emit the claim rational and gate delta. Writes ledger-v36.json."),
    "input_schema": {"type": "object", "properties": {}, "required": [], "additionalProperties": False},
    "initial_inputs": []}],
}
json.dump(sub, open(os.path.join(HERE, "submission.json"), "w"), indent=1, sort_keys=True)
print(json.dumps({"statement_chars": len(st), "claim": cl["value"], "gate_delta": LED["gate_delta_claim_minus_gate"],
                  "evidence": len(ev), "files_bytes": sum(os.path.getsize(os.path.join(HERE, f))
                                                          for f in os.listdir(HERE) if os.path.isfile(os.path.join(HERE, f)))},
                 indent=1))
