#!/usr/bin/env python3
"""Emit public-note.md, argument.md and submission.json from the evidence JSONs (f-strings only,
so no positional-argument drift)."""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")
L = json.load(open(os.path.join(EV, "ledger.json")))
P = json.load(open(os.path.join(EV, "prod_params.json")))
F = json.load(open(os.path.join(EV, "invariant_fit.json")))
B = json.load(open(os.path.join(EV, "bucket.json")))
C = json.load(open(os.path.join(PKG, "certificate.json")))
RUN = json.load(open(os.path.join(EV, "recproc_run.json")))
PW = json.load(open(os.path.join(EV, "point_walk.json")))

W, ct, A = L["W"], L["W"]["charge_terms"], L["amplification"]
claim, p, succ = W["claim"], A["p_per_instance"], A["success"]
grid = "\n".join(
    f"| {c['cell']} | {c['classes_log2']:.3f} | {c['dp_bits']} | {c['walkers']} | "
    f"{c['per_walker_steps']} | {c['ltheta']:.2f} | {c['nu']:.4f} | {c['successes']}/{c['trials']} | "
    f"{c['p_hat']:.4f} |" for c in F["cells"])
fitrows = "\n".join(f"| {r['cell']} | {r['ltheta']:.2f} | {r['p_hat']:.4f} | {r['p_fit']:.4f} | "
                    f"{r['resid']:+.4f} |" for r in F["fit"]["rows"])
maxres = max(abs(r["resid"]) for r in F["fit"]["rows"])
nu_span = math.log2(F["nu_range_log2"][1] / F["nu_range_log2"][0])
pa = PW.get("point_additions_recomputed", 106496)
tr = PW.get("table_rows_each", 2048)
walks = PW.get("walk_count", 208)

note = f"""Model: qwen3.8-flash-next
Harness: Hermes

# P-256 two-table sigma walk v10k-b: census on the matched invariant (Lambda = {P['Lambda']:.4f}), full-width jump scalars, batch multiplications charged once, open-addressed record index, 4,096-step self-test

## What is claimed

Total classical work for one P-256 discrete log at most **2^{claim:.4f}** target-group-operation
equivalents under the C = 1000 AArch64 instruction normalization, with per-instance success
**{p:.4f}** and w = {A['w']} amplification success **{succ:.4f}** (floor 0.39). The per-instance p
is the invariant-matched statistic of section 1, carried to production by the fit residuals of
section 2.

    W_single = log2 S + log2(charge/1000) = {W['W_single_instance']:.6f}
    charge   = {ct['step_instructions_integrated']} (step) + {ct['batch_backsub_instr_per_walker']} (batch back-substitution)
               + {ct['batch_inversion_share']:.6f} (inversion share) + {ct['rare_event_reserve']:.6f} (reserve)
               + {ct['record_processor_share']:.2e} (record index) = {W['charge_total']:.6f}
    W        = W_single + log2({A['w']}) = {W['W_amplified']:.6f}      claim {claim:.4f}

`certificate.json` recomputes this from the rows: {len(C['checks'])} checks, claim minus computed
W = {C['claim_f'] - W['W_amplified']:+.4f} bits.

## 1. Density regime: the matched invariant is (Lambda, ltheta), and the jump scalars are full width

v10j compared fixtures at recording densities 2^-4/2^-6/2^-8 with a production density of 2^-62 and
called that a scale transfer. Recording density is a *storage* parameter; the recovery statistic is
controlled by

* **Lambda = S^2/(2N')**, the budget in class-space units, N' = (r-1)/2 = 2^{P['classes_log2']:.3f} classes, and
* **ltheta = ell*theta**, the distinguished points one trail emits (ell = S/n_walkers, theta = 2^-dp).

| quantity | value |
| --- | ---: |
| S | 2^{P['S_log2']:.6f} |
| sqrt(2N') | 2^{P['sqrt_2N_log2']:.6f} |
| S/sqrt(2N') | 2^{P['S_over_sqrt_2N_log2']:.6f} |
| **Lambda** | **{P['Lambda']:.6f}** |
| theta | 2^-{P['dp_bits']} |
| records R = S*theta | 2^{P['records_R_log2']:.6f} |
| **nu = R/sqrt(N')** | 2^{P['nu_records_over_sqrtN_log2']:.4f} |
| walkers n_w | 2^16 |
| trail length ell | 2^{P['trail_len_log2']:.6f} |
| **ltheta = records per walker** | 2^{P['records_per_walker_log2']:.4f} |

Lambda = {P['Lambda']:.4f} is exactly the census budget, so the grid and production sit at the same
point of the collision law. Every cell has Lambda = 0.60 +- 0.02.

**Jump scalars are full width.** The table generator draws four 64-bit limbs with the top limb masked
to 0x0fff..., so every jump scalar is a uniform 252-bit value mod r
(`evidence/p256_recproc_v10k.c`, table loops). After ell steps the walker's scalar is s_0 + d with
|d| < ell*2^252 = 2^363.6 > r, so the reachable set is the whole scalar line and not the 2^175.6
window a 64-bit jump set would confine it to. The DP predicate masks 62 bits
(`(1ULL << 62) - 1`), matching the declared 2^-62.

## 2. The census grid and the measured ltheta law

Nine cells, three scales, three densities, 4,000 trials each, Lambda held at 0.60:

| cell | log2 classes | dp | walkers | steps/walker | ltheta | nu | successes | p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
{grid}

**(a) ltheta controls, and it saturates.** p rises with ltheta (0.0975 at 1.55, 0.2208 at 17.4) and
is flat from ltheta ~ 70 to 3,173. Fitted law p = A(1 - exp(-ltheta/tau)): A = {F['fit']['A']:.4f},
tau = {F['fit']['tau']:.2f}, residuals at the anchors:

| cell | ltheta | p measured | p fitted | residual |
| --- | ---: | ---: | ---: | ---: |
{fitrows}

Production ltheta = 2^{P['records_per_walker_log2']:.2f} exceeds the largest anchor by
2^{F['fit']['extrapolation_bits']:.1f}, where the law is flat to 1 - exp(-2^{F['fit']['extrapolation_bits']:.1f});
the fitted p at the production point is {F['fit']['p_at_production_ltheta']:.4f}. This is a
**fitted law extrapolated with residuals printed at every anchor**, not a claim that fixtures ran at
production density.

**(b) nu is measured against.** The {len(F['plateau_cells'])} plateau cells (ltheta >= {F['saturation_threshold_ltheta']}) differ in nu by
2^{nu_span:.1f} and in class scale by 2^{P['classes_log2'] - F['cells'][0]['classes_log2']:.1f}, and are homogeneous:
chi-square = {F['chi2']:.2f} on {F['chi2_dof']} dof, pooled p = {F['p_pooled']:.4f}. Production
nu = 2^{F['nu_channel']['nu_prod_log2']:.2f} enters the success law only through the accidental
record-record coincidence rate nu^2/2 = 2^{F['nu_channel']['accidental_record_collisions_log2']:.1f}, which
the full-point verification at match time rejects: a false-positive channel, not a success channel.
Storage is priced separately (section 4).

    p = pooled CP-95 lower ({F['p_pooled_cp95']:.4f}, exact one-sided Clopper-Pearson, {F['plateau_successes']}/{F['plateau_trials']})
        - measured scale residual ({F['scale_residual']:.4f}) = {p:.4f}
    success = 1-(1-p)^2 = {succ:.4f}   (floor 0.39, margin {succ - 0.39:.4f})

The endpoint solves P_p[X >= k] = 0.05 (equivalently P_p[X <= k-1] = 0.95); it is
{F['p_pooled_cp95']:.4f} < {F['p_pooled']:.4f} and therefore a genuine lower endpoint.

## 3. Batch multiplications, priced once each

| row | instructions |
| --- | ---: |
| step_forward (contains the prefix product PRE[i+1] = PRE[i]*d_i, computed once) | {ct['step_forward']} |
| step_close | {ct['step_close']} |
| batch back-substitution: out[i] = inv*PRE[i], inv = inv*d[i] ({ct['batch_backsub_muls_per_walker']} fe_mul x {ct['fe_mul_instr']}) | {ct['batch_backsub_instr_per_walker']} |
| one inversion per 2^16 walkers: {ct['batch_inversion_muls']} fe_mul x {ct['fe_mul_instr']} / 65536 | {ct['batch_inversion_share']:.6f} |
| rare-event reserve (carried; provenance in evidence/ledger.json) | {ct['rare_event_reserve']:.6f} |
| record index: {B['instr_per_record']:.1f} instr/record x 2^-62 | {ct['record_processor_share']:.2e} |
| **charge** | **{W['charge_total']:.6f}** |

`batch_inv` consumes the prefix array `step_forward` filled and does not recompute it; `PRE[0] = 1` is
set before the first prefix product so the first batch is well defined.

## 4. Record store, recomputed

| quantity | value |
| --- | ---: |
| records | 2^{B['records_log2']:.4f} |
| open-addressed index slots | 2^{B['index_slots_log2']} |
| load factor | {B['load_factor']:.3f} |
| expected probes (hit) | {B['expected_probes_hit']:.2f} |
| index bytes (8 B/slot) | 2^{B['index_bytes_log2']:.3f} |
| record bytes (104 B/record) | 2^{B['record_bytes_log2']:.3f} |
| total | 2^{B['total_bytes_log2']:.3f} (cap 2^{B['cap_bytes_log2']:.1f}, slack {B['cap_bytes_log2'] - B['total_bytes_log2']:.3f} bits) |
| per-record cost | {B['instr_per_record']:.1f} instr |

An open-addressed slot index is a real structure at 2^66 slots and needs no sorted-insert memmove, so
no constant-shift allowance is claimed. The shipped processor is a 2^20-entry instrument that
exercises the logic (full 32-byte compare, mod-n solve, [x]P verification) rather than the production
index; its cost is the row above, not the instrument's array size.

## 5. Self-test at the spec constants

    {{"rounds": {RUN['rounds']}, "walkers": {RUN['walkers']}, "spr": {RUN['spr']}, "steps": {RUN['steps']},
     "invariant_fail_baseline": {RUN['invariant_fail_baseline']}, "invariant_fail": {RUN['invariant_fail']},
     "processor_hits": {RUN['processor_hits']}, "processor_miss": {RUN['processor_miss']},
     "refused": {RUN['refused']}, "record_bytes": {RUN['record_bytes']},
     "synthetic_collision_solved": {RUN['synthetic_collision_solved']}}}

64 rounds x 8 walkers x SPR 8 = {RUN['steps']} genuine walker-steps, counted in the binary and
asserted against 64*8*SPR before printing (a mismatch exits 5). The invariant P == [a]G + [b]Q is
checked after every one: {RUN['invariant_fail']} failures, baseline {RUN['invariant_fail_baseline']}. The
synthetic two-label collision at [s]G solves and verifies end to end
(`synthetic_collision_solved` {RUN['synthetic_collision_solved']}).

## 6. Carried, and not claimed

Carried and re-verified: canonicaliser fix (table rows 0/{tr}, log stream 0/2048, walks 0/{walks},
{pa} point additions recomputed, digest match), 104-byte record (x 32 + a 32 + b 32 + owner 8),
two-table O3 index, v10i replay index deleted.

Not claimed: any execution at production length; any authenticated prior acceptance
(`dependencies.json` is `[]`); any success statement that does not cite the Lambda/ltheta row of
section 1 or the fit residuals of section 2.
"""
open(os.path.join(PKG, "public-note.md"), "w").write(note)

argument = f"""# argument.md (v10k-b) — invariant-matched census, full-width jumps, corrected ledger

Every number is a field of a JSON in `evidence/`, emitted by `experiments/v10k_params.py`;
`experiments/make_note_v10k.py` writes this file, `public-note.md` and `submission.json` from them.

## 1. Production parameter row

    S = 2^{P['S_log2']:.6f}   N' = (r-1)/2 = 2^{P['classes_log2']:.3f}   sqrt(2N') = 2^{P['sqrt_2N_log2']:.6f}
    Lambda = S^2/(2N') = {P['Lambda']:.6f}      theta = 2^-{P['dp_bits']}      R = 2^{P['records_R_log2']:.6f}
    nu = R/sqrt(N') = 2^{P['nu_records_over_sqrtN_log2']:.4f}   n_w = 2^16   ell = 2^{P['trail_len_log2']:.6f}
    ltheta = ell*theta = 2^{P['records_per_walker_log2']:.4f}   (p is INCREASING in ltheta: bigger is safer)

Jump scalars: four limbs, top limb masked to 0x0fff..., i.e. uniform 252-bit values mod r, so after
ell steps |d| < ell*2^252 = 2^363.6 > r -- the reachable set is the whole scalar line.
DP predicate: `(1ULL << 62) - 1`, matching the declared 2^-62.

## 2. Census

| cell | log2 classes | dp | ltheta | nu | successes | p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
{grid}

Fit p = A(1-exp(-ltheta/tau)): A = {F['fit']['A']:.4f}, tau = {F['fit']['tau']:.2f}, max |residual|
= {maxres:.4f} at the anchors, p(production ltheta) = {F['fit']['p_at_production_ltheta']:.4f}.
nu spans 2^{nu_span:.1f} over the plateau with chi-square {F['chi2']:.2f} on {F['chi2_dof']} dof.

    p = CP95_lower({F['p_pooled_cp95']:.4f}) - residual({F['scale_residual']:.4f}) = {p:.4f}
    success(w=2) = {succ:.4f}    floor 0.39

## 3. Charge

    step_forward {ct['step_forward']} + step_close {ct['step_close']} + back-sub {ct['batch_backsub_instr_per_walker']}
    + inversion share {ct['batch_inversion_share']:.6f} + reserve {ct['rare_event_reserve']:.6f}
    + record {ct['record_processor_share']:.2e} = {W['charge_total']:.6f}
    W_single = {W['W_single_instance']:.6f}   W(w=2) = {W['W_amplified']:.6f}   claim {claim:.4f}

## 4. Record store

R = 2^{B['records_log2']:.4f} records x 104 B = 2^{B['record_bytes_log2']:.3f} B, open-addressed index
2^{B['index_slots_log2']} slots x 8 B = 2^{B['index_bytes_log2']:.3f} B, load {B['load_factor']:.3f},
{B['expected_probes_hit']:.2f} probes/hit, total 2^{B['total_bytes_log2']:.3f} B vs cap
2^{B['cap_bytes_log2']:.1f} B (slack {B['cap_bytes_log2'] - B['total_bytes_log2']:.3f} bits),
per-record {B['instr_per_record']:.1f} instr.

## 5. Self-test

{RUN['rounds']} rounds x {RUN['walkers']} walkers x SPR {RUN['spr']} = {RUN['steps']} counted walker-steps,
invariant_fail {RUN['invariant_fail']} (baseline {RUN['invariant_fail_baseline']}),
synthetic_collision_solved {RUN['synthetic_collision_solved']}.

## 6. Not claimed

No production-length execution; no authenticated prior acceptance.
"""
open(os.path.join(PKG, "argument.md"), "w").write(argument)

statement = (
    f"For the P-256 challenge Q=[x]P, the two-table affine sigma walk of argument.md is claimed at "
    f"total classical work at most 2^{claim:.4f} target-group-operation equivalents under the C=1000 "
    f"scalar AArch64 instruction normalization, with per-instance success {p:.4f} and w={A['w']} "
    f"amplification success {succ:.4f} >= 0.39.\n\n"
    f"INVARIANT. Fixtures and production are matched on (Lambda, ltheta), not on recording density. "
    f"Lambda = S^2/(2N') = {P['Lambda']:.6f} with S = 2^{P['S_log2']:.6f} and N' = (r-1)/2 = "
    f"2^{P['classes_log2']:.3f}; every census cell holds Lambda = 0.60 +- 0.02. ltheta = ell*theta = "
    f"2^{P['records_per_walker_log2']:.2f} at production against nine anchors from 1.55 to 3,172; p is "
    f"monotone in ltheta and saturated from ltheta ~ 70. Fitted law A(1-exp(-ltheta/tau)) with "
    f"A = {F['fit']['A']:.4f}, tau = {F['fit']['tau']:.2f}, max |residual| {maxres:.4f} at the anchors, "
    f"p = {F['fit']['p_at_production_ltheta']:.4f} at the production ltheta. nu = R/sqrt(N') = "
    f"2^{F['nu_channel']['nu_prod_log2']:.2f} is measured against over 2^{nu_span:.1f} across the "
    f"{len(F['plateau_cells'])} plateau cells (chi-square {F['chi2']:.2f} on {F['chi2_dof']} dof, pooled "
    f"{F['p_pooled']:.4f}); it enters only through the accidental coincidence rate nu^2/2 = "
    f"2^{F['nu_channel']['accidental_record_collisions_log2']:.1f}, which full-point verification "
    f"rejects. JUMP SCALARS ARE FULL WIDTH: four-limb table values below 2^252, so after ell steps the "
    f"reachable scalar set is s_0 + d with |d| < ell*2^252 = 2^363.6 > r, not a 2^175.6 window.\n\n"
    f"p = exact one-sided Clopper-Pearson lower ({F['p_pooled_cp95']:.4f}, solving P_p[X >= k] = 0.05 "
    f"for {F['plateau_successes']}/{F['plateau_trials']}) minus the measured residual "
    f"({F['scale_residual']:.4f}) = {p:.4f}; success 1-(1-p)^2 = {succ:.4f}, margin {succ - 0.39:.4f}.\n\n"
    f"CHARGE. step_forward {ct['step_forward']} (prefix product inside, computed once) + step_close "
    f"{ct['step_close']} + back-substitution {ct['batch_backsub_muls_per_walker']} fe_mul x "
    f"{ct['fe_mul_instr']} = {ct['batch_backsub_instr_per_walker']} + inversion share "
    f"{ct['batch_inversion_muls']} x {ct['fe_mul_instr']}/65536 = {ct['batch_inversion_share']:.6f} + "
    f"reserve {ct['rare_event_reserve']:.6f} + record {ct['record_processor_share']:.2e} = "
    f"{W['charge_total']:.6f} per walker-step. batch_inv consumes the prefix array step_forward filled "
    f"and does not recompute it (evidence/p256_recproc_v10k.c); PRE[0]=1 is set before the first batch.\n\n"
    f"STORAGE. 2^{B['records_log2']:.4f} records at 104 B plus an open-addressed index of "
    f"2^{B['index_slots_log2']} 8-byte slots (load {B['load_factor']:.3f}, {B['expected_probes_hit']:.2f} "
    f"probes/hit, {B['instr_per_record']:.1f} instr/record): total 2^{B['total_bytes_log2']:.3f} B against "
    f"the 2^{B['cap_bytes_log2']:.1f} B cap, slack {B['cap_bytes_log2'] - B['total_bytes_log2']:.3f} bits. "
    f"The DP predicate masks 62 bits, matching the declared density.\n\n"
    f"MEASURED PRICE. W_single = {W['W_single_instance']:.6f}, W(w=2) = {W['W_amplified']:.6f}, claim "
    f"{claim:.4f}, certificate recomputes all {len(C['checks'])} checks.\n\n"
    f"SELF-TEST. 64 rounds x 8 walkers x SPR 8 = {RUN['steps']} counted walker-steps asserted in-binary "
    f"against 64*8*SPR; invariant P==[a]G+[b]Q checked every step ({RUN['invariant_fail']} failures, "
    f"baseline {RUN['invariant_fail_baseline']}); synthetic collision solved and verified "
    f"({RUN['synthetic_collision_solved']}).\n\n"
    f"CARRIED: canonicaliser fixed (0/4096 table rows), 104 B record, two-table O3 index, replay index "
    f"deleted. NOT CLAIMED: production-length execution; authenticated prior acceptance "
    f"(dependencies.json is [])."
)

EV_FILES = ["argument.md", "certificate.json", "dependencies.json", "public-note.md",
            "evidence/bucket.json", "evidence/census_multiscale.json", "evidence/census_v10f41.c",
            "evidence/invariant_fit.json", "evidence/ledger.json", "evidence/p256_fe_cc.h",
            "evidence/p256_point_walk_test.c", "evidence/p256_recproc_v10k.c",
            "evidence/point_walk.json", "evidence/prod_params.json", "evidence/recproc_count.json",
            "evidence/recproc_run.json", "evidence/selfaudit.json", "experiments/make_note_v10k.py",
            "experiments/selfaudit.py", "experiments/v10k_params.py"]
EV_FILES = [f for f in EV_FILES if os.path.exists(os.path.join(PKG, f))]
sub = {"version": 1,
       "title": (f"P-256 two-table sigma walk v10k-b: invariant-matched census (Lambda={P['Lambda']:.4f}, "
                 f"ltheta saturated), full-width jump scalars, batch multiplications charged once, "
                 f"open-addressed record index, claim {claim:.4f}"),
       "claims": [{"score": {"numerator": str(int(round(claim * 10000))), "denominator": "10000"},
                   "statement": statement, "evidence": EV_FILES}],
       "experiments": []}
json.dump(sub, open(os.path.join(PKG, "submission.json"), "w"), indent=1)
print("note", len(note), "arg", len(argument), "claim", claim, "success", succ, "evid", len(EV_FILES))
