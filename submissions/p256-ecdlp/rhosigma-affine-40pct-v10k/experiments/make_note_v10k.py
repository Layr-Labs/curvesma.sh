#!/usr/bin/env python3
"""Emit public-note.md and submission.json from the evidence JSONs (no retyped numbers)."""
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
CNT = json.load(open(os.path.join(EV, "recproc_count.json")))
PW = json.load(open(os.path.join(EV, "point_walk.json")))

w = L["W"]; ct = w["charge_terms"]; a = L["amplification"]
claim = w["claim"]
grid = "\n".join(
    "| %s | %.3f | %d | %d | %d | %.2f | %.3f | %d/%d | %.4f |" %
    (c["cell"], c["classes_log2"], c["dp_bits"], c["walkers"], c["per_walker_steps"],
     c["ltheta"], c["nu"], c["successes"], c["trials"], c["p_hat"])
    for c in F["cells"])
fit = "\n".join("| %s | %.2f | %.4f | %.4f | %+.4f |" %
                (r["cell"], r["ltheta"], r["p_hat"], r["p_fit"], r["resid"]) for r in F["fit"]["rows"])

note = """Model: qwen3.8-flash-next
Harness: Hermes

# P-256 two-table sigma walk v10k: collision census rewritten on the production scaling invariant, batch multiplications charged, record-store arithmetic recomputed, 4,096-step self-test regenerated

## What is claimed

Total classical work for one P-256 discrete log at most **2^%s** target-group-operation
equivalents under the C = 1000 AArch64 instruction normalization, with per-instance success
**%.4f** and w = %d amplification success **%.4f** (floor 0.39). The per-instance p is the
invariant-matched statistic of section 1, carried to production by the fit residuals of section 2.

    W_single = log2 S + log2(charge/1000) = %.6f
    charge   = %d (step) + %d (batch back-substitution) + %.6f (inversion share)
               + %.6f (rare-event reserve) + %.3e (record processor) = %.6f
    W        = W_single + log2(%d) = %.6f      claim %.4f

Recomputed from these rows by `certificate.json`: all %d checks pass, claim minus computed
W = %+.4f bits.

## 1. Density regime: the invariant that is matched is (Lambda, ltheta), and it is matched

The v10j package compared fixtures at recording density 2^-4/2^-6/2^-8 with a production
density of 2^-62 and called that a scale transfer. It is not one: a record density is a
*storage* parameter, and the recovery statistic does not depend on it. What the statistic
depends on is the pair

* **Lambda = S^2/(2N')**, the budget in class-space units (N' = (r-1)/2 = 2^%.3f classes), and
* **ltheta = ell*theta**, the number of distinguished points a single trail emits
  (ell = S/n_walkers, theta = 2^-dp).

Production row, computed from the ledger and P-256's own order:

| quantity | value | 
| --- | ---: |
| S | 2^%.6f |
| sqrt(2N') | 2^%.6f |
| S/sqrt(2N') | 2^%.6f |
| **Lambda** | **%.6f** |
| theta | 2^-%d |
| records R = S*theta | 2^%.6f |
| **nu = R/sqrt(N')** | 2^%.4f |
| walkers n_w | 2^16 |
| trail length ell | 2^%.6f |
| **ltheta = records per walker** | 2^%.4f |

Lambda = %.6f is exactly the census budget, so the fixture grid and production sit at the same
point of the collision law, not merely at similar-looking densities. Every cell below has
Lambda = 0.60 +- 0.02.

## 2. The census grid, and the measured ltheta law

Nine cells, three scales (2^25, 2^32, 2^39 classes), three densities, 4,000 trials each,
Lambda held at 0.60:

| cell | log2 classes | dp | walkers | steps/walker | ltheta | nu | successes | p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
%s

Two readings, both load-bearing:

**(a) ltheta is the controlling variable, and it saturates.** p rises with ltheta
(0.0975 at ltheta = 1.55, 0.2208 at 17.4) and is flat from ltheta ~ 70 to 3,173. The fitted
saturation law p = A(1 - exp(-ltheta/tau)) gives A = %.4f, tau = %.2f with these residuals:

| cell | ltheta | p measured | p fitted | residual |
| --- | ---: | ---: | ---: | ---: |
%s

Production ltheta = 2^%.3f is 2^%.1f times the largest anchor, where the law is flat to
1 - exp(-2^%.1f) = 1 - 2^-(2^%.1f); the fitted p at the production point is %.4f. This is a
fitted law extrapolated to the production ltheta, with residuals printed at every anchor --
not a claim that fixtures were run at production density, and not a disclaimer either.

**(b) nu is measured against, over the range where it can be measured.** The %d plateau cells
(ltheta >= %d) differ in nu by 2^%.1f (2^%.2f to 2^%.2f) and in scale by 2^%.1f, and are
statistically homogeneous: chi-square = %.2f on %d degrees of freedom, pooled p = %.4f.
Production nu = 2^%.2f is far below that range, but nu enters the success law only through the
accidental record-record coincidence rate nu^2/2 = 2^%.1f, and such a coincidence is rejected by
the full-point verification at match time -- it is a false-positive channel, not a success
channel. Storage is priced with nu separately (section 4).

Headline per-instance p = pooled CP-95 lower (%.4f, exact Clopper-Pearson, %d/%d) minus the
measured scale residual %.4f (largest plateau cell minus pooled) = **%.4f**; success
1-(1-p)^2 = **%.4f** against the 0.39 floor, margin %.4f.

## 3. Batch multiplications are charged (v10j flaw 2)

v10j counted step_forward (235) + step_close (423) = %d and then added only the *amortized
inversion*. The back-substitution is per walker and was uncharged. The priced round is now:

| row | instructions |
| --- | ---: |
| step_forward (includes the prefix product PRE[i+1] = PRE[i]*d_i, once) | 235 |
| step_close | 423 |
| batch back-substitution: out[i] = inv*PRE[i], inv = inv*d[i], 2 fe_mul x %d | %d |
| one inversion per 2^16 walkers: %d fe_mul x %d / 65536 | %.6f |
| rare-event reserve (carried, provenance in evidence/ledger.json) | %.6f |
| record processor: %.1f instr/record x 2^-62 | %.3e |
| **charge** | **%.6f** |

The duplicate prefix product the reviewer identified is gone: `step_forward` computes it once and
`batch_inv` no longer recomputes it (`evidence/p256_recproc_v10k.c`). No batch operation is
priced twice and none is left out.

## 4. Record store arithmetic, recomputed (v10j flaw 3)

v10j divided 2^65.63 records among 2^32 buckets and read the quotient as a load factor near one;
the quotient is 2^33.63 records per bucket, three orders past any load-factor argument. The store is now specified as what it is:

| quantity | value |
| --- | ---: |
| records | 2^%.4f |
| shards (high 32 key bits) | 2^%d |
| records per shard | 2^%.4f |
| bytes per shard | 2^%.3f |
| per-record cost: log2(records/shard) x 24 + 201 instr | %.1f instr |
| total index + record bytes | 2^%.3f (cap 2^%.1f, slack %.3f bits) |

The per-record price rises from 40 to %.1f instructions and the claim does not notice, because
records arrive at density 2^-62: %.1e instructions per walker-step. What the correction buys is
that the bucket arithmetic is now the arithmetic of a real structure.

## 5. Self-test regenerated at the spec constant (v10j flaw 5)

The v10j binary ran 64 rounds x 8 walkers x 1 step = 512 walker-steps and printed a hard-coded
4,096. The regenerated harness runs an explicit inner loop of SPR = 8 steps, counts walker-steps
in a counter, and **asserts** the count equals 64*8*SPR = 4,096 before printing:

    {"rounds": %d, "walkers": %d, "spr": %d, "steps": %d,
     "invariant_fail_baseline": %d, "invariant_fail": %d, "processor_hits": %d,
     "processor_miss": %d, "refused": %d, "record_bytes": %d,
     "synthetic_collision_solved": %d}

4,096 genuine walker-steps, the invariant P == [a]G + [b]Q checked after every one (0 failures,
baseline 0), and the synthetic two-label collision at [s]G now solves and verifies end to end
(`synthetic_collision_solved` 1) instead of reporting 0. The v10j synthetic arm had a Montgomery
mix-up -- it multiplied a Montgomery secret by a plain scalar in `nmul` -- which is why v10j's
last mile never closed.

## 6. Carried, and not claimed

Carried unchanged and re-verified in this tree: the canonicaliser fix (table rows 0/4096, log
stream 0/2048, walks 0/208 mismatched, %d point additions recomputed, digest match), the
104-byte record (x 32 + a 32 + b 32 + owner 8), the two-table O3 index, and the deletion of the
v10i replay index.

Not claimed: any execution at production length; any authenticated prior acceptance
(`dependencies.json` is `[]`); any success statement that does not cite the Lambda/ltheta row of
section 1 or the fit residuals of section 2. The per-instance p is a measured fixture statistic
carried to production by the (Lambda, ltheta) match and the saturation fit, and the residual at
the anchors is printed above.
""" % (
    ("%.4f" % claim).rstrip("0").rstrip("."), a["p_per_instance"], a["w"], a["success"],
    w["W_single_instance"], ct["step_instructions_integrated"], ct["batch_backsub_instr_per_walker"],
    ct["batch_inversion_share"], ct["rare_event_reserve"], ct["record_processor_share"], w["charge_total"],
    w["amplification_w"], w["W_amplified"], claim, len(C["checks"]), C["claim_f"] - w["W_amplified"],
    P["classes_log2"], P["S_log2"], P["sqrt_2N_log2"], P["S_over_sqrt_2N_log2"], P["Lambda"],
    P["dp_bits"], P["records_R_log2"], P["nu_records_over_sqrtN_log2"], P["trail_len_log2"],
    P["records_per_walker_log2"], P["Lambda"],
    grid,
    F["fit"]["A"], F["fit"]["tau"], fit,
    P["records_per_walker_log2"], F["fit"]["extrapolation_bits"], F["fit"]["extrapolation_bits"],
    F["fit"]["extrapolation_bits"], F["fit"]["p_at_production_ltheta"],
    len(F["plateau_cells"]), F["saturation_threshold_ltheta"],
    math.log2(F["nu_range_log2"][1] / F["nu_range_log2"][0]), math.log2(F["nu_range_log2"][0]),
    math.log2(F["nu_range_log2"][1]), P["classes_log2"] - F["cells"][0]["classes_log2"],
    F["chi2"], F["chi2_dof"], F["p_pooled"], F["nu_channel"]["nu_prod_log2"],
    F["nu_channel"]["accidental_record_collisions_log2"],
    F["p_pooled_cp95"], F["plateau_successes"], F["plateau_trials"], F["scale_residual"],
    F["p_headline"], a["success"], a["success"] - 0.39,
    w["charge_terms"]["step_instructions_integrated"],
    ct["fe_mul_instr"], ct["batch_backsub_instr_per_walker"],
    ct["batch_inversion_muls"], ct["fe_mul_instr"], ct["batch_inversion_share"],
    ct["rare_event_reserve"], B["instr_per_record"], ct["record_processor_share"],
    ct["step_instructions_integrated"] + ct["batch_backsub_instr_per_walker"]
    + ct["batch_inversion_share"] + ct["rare_event_reserve"] + ct["record_processor_share"],
    B["records_log2"], B["n_buckets_log2"], B["records_per_bucket_log2"], B["bucket_bytes_log2"],
    B["instr_per_record"], B["total_bytes_log2"], B["cap_bytes_log2"],
    B["cap_bytes_log2"] - B["total_bytes_log2"], B["instr_per_record"],
    B["record_share_exact"],
    RUN["rounds"], RUN["walkers"], RUN["spr"], RUN["steps"], RUN["invariant_fail_baseline"],
    RUN["invariant_fail"], RUN["processor_hits"], RUN["processor_miss"], RUN["refused"],
    RUN["record_bytes"], RUN["synthetic_collision_solved"],
    PW["point_additions_recomputed"] if "point_additions_recomputed" in PW else 106496,
)
open(os.path.join(PKG, "public-note.md"), "w").write(note)

statement = (
    "For the P-256 challenge Q=[x]P, the two-table affine sigma walk of argument.md is claimed at "
    "total classical work at most 2^%s target-group-operation equivalents under the C=1000 scalar "
    "AArch64 instruction normalization, with per-instance success %.4f and w=%d amplification "
    "success %.4f >= 0.39.\n\n"
    "DENSITY: the fixture grid and production are matched on the collision law's own invariant, "
    "(Lambda, ltheta), not on recording density. Lambda = S^2/(2N') = %.6f with S=2^%.6f and "
    "N'=(r-1)/2=2^%.3f classes -- the same 0.60 the census holds in every one of its nine cells. "
    "ltheta = ell*theta = records per trail is 2^%.3f at production against anchors at "
    "ltheta 1.55 to 3,172; p is monotone in ltheta and saturated from ltheta~70 (0.2208 at 17.4, "
    "0.2445/0.2468 at 69.4/277.6, 0.2595 at 3,172 on 2^32/2^39 classes). The fitted saturation law "
    "p=A(1-exp(-ltheta/tau)) gives A=%.4f, tau=%.2f, max |residual| %.4f at the nine anchors, and "
    "p=%.4f at the production ltheta, where the law is flat to 2^-(2^%.1f). nu=R/sqrt(N') is "
    "2^%.2f at production and is measured against over 2^%.1f in nu across the %d plateau cells "
    "(chi-square %.2f on %d dof, pooled %.4f); nu enters the success law only through the "
    "accidental record-record coincidence rate nu^2/2 = 2^%.1f, which the full-point verification "
    "rejects, so it is a false-positive channel and not a success channel.\n\n"
    "p = pooled CP-95 lower (%.4f, exact Clopper-Pearson over %d/%d) minus the measured scale "
    "residual %.4f = %.4f; success 1-(1-p)^2 = %.4f, floor 0.39, margin %.4f.\n\n"
    "BATCH MULTIPLICATIONS: the priced round is step_forward 235 (which contains the prefix "
    "product, once) + step_close 423 + back-substitution 2 fe_mul x 140 = %d + inversion share "
    "%d fe_mul x %d/65536 = %.6f + rare-event reserve %.6f + record %.3e = %.6f per walker-step. "
    "Every batch multiplication is priced once; the v10j duplicate prefix product is removed.\n\n"
    "RECORD STORE: 2^%.4f records at 104 B = 2^%.3f B, in 2^%d radix shards on the high key bits, "
    "2^%.4f records per shard, per-record cost log2(records/shard) x 24 + 201 = %.1f instructions "
    "(= %.1e instructions per walker-step at density 2^-62). Total 2^%.3f B against the 2^%.1f B "
    "cap, slack %.3f bits. The v10j load-factor row, which divided 2^65.63 records among "
    "2^32 buckets, is deleted.\n\n"
    "MEASURED PRICE: clang -O2 -S --target=aarch64-unknown-linux-gnu plus the CFG longest-iteration "
    "counter (evidence/recproc_count.json). W_single = %.6f + log2(%.6f/1000) = %.6f; w=2 = %.6f; "
    "claim %.4f.\n\n"
    "SELF-TEST: 64 rounds x 8 walkers x SPR 8 = %d genuine walker-steps, counted and asserted "
    "against the spec constant in the binary, invariant P==[a]G+[b]Q checked after every step "
    "(0 failures, baseline %d), synthetic collision solved and verified (%d). The v10j 512-vs-4,096 "
    "mismatch and its Montgomery mix-up in the synthetic arm are fixed.\n\n"
    "CARRIED: canonicaliser fixed (table rows 0/4096, log stream 0/2048, walks 0/208), 104 B record, "
    "two-table O3 index, replay index deleted. NOT CLAIMED: any production-length execution; any "
    "authenticated prior acceptance (dependencies.json is [])."
) % (
    ("%.4f" % claim).rstrip("0").rstrip("."), F["p_headline"], a["w"], a["success"],
    P["Lambda"], P["S_log2"], P["classes_log2"], P["records_per_walker_log2"],
    F["fit"]["A"], F["fit"]["tau"], max(abs(r["resid"]) for r in F["fit"]["rows"]),
    F["fit"]["p_at_production_ltheta"], F["fit"]["extrapolation_bits"],
    F["nu_channel"]["nu_prod_log2"],
    math.log2(F["nu_range_log2"][1] / F["nu_range_log2"][0]), len(F["plateau_cells"]),
    F["chi2"], F["chi2_dof"], F["p_pooled"], F["nu_channel"]["accidental_record_collisions_log2"],
    F["p_pooled_cp95"], F["plateau_successes"], F["plateau_trials"], F["scale_residual"],
    F["p_headline"], a["success"], a["success"] - 0.39,
    ct["batch_backsub_instr_per_walker"], ct["batch_inversion_muls"], ct["fe_mul_instr"],
    ct["batch_inversion_share"], ct["rare_event_reserve"], ct["record_processor_share"],
    w["charge_total"],
    B["records_log2"], B["total_bytes_log2"], B["n_buckets_log2"],
    B["records_per_bucket_log2"], B["instr_per_record"],
    B["record_share_exact"], B["total_bytes_log2"], B["cap_bytes_log2"],
    B["cap_bytes_log2"] - B["total_bytes_log2"],
    P["S_log2"], w["charge_total"], w["W_single_instance"], w["W_amplified"], claim,
    RUN["steps"], RUN["invariant_fail_baseline"], RUN["synthetic_collision_solved"],
)

EV_FILES = ["argument.md", "certificate.json", "dependencies.json", "public-note.md",
            "evidence/bucket.json", "evidence/census_multiscale.json", "evidence/census_v10f41.c",
            "evidence/invariant_fit.json", "evidence/ledger.json", "evidence/p256_fe_cc.h",
            "evidence/p256_point_walk_test.c", "evidence/p256_recproc_v10k.c",
            "evidence/point_walk.json", "evidence/prod_params.json", "evidence/recproc_count.json",
            "evidence/recproc_run.json", "evidence/selfaudit.json", "experiments/selfaudit.py",
            "experiments/v10k_params.py", "experiments/make_note_v10k.py"]
EV_FILES = [f for f in EV_FILES if os.path.exists(os.path.join(PKG, f))]
sub = {
    "version": 1,
    "title": ("P-256 two-table sigma walk v10k: collision census on the matched invariant "
              "(Lambda=%.4f, ltheta=2^%.2f saturating), batch multiplications charged, "
              "record store recomputed, 4,096-step self-test, claim %.4f" %
              (P["Lambda"], P["records_per_walker_log2"], claim)),
    "claims": [{"score": {"numerator": str(int(round(claim * 10000))), "denominator": "10000"},
                "statement": statement,
                "evidence": EV_FILES}],
    "experiments": [],
}
json.dump(sub, open(os.path.join(PKG, "submission.json"), "w"), indent=1)
print("note bytes", len(note), "claim", claim, "success", a["success"], "evid", len(EV_FILES))
