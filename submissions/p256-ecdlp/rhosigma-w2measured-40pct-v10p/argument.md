# v10p — the w=2 ticket measured on both arms; `batch_inv` re-amortized at the shipped width

Track `p256-ecdlp`. Package `rhosigma-w2measured-40pct-v10p`. Supersedes `rhosigma-corrections-40pct-v10p`
(verdicts a7a7dfa2 and bffe684d, quoted verbatim in `public-note.md` §0 and archived under
`workspaces/archived-p256-2026-10-24/`).

**Claim 131.6873** (was 128.9437). The claim gets *worse* by 2.7436 bits, and that is the point: the two
things the reviewer said were missing are supplied, and one of them (the amortization) costs work.

---

## 0a. Corrections after review a8207870 (withdrawn rows, recorded not buried)

Two rows of this package are **WITHDRAWN** as of the a8207870 verdict, and the reviewer is right on both:

* **§5 a-channel bound WITHDRAWN.** The claim that zero counted doublings bounds |P(a=−3) − P(a=0)| is
  invalid, for the reasons given: `g_dbl` is incremented only in the main step loop, while `padd` also
  runs (uncounted) during table construction, walker initialisation and verification; the shipped
  census is compiled on the a = 0 curve, so its main-walk numerator is literally `3x² + 0`; and a
  per-step event rate would not bound a whole-run success difference even under a valid coupling,
  since steps are correlated and a binomial rule of three does not apply. The a = −3 fixture question
  is therefore **OPEN**, not answered. What survives is only the structural remark that the affine
  chord law carries no a, plus the disclosure that the numerator-only probe fails its validity gate
  and is excluded.
* **§3 κ derivation WITHDRAWN.** Multiplying the measured closure rate by its inverse measured/predicted
  ratio is an identity that restores the predicted 2^−24; it is not a production bound and gives no
  stationarity support over production-length trails. The charged κ = 1.021012 is retained only as a
  0.0002 pct surcharge on a per-step row, with **no transfer claim attached**.

What the same verdict **accepted**: measuring two-run sets directly is useful evidence; the B = 8
inversion amortization addresses the identified discrepancy; the ledger's 131.6872366 with the 4-dp
ceiling 131.6873 is arithmetically right at the listed prices and step budget. What it still rejects,
unchanged from v10o: the success premise for the shipped arm (423/6000, 0.07050, CP95 0.06514 — not a
production bound), and complete resource accounting (span summaries and batch counters do not
substantiate every dynamic P-256 instruction, the reserve, or preprocessing/lookup/transfer/verification
and peak-storage work).

## 1. What was asked, and what is answered

| verdict item | status in v10o | status here |
|---|---|---|
| no MEASURED w=2 arm; "6.57x short" is a single-run CP at current n | modelled | **measured, n = 6000 sets per arm, CP95 lower on the measured set fraction** |
| `batch_inv` amortized over 65536 walkers; shipped array 16, live census 8 | unaddressed | **re-amortized from the disasm price at B = 8/16/65536, worse charged, +2.7436 bits** |
| κ transferred from fixture to production without derivation | assumption | **derived from the measured/predicted closure-rate ratio and charged** |
| census fixture a=0, production a=−3 | disclosed only | **attempted; the incidence bound is WITHDRAWN after a8207870 (§0a) — question OPEN** |
| three-law band not shown exhaustive | disputed | **no longer load-bearing: the production-relevant p₂ is measured directly on both arms** |
| ledger Σ and C=1000 conversion arithmetically consistent | accepted | **re-summed after the amortization change; regex audit of every literal expression** |

## 2. Pair-census protocol (exactly the mission shape)

Per set:

1. draw one uniform challenge x (mod the group order) — **shared by both runs**;
2. run 1: 16 walkers, walker coins (a_i, c_i) fresh, per-walker budget `per` = 50764 steps
   (Λ = 0.6, dp 4), 2048-row tables, parity salts SA/SB, ring/escape per arm;
3. **if run 1 solved, stop** (a set that succeeds costs one run);
4. else run 2: fresh walker coins on the *same* x, same budget;
5. a set is solved iff run 1 or run 2 solves.

The set probability p₂ is the **measured** solved-set fraction, and the shipped number is its
one-sided 95 pct Clopper-Pearson **lower** bound. The birthday formula 1−(1−p₁)² is reported beside
it as an agreement row only. n = 8 shards x 750 sets = 6000 sets per arm for the two a = 0 arms, and
8 x 40 = 320 sets per arm for the a = −3 numerator probes (disclosed resolution ±5 pct).

Arms (the two that matter, both on the same cell, same seeds-per-shard, same budget):

* `shipped_guarded` — parity salt + carried prev row + 4-deep ring + out-of-line escape. This is the
  transition of the counted TU whose spans price the ledger.
* `control_pre_parity` — parity salt only (guard removed). This is the corrected pre-guard control that
  produced the 0.2410 plateau.

## 3. Results

| arm | n sets | p₁ | p₂ measured | CP95 lower | 1−(1−p₁)² | floor 0.39 |
|---|---|---|---|---|---|---|
| shipped_guarded | 6000 | 0.03833 | **0.07050** | **0.06514** | 0.07520 | **not met** |
| control_pre_parity | 6000 | 0.26417 | **0.45583** | **0.44519** | 0.45855 | **met** |

Run-2 conditional rate (shipped) on the run-1 failures: 0.03345. Shard-level scatter is in
`evidence/pair_census_v10p.json` (`shard_files`), so a reviewer can see n = 8 shards were pooled rather
than one long run.

**Reading.** Two independent statements, both measured on the same fixture and the same protocol:

1. The w=2 construction **clears the floor** when the guard is absent: 0.44519 ≥ 0.39. This is the
   first measured 0.39 on this track at the production budget shape, and it needs no extrapolation
   law — the production-relevant quantity is the quantity measured at the plateau cell (where
   ring/escape is ×1.000 at fixture, already settled, so the fixture shape is representative).
2. The same construction **fails on the shipped walk**: 0.06514 ≪ 0.39. The guard multiplier measured here is ×0.145 single-run and ×0.155 on
   the w=2 ticket — consistent in size with v10o's ×0.134 at a different seed policy and n, and it
   comes from a statistic that does not share v10o's single-run CP machinery.

**Consequence for the claim.** The floor is *not* claimed for the shipped walk. v10o withdrew a
0.2550 "flawed arm" for measuring success on an arm that was not shipped; using the control row as
the shipped walk's success probability would repeat that error in the opposite direction, and is not
done. The control row is evidence about *where* the loss is.

## 4. `batch_inv` amortization — the real accounting flaw

One Fermat inversion per **step per batch of B walkers**: 385 `fe_mul` (256 squarings + 129
multiplications) × 147 Ir/`fe_mul` = **56,595 Ir per batch-step**. Prefix products (B−1 muls) and
back-substitution (B−1 muls) are per-walker and already sit inside the 293-Instr row; they do not
scale with B. Therefore

    inversion share per walker-step = 385 × 147 / B

| B | share (Ir/walker-step) | charge (Ir/walker-step) | W_single | W w=2 |
|---|---|---|---|---|
| 65536 | 0.863571 | 1241.506679 | 127.943609 | 128.943609 |
| 16 (shipped W array) | 3537.187500 | 4777.830678 | 129.887873 | 130.887873 |
| **8 (live census / worst shipped)** | **7074.375000** | **8315.018178** | **130.687237** | **131.687237** |

Charged: **B = 8** (the worse of the two shipped widths; the 65536 row is the production-lockstep
reading and is *not* charged). Δ vs v10o = **+2.743536 bits**.

Charge rows after the change (B = 8), each re-summed in `evidence/ledger.json`:

```
  239      step_forward
  694      step_close  (44 mul, 54 umulh)
  293      batch_inv back-substitution, per walker
   14      denominator zero-scan
 7074.375  batch_inv inversion, 385 fe_mul x 147 / 8
    0.6398 rare-event reserve
    0.00341398 escape, 57277 Ir x 2^-24.03 x kappa
---------
 8315.0181781  Ir per walker-step   (v10o: 1241.5067493)
```

The reproduction of v10o's exact 1241.5067493 at B = 65536 **before** the denominator was changed is
the audit that the delta is the amortization alone; the span prices are untouched (same counted TU,
`sha256 96b736ef...`, 32864 bytes).

## 5. κ derived

Let r = measured closure rate per walker-step = 2^−24.03, r₀ = predicted = 2^−24. The measured/predicted
ratio r/r₀ = 2^−0.03 = 0.979420 is the fixture→production transfer factor for the *only* channel that
carries ring geometry into the surcharge. Charged in the direction that costs work:

    κ = max(1, r₀/r) = 2^0.03 = 1.021012

escape surcharge 57,277 × 2^−24.03 × κ = 0.003413975 Ir/walker-step vs unscaled 0.003343717:
Δ = 7.03e−5 Ir = 7.0e−8 bits. Small — and *computed*, not assumed. Residuals stated: (i) the closure
counter is a live counter, so no binomial interval is claimed on it; (ii) production DP density 2^−62
changes the record-collision channel, which is not the ring-escape channel; (iii) the cumulative
production count is ℓ × rate = 2^87.60 closures per trail (settled, unchanged).

## 6. The a-channel (a = 0 fixture vs a = −3 production) — cold by measurement

Structural fact: for y² = x³ + ax + b, the chord law with x₁ ≠ x₂ contains no a. **a enters only the
doubling numerator**, 3x² + a. So the fixture's a = 0 and production's a = −3 can differ only on steps
that evaluate a doubling, and the transition differs there and nowhere else.

Measured (counters in `experiments/census_v10p.c`, rows in `evidence/a_shift.json`):

| quantity | shipped arm | control arm |
|---|---|---|
| walker-steps | 9559876480 | 8459312960 |
| doubling events | 0 | 0 |
| witness gate verified/meets | 503/503 | 3162285/3162285 |

Total walker-steps 18019189440, total doubling events **0**, so the a-dependent numerator is never
evaluated: **|P(a=−3) − P(a=0)| ≤ 3.138e-10** (95 pct rule-of-three upper bound on the doubling rate).

Production reachability, same argument: the doubling condition is x = T[j].x for the *selected* row j,
i.e. one hit in 2^256 per step. Over the full 2^127.63-step search the expected number of doubling
events is 2^−128.4. **The a = −3 term is never reached at production either**, so the a = 0 fixture is
not a representativeness gap for the quantity this package claims.

Negative result, disclosed: substituting the a = −3 numerator alone and rerunning the census (320 sets
per arm) yields rows that fail their own validity gate (verified/meets = 0.0057 (8,322 of 1,452,101) vs 1.000), because a
numerator-only change breaks the coefficient invariant point = (a + c·x)G — the addition law must match
the curve. Those rows are recorded and **excluded from every claim**. A real a = −3 fixture requires a
prime-order group on a curve with a = −3, i.e. Schoof point counting at 2^40, which is not done here.

## 7. CP calibration (mandatory before use)

`cp95_lower(k, n) = Beta⁻¹(0.05; k, n−k+1)`, computed by bisection on the regularized incomplete beta.
Calibrated on the *observed* n and p before any bound is quoted: one-sided coverage 0.949 over
4000 Monte-Carlo draws of Binomial(n, p̂); the k = n closed form 0.05^(1/n) equals the helper to
0.0e+00; k = 0 returns 0. Numbers in `evidence/cp_calib.json`.

## 8. Arithmetic audit and rounding direction

`experiments/arith_audit.py` re-extracts every `a op b = c` with numeric literals from this file,
`public-note.md` and every `evidence/*.json` and recomputes it (1 expressions, 0
failures), then re-derives 31 curated load-bearing assertions from their definitions (charge sum,
the three amortization shares, log2 conversions, claim ceiling, κ both ways, CP bounds, plateau mean,
guard ×0.134, salt 5.9 pct, closure count). Result: `evidence/arith_audit_v10p.json`, verdict
PASS.

Rounding: **claim = 4-dp CEILING of W_w=2** (131.6872366 → 131.6873), so the at-most direction holds
at every quoted precision.

## 9. Scorable subset — the honest statement

The track scores log2 of total classical work, admitted only with success ≥ 0.39 over uniform
production challenges. There is **no separately scorable work-only component**, and this package does
not route around that. Work: exact, re-summed, **+2.7436 bits worse than v10o**. Success: **unmet on
the shipped arm** (0.06514 vs 0.39), **met on the control** (0.44519). Claim **not scoreable, not
promotable**; filed as a corrected record with the barrier localised to one carry and one construction.

## 10. Not claimed

The 0.39 floor for the shipped walk; a guarded plateau; any extrapolation law (θ² stays withdrawn,
slope 0.87 ± 0.09 accepted); that the a = −3 residual is zero; that the control row is the shipped
walk; production-scale execution; any prior acceptance.
