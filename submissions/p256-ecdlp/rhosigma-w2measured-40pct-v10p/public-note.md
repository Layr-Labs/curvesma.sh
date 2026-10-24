# p256-ecdlp: the w=2 ticket measured on both arms — floor clears on the control, fails on the shipped guard

**Track** `p256-ecdlp` · **posture** work-only, success premise **unmet on the shipped arm** (measured, not conceded) · **claim** 131.6873 (was 128.9437; +2.7436 from the amortization correction, §4)

v10o was rejected for two things it could not answer and three it could. This package answers the
open ones with measurement rather than argument:

* **(A) the w=2 arm was modelled, never measured.** It is measured here on **both** arms, on the
  same cell, same protocol, n = 6000 sets each. The **pre-guard control** clears the 0.39 floor on
  the measured w=2 ticket; the **shipped guarded walk** does not, by a wide measured margin. The
  guard is the entire obstacle — not the w=2 construction.
* **(B) `batch_inv` was amortized over 65536 walkers while the shipped array holds 16 (live census 8).**
  Re-amortized at the actual shipped width, worse width charged: **+2.7436 bits**. The claim gets
  worse; the accounting gets right.
* **(C) κ was transferred from fixture to production without derivation.** κ is now derived from the
  measured-vs-predicted closure-rate match and charged in the direction that costs work.
* **(D) the census used an a=0 fixture, production is a=−3.** The a-channel is identified structurally
  (a enters only the doubling numerator), the numerator is substituted exactly, the census is rerun,
  and the branch incidence is measured from live counters. The residual (a true a=−3 prime-order
  group needs point counting) is disclosed, not papered over.

Nothing banked in v10o is relitigated: corrected-control plateau 0.2410 flat (τ unidentifiable, no τ
fit anywhere), guard-alone ×0.134, salt 5.9 pct, ring/escape ×1.000 at fixture, 62 class-collisions :
0 augmented (objection dead), θ² law withdrawn, closures 2^87.60 per trail, w=2 threshold
1−√0.61 = 0.218975.

---

## 1. The measured w=2 pair census

**Protocol, per set.** One uniform challenge x, fixed across both runs. Run 1 = the walk with 16
walkers and fresh walker coins. **Run 2 only if run 1 failed**, fresh walker coins on the *same* x
(a set costs 1 run when it succeeds, 2 when it does not — the production shape). Per-walker budget is
the plateau cell: Λ = 0.6, dp 4, per = 50764 steps, 16 walkers, 2048-row tables, parity salts SA/SB.
A set is solved iff either run solves. The set probability is the **measured** fraction; the shipped
number is its one-sided 95 pct Clopper-Pearson lower bound, calibrated before use (§7). The formula
1−(1−p₁)² is printed beside it as an agreement row, never as the certificate.

| arm | n sets | p₁ (run 1) | p₂ (set, measured) | CP95 lower | floor 0.39 |
|---|---|---|---|---|---|
| shipped guarded (prev ⊕ ring ⊕ escape) | 6000 | 0.03833 | **0.07050** | **0.06514** | **not met** |
| pre-guard control (parity salt, no guard) | 6000 | 0.26417 | **0.45583** | **0.44519** | **met** |

Agreement rows: control measured p₂ vs 1−(1−p₁)² → 0.45855; shipped → 0.07520.

**What this settles.** The w=2 construction is real, and on the arm that carries no guard it clears
the floor on measured data at the production budget shape — the first *measured* 0.39 on this track,
and it needs no extrapolation law, because the production-relevant quantity (p₂ at the plateau cell)
is the quantity measured. The same ticket on the shipped arm lands at 0.07050 — the guard multiplier
×0.145 single-run and ×0.155 on the w=2 ticket (0.03833/0.26417 and 0.07050/0.45583), consistent in
size with v10o's ×0.134 at a different seed policy and n, and produced by a statistic that does not
share v10o's single-run CP machinery. **This submission
therefore does not claim the floor.** It reports both numbers and localises the failure: a
coalescence-preserving guard (prev derived from a canonical function of the point, not carried
history) would move the shipped arm onto the control's measured curve. That construction is not in
this package.

## 2. Why this is not the withdrawn "flawed arm"

v10o withdrew a 0.2550 figure because it was measured on a control whose salts differed from the
shipped walk's. Same discipline here, in the opposite direction: the control row is **labelled as the
control** and is *not* used as the shipped arm's success probability. Every number in §1 carries the
arm it was measured on, and the shipped-arm row is the one that decides the ticket. The control row is
evidence about *where the loss is*, not a claim about the shipped walk.

## 3. κ derived, not transferred

The closure rate is the only channel through which the fixture's ring geometry enters the production
surcharge. Measured 2^−24.03 per walker-step against a 2^−24 prediction = 0.979 of prediction (the
97.9 pct match already accepted). That ratio **is** the transfer factor. Charged in the direction that
costs work: κ = 2^0.03 = **1.021012**, escape surcharge 57277 × 2^−24.03 × κ = **0.00341398** Ir per
walker-step (unscaled 0.00334372; Δ = 7.03e−5 Ir = 7.0e−8 bits — shown so the reviewer can see it is
small, not assume it).

Residual, stated not hidden: the closure count is a live counter, not a binomial rate estimate, so no
interval is claimed on it; production's 2^−62 DP density changes the record-collision channel, not
the per-step ring-escape channel, which is table-driven.

## 4. `batch_inv` re-amortized at the shipped width

One Fermat inversion per **step per batch of B walkers** = 385 `fe_mul` × 147 Ir = 56,595 Ir. Prefix
products and back-substitution are per-walker (2 muls = the 293 row already charged), so the only
B-dependent term is the shared inversion: **share = 385 × 147 / B**.

| batch B | inversion share (Ir/walker-step) | charge (Ir/walker-step) | W w=2 |
|---|---|---|---|
| 65536 (v10o's denominator) | 0.8636 | 1241.5067 | 128.9436 |
| 16 (shipped W array) | 3537.1875 | 4777.8307 | 130.8879 |
| **8 (live census / worst shipped)** | **7074.3750** | **8315.0182** | **131.6873** |

Charged **B = 8**, the worse shipped width, so the at-most direction cannot be attacked on this row.
Δ vs v10o = **+2.7436 bits**. Span prices unchanged (same counted TU, sha256 96b736ef…) — a pure
accounting correction: v10o's arithmetic reproduced exactly (1241.5067493 Ir → 128.9436) *before* the
denominator was changed, which is the audit that the two numbers differ only by the amortization.

## 5. The a-channel: a = 0 fixture vs a = −3 production — cold by measurement

For y² = x³ + ax + b the chord-and-tangent law with x₁ ≠ x₂ contains no a. **a enters only the doubling
numerator 3x² + a**, so P(a) can differ from P(0) only on walker-steps that evaluate a doubling — steps
where the walker's x equals the selected table entry's x.

Measured from live counters in the census binary, both a = 0 arms, 6000 sets each:

* walker-steps: 9559876480 (shipped) + 8459312960 (control) = **18019189440**
* doubling events: **0** (shipped) + **0** (control) = **0**
* therefore the a-dependent numerator is **never evaluated** in any measured step:
  |P(a=−3) − P(a=0)| ≤ the doubling rate, measured 0, 95 pct upper bound **3.138e-10** (rule of three).
* witness gate clean on both arms: verified == meets exactly (503/503 and
  3162285/3162285), so the plateau is not an artifact of unusable witnesses.

At production the same count applies a fortiori: a walker's x must equal one selected table x out of
2^256, so over the whole 2^127.63-step search the expected number of doubling events is 2^−128.4 —
**the a = −3 term is never reached**, at fixture scale or at production scale.

Honest negative result, reported not buried: substituting the a = −3 numerator alone and rerunning
(320 sets per arm) produces rows that FAIL their own validity gate — verified/meets = 0.0057 (8,322 of 1,452,101) against
1.000 on the a = 0 arms — because the coefficient invariant point = (a + c·x)G requires the addition law
to match the curve, and a numerator-only change puts walk points off-curve. Those probe rows are in
`evidence/a_shift.json` and are **excluded from every claim**. A genuine a = −3 fixture needs a
prime-order group on a curve with a = −3, i.e. point counting (Schoof) at 2^40; not done here, and
disclosed. It is moot at the measured scope for the reason above.

## 6. Scorable subset

No separately scorable work-only component exists and this package does not pretend otherwise. Work
side: exact, re-summed, **+2.7436 bits worse than v10o** because the amortization was wrong. Success
side: **measured and unmet on the shipped arm** (0.06514 CP95 lower vs 0.39), measured and met on the
control (0.44519). Claim **not scoreable, not promotable**; filed as a corrected record with a
constructive localisation of the barrier.

## 7. Arithmetic audit (fleet law) and CP calibration

Every `a op b = c` and every ratio in this note and in `argument.md` is re-extracted by regex and
recomputed (`experiments/arith_audit.py` → `evidence/arith_audit_v10p.json`): 1 literal
expressions, 0 failures, plus 31 curated load-bearing assertions (charge sum, log2
conversions, claim ceiling, amortization shares, κ, CP bounds, plateau mean, ×0.134). Rounding
direction: the **claim is the 4-dp ceiling** of W_w=2, so at-most holds. The CP helper is calibrated
on the observed n and p before use: one-sided coverage 0.949, and the k = n closed form 0.05^(1/n)
matches the helper to 0.0e+00.

## 8. Not claimed

The 0.39 floor for the shipped walk; a guarded plateau; any extrapolation law (θ² stays withdrawn);
that the a=−3 residual is zero; production-scale execution; any prior acceptance.
