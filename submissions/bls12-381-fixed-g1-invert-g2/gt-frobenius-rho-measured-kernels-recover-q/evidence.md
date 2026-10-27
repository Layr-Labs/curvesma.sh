# Evidence: the experiments of `experiment.py`

Section 7 of `proof.md`; section numbers refer to `proof.md`.

## 7. Experiments (`experiment.py`, standard library only)

* **Part A (BLS12-381):** GT membership; Frobenius formula = p-th power; 4-multiplication key
  formula = keys of the 12 explicit conjugates, all distinct; kappa picks one identical canonical
  element from every conjugate; p^6 inverts; ord_r(p) = 12; 1+p^4+p^8 = 0; r-1 factorisation;
  ord_r(2); the Granger-Scott formula equals y*y with 9 Fp2 squarings and 34 Fp2 additions.
* **Part C (BLS12-381, R = 16):** a real walk of 2,000 steps per rule in GT with the branch rule
  T mod 64 and the table index from the 4 bits above the 6-bit branch field (keys only on multiply
  steps, twisted table, squaring, 2-cycle escape; no distinguished-point predicate and no R = 1024
  table). Checks: every one of the first 150 steps lands in the class of y^ M_h or y^^2 as computed
  with explicit Frobenius; the exponent invariant y = g^alpha z^beta at the end; the squaring
  fraction (37/64 expected) over the branch-rule decisions with its binomial standard deviation and
  z-score, chi-square of the low 6 bits and of the next 4 bits of Re(c_0), and the distinct T values,
  the last step with a new T and the period of T over the last 64 steps. The declared run's readings for both
  rules are quoted in proof.md section 0a.
* **Part S (bounded sampler):** exact reduction on a small modulus, the R_MAX = 256 all-reject cap returning
  the residue 1 after 256 attempts, and 20,000 seeded draws in range; prints the fixed maxima (2,320 and
  593,920 instructions; 29.696 units), acceptance r/2^255 > 1/2, exhaustion 2^-256 per sample and the
  2^-218.99 exhaustion allowance over a run under the ideal-source model (no distributional bound on the
  ChaCha20 stream is claimed).
* **Part D (ledger and measurement record):** reads `measure/results.tsv` and the six raw Callgrind
  records, asserts the three blst marginals (383, 22,996, 9,803) and each raw record's total, charges the two
  straight-line kernels shipped under `kernels/` (multiply 12,735, squaring 6,786; `kernels/check_kernels.py`
  regenerates both and checks a reduced sample against a big-integer reference), and recomputes the ledger for
  three price vectors (asserting the first two), the cap and its integer inequality, the drift table, the
  matcher's worst case (54 x 2^55 + 57 levels) and stopping allowance, the walk budget, threshold and every-run
  inequality on integers, the word widths, the report, meter and switch charges, the restart counter,
  preprocessing and the non-walk sum, and the sampler's fixed charge. The key combination's 154 instructions are
  the length of the listing `flags_check.py` executes (key_combination: 1,200 cases, every key correct mod p), and the full-width key minimum's 192 are the length of
  KEY_MIN (key_minimum: 336 cases: keys equal above limb d with a near pair at the minimum +- 2^(64e), e <= d, for every
  d, exact ties, random keys and fixed orders; index and limbs equal Python's minimum, lowest index on a tie; three
  mutants fail).
* **Part E (pair-collision law, toy fields):** sample m distinct classes, apply one step F to each and
  count coincidences of the images; N times the pair-collision rate is the quantity H1 predicts as
  1 - s^2 - (1-s)^2/(12R). Each row is one offline run (outside the runner's input schema: `trials` 0, 60,000 classes) of
  `echo '{"q":67,"trials":0,"R":64,"branch_bits":6,"sq_cut":37,"dp_bits":3,"walks":4,"seed":2603,\
  "production_checks":0,"production_steps":0,"pair_samples":60000}' | python3 -I experiment.py`
  with q, R and seed as in the row. The relative standard error is 1/sqrt(pairs).
  The declared run's part E (5,000 classes, q = 67) is a smoke test only.

  | q, R, s, seed | observed | exact prediction | conservative value used | standard error |
  | --- | ---: | ---: | ---: | ---: |
  | 67, 64, 37/64, seed 2603 | 0.6634 | 0.6655 | 0.6630 | 0.0069 |
  | 211, 16, 37/64, seed 2604 | 0.6724 | 0.6648 | 0.6546 | 0.0164 |

  At the production R = 1024 the exact prediction is 0.665757 and the
  conservative value 1 - s^2 - (1-s)^2/R is 0.665598, the latter used by section 4. Both observed rows
  are within 0.5 and 1.1 standard errors of their settings' exact and conservative values.
* **Part B (toy fields, q = 67, 139, 163, 211, 379):** the complete algorithm (keys, kappa trick, twisted
  table, squaring, 2-cycle and checkpoint escapes, DPs in a report log, section 0-ii's matcher: whole-record
  heapsort with its level count checked against the worst-case formula, adjacent dedup, pairs under the 2^16
  budget, halt at the first verified candidate, solving and g^x = z). The toy's walk phase ends a detection
  window of walks x 20 x 2^dp_bits steps after the first class collision, its stand-in for the production
  budget. It runs section 0-i's schedule at toy scale: walks round-robin, one step per turn, a 600-instruction
  switch, and one pre-debit meter against a 2^40-instruction budget (never exceeded; highest meter 1.6323e-04
  of the budget in the 0.42 rule census at the earlier, higher step prices, 1.4478e-04 in the
  census for 0.405 and at most 1.436e-04 in census400 at their own earlier step prices, 1.797e-05 in the declared run; no stops). The loop runs
  exactly `trials` trials, with no time limit. Toy keys take
  few values, so the branch/DP bits are a hash of the canonical conjugate, a stand-in for the bits of the
  381-bit T; Part C tests those bits in the real field. Part B also reports the mean-charge fluctuation
  against the independent-steps prediction (section 3.4).

Earlier census (model 0.43, superseded for the success value by the census below; all logs recovered and
verified): offline runs of the command above with `pair_samples` 0, `dp_bits` 3, `walks` 4, R = 64, branch 6
bits at 37/64, at q = 67 (seed 2601, 3,000 runs) and q = 163 (seed 2602, 900 runs), at the earlier step prices.
Mean T_c / sqrt N 1.5565 and 1.5422 against the model's 1.5392. The 2-cycle rates (proof.md section 3.3) and
the mean-charge fluctuation (section 3.4) are from these runs. The toy hash's visit-weighted branch fraction is
0.5735 and 0.5774 there (0.5716 in the declared run), so at most 0.4675 points below 37/64 in the census and
0.6575 in the declared run, the spread of short walks (H1+ bounds the whole-run average, proof.md 3.4); Part C reads the production rule at 0.5824.

Census at model 0.42 (provenance only; S(0.4001) is a shorter budget): Part B runs at the toy configuration above (`pair_samples` 0, every call kept, seeds
contiguous and disjoint) at q = 67, 163, 211 and 379. Rule, fixed before any count was read, on the rule census
alone: at every q, p - 4 sqrt(p(1-p)/n) >= 0.401, p = k/n (lower k on T_c rounding ambiguity), else 0.43. It held
at all four q. Rule census: 3,345 calls of 200 trials, seeds 5,j00,000+i, j=1..6, from a snapshot of experiment.py
with sha256 32c368062197b77a924070eb31c257ca4bea25fdbbc643d633dfeb0008d7d811 (walk and random-number functions
AST-identical to this file's; part_B differs in a time guard and printed fields), no meter stop, sealed before
23:40Z. Own stratum, not part of the rule: our earlier runs at the same configuration (2,000 trials per call,
seeds 422,j00+i, j=0..6). Last column rounded down.

| q | corpus | runs n | k (T_c <= S(0.42)) | p | SE | p - 4 SE |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 67 | rule census | 146,400 | 62,006 | 0.42354 | 0.00129 | 0.41837 |
| 67 | own (stratum) | 24,000 | 10,334 | 0.43058 | 0.00320 | 0.41779 |
| 163 | rule census | 137,600 | 58,103 | 0.42226 | 0.00133 | 0.41693 |
| 163 | own (stratum) | 12,000 | 5,093 | 0.42442 | 0.00451 | 0.40636 |
| 211 | rule census | 221,600 | 93,857 | 0.42354 | 0.00105 | 0.41934 |
| 211 | own (stratum) | 16,000 | 6,739 | 0.42119 | 0.00390 | 0.40557 |
| 379 | rule census | 163,400 | 70,028 | 0.42857 | 0.00122 | 0.42367 |
| 379 | own (stratum) | 16,000 | 6,932 | 0.43325 | 0.00392 | 0.41758 |

Census for model 0.405 (provenance only; S(0.4001) is a shorter budget):
plan and counting as for census400 below, 320,000 runs: q = 163 and 211 two jobs of 40,000 each, q = 67 and
379 one job of 80,000 each, seeds 425,000 + 100 job + batch. Rule: at every q, p - 4 SE >= 0.396 at
S(0.405) = S(0.42) sqrt(ln(1/0.595)/ln(1/0.58)); an error, a missing call or a meter stop keeps 0.42.

| q | runs | lower k | p | SE | p - 4 SE | rule p - 4 SE >= 0.396 |
|---|---|---|---|---|---|---|
| 67 | 80000 | 32687 | 0.40859 | 0.00174 | 0.40163 | pass |
| 163 | 80000 | 32600 | 0.40750 | 0.00174 | 0.40055 | pass |
| 211 | 80000 | 32552 | 0.40690 | 0.00174 | 0.39995 | pass |
| 379 | 80000 | 32725 | 0.40906 | 0.00174 | 0.40210 | pass |

Census400 for model 0.4001 (section 4), a fresh plan fixed before any count was read and
pooled with nothing: Part B of an earlier revision's experiment.py, as in the census for 0.405 (this revision changes the kernel and
toy step prices, which the toy meter and the fluctuation statistic read but the walk and T_c do not), toy
configuration above, 2,000 trials per call, 240,000 runs: q = 163 and 211 two jobs of 40,000 each, q = 67 one
job of 80,000, seeds 426,000 + 100 job + batch; q = 379 was not in this census (registered for five cores
before any data). Rule: at every tested q, p - 4 SE >= 0.391, k the runs with T_c at or below S(0.40) =
S(0.42) sqrt(ln(1/0.60)/ln(1/0.58)) on the printed 4-decimal values (lower k); an error, a missing call or a
meter stop fails that q. S(0.4001) > S(0.40), so each p is a lower bound on the fraction within S(0.4001).
p - 4 SE is rounded down.

| q | runs | lower k | p | SE | p - 4 SE | rule p - 4 SE >= 0.391 |
|---|---|---|---|---|---|---|
| 67 | 80000 | 32332 | 0.40415 | 0.00174 | 0.39719 | pass |
| 163 | 80000 | 32246 | 0.40307 | 0.00174 | 0.39611 | pass |
| 211 | 80000 | 32017 | 0.40021 | 0.00174 | 0.39325 | pass |

Toy fields at R = 64; production R = 1024 is the H1 extrapolation of section 6.

* **Part F (BLS12-381, `reduction.py`):** samples P = [b]P0 as the challenge does, computes g = e(P, Q0)
  with the pinned pairing(Q0, P) operation by operation (10,331 counted operations), checks g^r = 1 and
  g != 1, takes x uniform in [1, r) for the logarithm, computes Q = [x]Q0 by the pinned projective
  double-and-add over Fp2 and one Fp2 inversion, and checks Q on the G2 curve, [r]Q = infinity, [r]Q0 =
  infinity, e(P, Q) = z with z = g^x (the verification pairing, 10,331 operations) and e(P, infinity) = 1.
  It counts each branch of the pinned addition (generic 25, equal points 28, opposite points 9, infinity
  operand 1) and doubling (20), so the hard maximum 254 x (20 + 28) + 617 = 12,809 units does not depend on
  the sampled x.

The registration prose has two typos: it says 280,000 trials (an arithmetic slip) and gives the seed range
as 426000-426419, but its own seed rule 426000 + 100 x job + batch over the 40-batch shards spans
426000-426439. The unchanged five-shard plan and rule, as run, yield 240,000 trials. The original
registration is preserved.
