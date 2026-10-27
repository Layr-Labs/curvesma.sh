# Evidence: the experiments of `experiment.py`

This file is section 7 of `proof.md`, moved here whole to keep `proof.md` under the per-file size limit; section
numbers below refer to `proof.md`.

## 7. Experiments (`experiment.py`, standard library only)

* **Part A (BLS12-381):** GT membership; Frobenius formula = p-th power; 4-multiplication key
  formula = keys of the 12 explicit conjugates, all distinct; kappa picks one identical canonical
  element from every conjugate; p^6 inverts; ord_r(p) = 12; 1+p^4+p^8 = 0; r-1 factorisation;
  ord_r(2); the Granger-Scott formula equals y*y with 9 Fp2 squarings and 34 Fp2 additions.
* **Part C (BLS12-381, R = 16):** a real walk of 2,000 steps per rule in GT with the branch rule
  T mod 64 and the table index from the 4 bits above the 6-bit branch field (keys only on multiply
  steps, twisted table, squaring, 2-cycle escape; no distinguished-point predicate and no R = 4096
  table). Checks: every one of the first 150 steps lands in the class of y^ M_h or y^^2 as computed
  with explicit Frobenius; the exponent invariant y = g^alpha z^beta at the end; the squaring
  fraction (36/64 expected) over the branch-rule decisions with its binomial standard deviation and
  z-score, chi-square of the low 6 bits and of the next 4 bits of Re(c_0), and the distinct T values,
  the last step with a new T and the period of T over the last 64 steps. The declared run reads:
  stored rule 0.5626 over 1,998 decisions (1,124 squarings), 0.01 points above 36/64, binomial
  standard deviation 0.0111, z-score 0.01; chi-square 51.67 (63 dof) for the low 6 bits and 21.56
  (15 dof) for the next 4; 2 two-cycle escapes, 1,998 distinct T, the last new one at step 1,999, no
  period, 0 class-check failures in the first 150 steps and the exponent invariant. The integer rule
  reads 0.5558, z-score -0.6, chi-squares 50.78 and 14.25, 3 two-cycle escapes, 1,997 distinct T, the
  last new one at step 1,999, no period; its key index differs from the stored rule in 811 of the
  stored rule's 874 multiply steps.
* **Part S (bounded sampler):** exact reduction on a small modulus, the R_MAX = 256 all-reject cap returning
  the residue 1 after 256 attempts, and 20,000 seeded draws in range; prints the fixed maxima (2,320 and
  593,920 instructions; 29.696 units), acceptance r/2^255 > 1/2, exhaustion 2^-256 per sample and the
  2^-218.99 exhaustion allowance over a run under the ideal-source model (no distributional bound on the
  ChaCha20 stream is claimed).
* **Part D (ledger and measurement record):** reads `measure/results.tsv` and the six raw Callgrind
  records, asserts the three blst marginals (383, 22,996, 9,803) and each raw record's total, charges the two
  straight-line kernels shipped under `kernels/` (multiply 14,563, squaring 7,226; `kernels/check_kernels.py`
  regenerates both and checks a reduced sample against a big-integer reference), and recomputes the ledger for
  three price vectors (asserting the first two), the cap and its integer inequality, the drift table, the
  matcher's worst case (54 x 2^55 + 57 levels) and stopping allowance, the walk budget, threshold and every-run
  inequality on integers, the word widths, the report, meter and switch charges, the restart counter,
  preprocessing and the non-walk sum, and the sampler's fixed charge.
* **Part E (pair-collision law, toy fields):** sample m distinct classes, apply one step F to each and
  count coincidences of the images; N times the pair-collision rate is the quantity H1 predicts as
  1 - s^2 - (1-s)^2/(12R). Each row is one offline run (outside the runner's input schema: `trials` 0, 60,000 classes) of
  `echo '{"q":67,"trials":0,"R":64,"branch_bits":6,"sq_cut":36,"dp_bits":3,"walks":4,"seed":2603,\
  "production_checks":0,"production_steps":0,"pair_samples":60000}' | python3 -I experiment.py`
  with q, R and seed as in the row, 60,000 classes per row. The relative standard error is 1/sqrt(pairs).
  The declared run's part E (5,000 classes, q = 67) is a smoke test only.

  | q, R, s, seed | observed | exact prediction | conservative value used | standard error |
  | --- | ---: | ---: | ---: | ---: |
  | 67, 64, 36/64, seed 2603 | 0.6796 | 0.6833 | 0.6806 | 0.0103 |
  | 211, 16, 36/64, seed 2604 | 0.6908 | 0.6826 | 0.6716 | 0.0241 |

  At the production R = 4096 the exact prediction 1 - s^2 - (1-s)^2/(12R) is 0.683620 and the
  conservative value 1 - s^2 - (1-s)^2/R is 0.683577, the latter used by section 4. Both observed rows
  are within 0.4 and 0.8 standard errors of their settings' exact and conservative values.
* **Part B (toy fields, q = 67, 139, 163, 211, 379):** the complete algorithm (keys, kappa trick, twisted
  table, squaring, 2-cycle and checkpoint escapes, DPs in a report log, section 0-ii's matcher: whole-record
  heapsort with its level count checked against the worst-case formula, adjacent dedup, pairs under the 2^16
  budget, halt at the first verified candidate, solving and g^x = z). The toy's walk phase ends a detection
  window of walks x 20 x 2^dp_bits steps after the first class collision, its stand-in for the production
  budget. It runs section 0-i's schedule at toy scale: walks round-robin, one step per turn, a 600-instruction
  switch, and one pre-debit meter against a 2^40-instruction budget (never exceeded; highest meter 5.981e-05,
  1.1687e-04 and 3.986e-05 of the budget in the two censuses and the declared run; no stops). The loop runs
  exactly `trials` trials, with no time limit. Toy keys take
  few values, so the branch/DP bits are a hash of the canonical conjugate, a stand-in for the bits of the
  381-bit T; Part C tests those bits in the real field. Part B also reports the mean-charge fluctuation
  against the independent-steps prediction (section 3.4).

Census of this revision, all logs recovered and verified; `dp_bits` 3, `walks` 4; S(x) is the step budget
at model success x. The two census rows are offline runs of the command above with `pair_samples` 0 and the
row's `q`, `trials` and `seed` (about 145 s and 96 s here, beyond the runner's 20 s limit; q = 67 also exceeds its 2,000-trial limit);
rerun with the fixed-count loop of this package, they reproduce every figure in the table:

| q, N classes, R, seed | branch | runs | mean T_c / sqrt N (model) | P[T_c <= S(0.40)] | P[T_c <= S(0.42)] | P[T_c <= S(0.43)] |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| 67, 129,145, 64, seed 2601, trials 3000 | 6 bits, 36/64 | 3,000 | 1.5065 (1.5192) | 0.409 | 0.431 | 0.441 |
| 163, 4,524,905, 64, seed 2602, trials 900 | 6 bits, 36/64 | 900 | 1.537 (1.5192) | 0.387 | 0.408 | 0.416 |
| **pooled** | 6 bits, 36/64 | 3,900 | 1.513 (1.5192) | - | 0.4257 | 0.4354 |
| 67, 129,145, 64, seed 1, trials 40 (declared run) | 6 bits, 36/64 | 40 | 1.2806 (1.5192) | 0.475 | 0.475 | 0.475 |

The pooled 0.4354 at S(0.43) has binomial standard error sqrt(0.4354 x 0.5646 / 3,900) = 0.0079, 0.57 standard
errors above the model 0.43 (0.4257 pooled at S(0.42), 0.33 standard errors above the model 0.42 used here);
the q = 67 and q = 163 rows have standard errors 0.0090 and 0.0164. The
2-cycle rate at R = 64 with squaring is 3.70e-4 per step in the q = 67 census, 2.45e-4 in the q = 163
census and 1.79e-4 in the declared run, against 2.49e-4 predicted at R = 64 and the 1.302e-3 bound; the
mean-charge fluctuation rms is 0.4078 and 0.4231 against 0.4028 predicted (0.4318 in the declared run,
where 40 trials is a single small sample). The toy hash's visit-weighted branch fraction is 0.559 and
0.5625 in the two census settings (0.5554 in the declared run), each printed to 4 decimals, so at most 0.35 points below 36/64
(0.5625) in the census and at most 0.71 in the declared run,
inside the 1.5095-point H1+ allowance; Part C reads the production branch rule directly at 0.5626.

* **Part F (BLS12-381, `reduction.py`):** samples P = [b]P0 as the challenge does, computes g = e(P, Q0)
  with the pinned pairing(Q0, P) operation by operation (10,331 counted operations), checks g^r = 1 and
  g != 1, takes x uniform in [1, r) for the logarithm, computes Q = [x]Q0 by the pinned projective
  double-and-add over Fp2 and one Fp2 inversion, and checks Q on the G2 curve, [r]Q = infinity, [r]Q0 =
  infinity, e(P, Q) = z with z = g^x (the verification pairing, 10,331 operations) and e(P, infinity) = 1.
  It counts each branch of the pinned addition (generic 25, equal points 28, opposite points 9, infinity
  operand 1) and doubling (20), so the hard maximum 254 x (20 + 28) + 617 = 12,809 units does not depend on
  the sampled x.
