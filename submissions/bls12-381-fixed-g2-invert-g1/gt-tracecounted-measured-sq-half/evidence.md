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
  steps, twisted table, squaring, 2-cycle escape; no distinguished-point predicate and no R = 1024
  table). Checks: every one of the first 150 steps lands in the class of y^ M_h or y^^2 as computed
  with explicit Frobenius; the exponent invariant y = g^alpha z^beta at the end; the squaring
  fraction (37/64 expected) over the branch-rule decisions with its binomial standard deviation and
  z-score, chi-square of the low 6 bits and of the next 4 bits of Re(c_0), and the distinct T values,
  the last step with a new T and the period of T over the last 64 steps. The declared run reads:
  stored rule 0.5824 over 1,997 decisions (1,163 squarings), 0.42 points above 37/64, binomial
  standard deviation 0.0111, z-score 0.38; chi-square 59.56 (63 dof) for the low 6 bits and 16.57
  (15 dof) for the next 4; 3 two-cycle escapes, 1,997 distinct T, the last new one at step 1,999, no
  period, 0 class-check failures in the first 150 steps and the exponent invariant. The integer rule
  reads 0.5877, z-score 0.86, chi-squares 57.26 and 15.43, 4 two-cycle escapes, 1,996 distinct T, the
  last new one at step 1,999, no period; its key index differs from the stored rule in 755 of the
  stored rule's 834 multiply steps.
* **Part S (bounded sampler):** exact reduction on a small modulus, the R_MAX = 256 all-reject cap returning
  the residue 1 after 256 attempts, and 20,000 seeded draws in range; prints the fixed maxima (2,320 and
  593,920 instructions; 29.696 units), acceptance r/2^255 > 1/2, exhaustion 2^-256 per sample and the
  2^-218.99 exhaustion allowance over a run under the ideal-source model (no distributional bound on the
  ChaCha20 stream is claimed).
* **Part D (ledger and measurement record):** reads `measure/results.tsv` and the six raw Callgrind
  records, asserts the three blst marginals (383, 22,996, 9,803) and each raw record's total, charges the two
  straight-line kernels shipped under `kernels/` (multiply 15,492, squaring 7,226; `kernels/check_kernels.py`
  regenerates both and checks a reduced sample against a big-integer reference), and recomputes the ledger for
  three price vectors (asserting the first two), the cap and its integer inequality, the drift table, the
  matcher's worst case (54 x 2^55 + 57 levels) and stopping allowance, the walk budget, threshold and every-run
  inequality on integers, the word widths, the report, meter and switch charges, the restart counter,
  preprocessing and the non-walk sum, and the sampler's fixed charge.
* **Part E (pair-collision law, toy fields):** sample m distinct classes, apply one step F to each and
  count coincidences of the images; N times the pair-collision rate is the quantity H1 predicts as
  1 - s^2 - (1-s)^2/(12R). This revision's runs are one per row of
  `echo '{"q":67,"trials":0,"R":64,"branch_bits":6,"sq_cut":37,"dp_bits":3,"walks":4,"seed":2603,\
  "production_checks":0,"production_steps":0,"pair_samples":60000}' | python3 -I experiment.py`
  for row 1; row 2 sets q 211, R 16 and seed 2604. Each samples 60,000 classes. The relative standard error is 1/sqrt(pairs).

  | q, R, s | observed | exact prediction | conservative value used | standard error |
  | --- | ---: | ---: | ---: | ---: |
  | 67, 64, 37/64 | 0.6634 | 0.6655 | 0.6630 | 0.0069 |
  | 211, 16, 37/64 | 0.6724 | 0.6648 | 0.6546 | 0.0164 |

  At the production R = 1024 the exact prediction 1 - s^2 - (1-s)^2/(12R) is 0.665757 and the
  conservative value 1 - s^2 - (1-s)^2/R is 0.665598, the latter used by section 4. Both observed rows
  are within 0.5 and 1.1 standard errors of their settings' exact and conservative values.
* **Part B (toy fields, q = 67, 139, 163, 211, 379):** the complete algorithm (keys, kappa trick, twisted
  table, squaring, 2-cycle and checkpoint escapes, DPs in a report log, section 0-ii's matcher: whole-record
  heapsort with its level count checked against the worst-case formula, adjacent dedup, pairs under the 2^16
  budget, halt at the first verified candidate, solving and g^x = z). The toy's walk phase ends a detection
  window of walks x 20 x 2^dp_bits steps after the first class collision, its stand-in for the production
  budget. It runs section 0-i's schedule at toy scale: walks round-robin, one step per turn, a 600-instruction
  switch, and one pre-debit meter against a 2^40-instruction budget (never exceeded; highest meter 5.024e-05,
  1.076e-04 and 2.027e-05 of the budget in the two censuses and the declared run; no stops). Toy keys take
  few values, so the branch/DP bits are a hash of the canonical conjugate, a stand-in for the bits of the
  381-bit T; Part C tests those bits in the real field. Part B also reports the mean-charge fluctuation
  against the independent-steps prediction (section 3.4).

Census of this revision, all logs recovered and verified; `dp_bits` 3, `walks` 4; S(x) is the step budget
at model success x; commands are the one above with `pair_samples` 0 and the named `q`, `trials` and `seed`. The trial count is fixed,
not time-bounded; the two census rows are offline runs above the 20-second runner limit (137 s and 90 s here):

| q, N classes, R, seed | branch | runs | mean T_c / sqrt N (model) | P[T_c <= S(0.40)] | P[T_c <= S(0.42)] | P[T_c <= S(0.43)] |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| 67, 129,145, 64, seed 2601, trials 3000 | 6 bits, 37/64 | 3,000 | 1.5565 (1.5392) | 0.402 | 0.423 | 0.4317 |
| 163, 4,524,905, 64, seed 2602, trials 900 | 6 bits, 37/64 | 900 | 1.5422 (1.5392) | 0.379 | 0.402 | 0.4133 |
| **pooled** | 6 bits, 37/64 | 3,900 | 1.5532 (1.5392) | - | - | 0.4274 |
| 67, 129,145, 64, seed 1, trials 40 (declared run) | 6 bits, 37/64 | 40 | 1.3675 (1.5392) | 0.475 | 0.500 | 0.500 |

The pooled 0.4274 has binomial standard error sqrt(0.4274 x 0.5726 / 3,900) = 0.0079, 0.33 standard
errors below the model 0.43; the q = 67 and q = 163 rows have standard errors 0.0090 and 0.0164. The
2-cycle rate at R = 64 with squaring is 3.27e-4 per step in the q = 67 census, 2.21e-4 in the q = 163
census and 2.65e-4 in the declared run, against 2.32e-4 predicted at R = 64 and the 1/(12R) = 1.302e-3 bound; the
mean-charge fluctuation rms is 0.4395 and 0.4199 against 0.4313 predicted (0.3899 in the declared run,
where 40 trials is a single small sample). The toy hash's visit-weighted branch fraction is 0.5735 and
0.5774 in the two census settings (0.5716 in the declared run), 0.47 points or less below 37/64 (0.66 in the declared run), inside
the 2.0398-point H1+ allowance; Part C reads the production branch rule (R = 16, no DP test) at 0.5824.
