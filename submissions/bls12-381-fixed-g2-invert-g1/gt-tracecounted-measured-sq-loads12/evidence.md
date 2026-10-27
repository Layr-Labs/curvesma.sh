# Evidence: the experiments of `experiment.py`

Section 7 of `proof.md`, moved here whole to keep that file under its size limit; section numbers refer to it.

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
  z-score, chi-square of the low 6 and of the next 4 bits of Re(c_0), and the distinct T values, the last
  step with a new T and the period of T over the last 64 steps. The declared run reads:
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
  2^-218.99 allowance over a run under the ideal-source model (no distributional bound on ChaCha20 is
  claimed).
* **Part D (ledger and measurement record):** reads `measure/results.tsv` and the six raw Callgrind
  records, asserts the three blst marginals (383, 22,996, 9,803) and each raw record's total, charges the two
  straight-line kernels shipped under `kernels/` (multiply 14,563, squaring 6,997; `kernels/check_kernels.py`
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
  for row 1; row 2 sets q 211, R 16 and seed 2604. Each samples 60,000 classes.

  | q, R, s | observed | exact prediction | conservative value used | standard error |
  | --- | ---: | ---: | ---: | ---: |
  | 67, 64, 37/64 | 0.6634 | 0.6655 | 0.6630 | 0.0069 |
  | 211, 16, 37/64 | 0.6724 | 0.6648 | 0.6546 | 0.0164 |

  At the production R = 1024, to six places, the exact prediction 1 - s^2 - (1-s)^2/(12R) is 0.665757 and the
  conservative value 1 - s^2 - (1-s)^2/R is 0.665598, the latter used by section 4. Both rows are within 0.5
  and 1.1 standard errors of their exact and conservative values.
* **Part B (toy fields, q = 67, 139, 163, 211, 379):** the complete algorithm (keys, kappa trick, twisted
  table, squaring, 2-cycle escape, DPs in a report log, section 0-ii's matcher with its level
  count checked against the worst-case formula, adjacent dedup, pairs under the 2^16 budget, halt at the first
  verified candidate, solving and g^x = z). The toy's walk phase ends a detection window of walks x 20 x
  2^dp_bits steps after the first class collision, its stand-in for the production budget. It runs section 0-i's schedule at toy scale: walks round-robin, one step per turn, a 600-instruction
  switch, and one pre-debit meter against a 2^40-instruction budget (never exceeded; highest meter 1.937e-05 in the
  declared run, below 1.45e-4 in the 0.405 census kept as history; no stops). Toy keys take few values, so the
  branch/DP bits are a hash of the canonical conjugate, a stand-in for the bits of the 381-bit T; Part C
  tests those bits in the real field, and part B reports the mean-charge fluctuation.

Two earlier census rows, re-run here and reproducing every column; `dp_bits` 3, `walks` 4; S(x) is the step budget
at model success x; commands are the declared one with `pair_samples` 0 and the named `q`, `trials` and `seed`,
fixed counts, offline above the runner's 20 s (137 s, 90 s):

| q, seed | runs | T_c/sqrt N (model) | S(0.40) | S(0.405) | S(0.42) | S(0.43) |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| 67, 2601, 3000 trials | 3,000 | 1.5565 (1.5392) | 0.402 | 0.406 | 0.423 | 0.4317 |
| 163, 2602, 900 trials | 900 | 1.5422 (1.5392) | 0.379 | 0.387 | 0.402 | 0.4133 |
| **pooled** | 3,900 | 1.5532 (1.5392) | - | about 0.402 | - | 0.4274 |

Census at the claim's own length S(0.405), fresh seeds, registered before the runs (`dp_bits` 3, `walks` 4, `R` 64,
`branch_bits` 6, `sq_cut` 37, `trials` 200 per seed, fixed counts, no early stop, no selection). Every registered seed
is reported; a failed or short call would count its missing trials as failures (none did: 170/170 seeds, 0 errors, no
meter stops, highest meter 1.4431e-4 of the budget). Offline runs outside the runner, 12.2-37.2 s per call:

| q, seed range | runs | T_c/sqrt N (model) | S(0.40) | S(0.405) | S(0.42) | S(0.43) | 2-cycle |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 67, 74010000-74010039 | 8,000 | 1.5347 (1.5392) | 0.4088 | 0.4131 | 0.4301 | 0.4383 | 2.95e-4 |
| 67, 74020000-74020039 | 8,000 | 1.5235 (1.5392) | 0.3977 | 0.4020 | 0.4148 | 0.4240 | 2.73e-4 |
| 163, 74030000-74030029 | 6,000 | 1.5452 (1.5392) | 0.3938 | 0.3973 | 0.4123 | 0.4240 | 2.23e-4 |
| 163, 74040000-74040029 | 6,000 | 1.5287 (1.5392) | 0.4006 | 0.4061 | 0.4216 | 0.4315 | 2.30e-4 |
| 163, 74050000-74050029 | 6,000 | 1.5281 (1.5392) | 0.4035 | 0.4085 | 0.4228 | 0.4338 | 2.32e-4 |
| **q = 67** | 16,000 | 1.5291 (1.5392) | 0.4033 | 0.4075+-0.0039 | 0.4225 | 0.4311 | 2.84e-4 |
| **q = 163** | 18,000 | 1.5340 (1.5392) | 0.3993 | 0.4040+-0.0037 | 0.4189 | 0.4297 | 2.28e-4 |
| **pooled** | 34,000 | 1.5317 (1.5392) | 0.4012 | 0.4056+-0.0027 | 0.4206 | 0.4304 | 2.40e-4 |

The pooled 0.4056 is 0.25 binomial standard errors above the model 0.405 and 5.8 above the required 0.39 (the q = 67
and q = 163 pools +0.65 and -0.28); both pools and the pool clear 0.39. The charge rms is 0.433 (q = 67) and 0.427
(q = 163) against 0.4296 predicted, the branch fraction 0.5742 and 0.5777 against 37/64 = 0.578125.

Those rows pool to 0.4274 at S(0.43) (SE 0.0079) and a 0.405 column of about 0.402 from the rounded fractions;
per-row SEs 0.0090 and 0.0164. The 2-cycle rate at R = 64 with squaring is 2.84e-4 per step in the q = 67 census,
2.28e-4 in the q = 163 census and 2.65e-4 in the declared run, against 2.32e-4 predicted and the 1/(12R) bound
about 1.302e-3. The toy hash's branch fraction is 0.5742 and 0.5777 there (0.5716 in the declared run), 0.40
points or less below 37/64 (0.66 in the declared run), inside the 2.0053-point H1+ allowance; Part C reads the
production rule (R = 16, no DP test) at 0.5824. The 0.405 census above is history; this package runs at S(0.4003), and the toy prints budgets at 0.405/0.4005/0.4003. Support is the public frontier's registered run at S(0.4005): 670 seeds at four sizes (134,000 runs), every seed reported, 54,071 pooled (0.4035, 10+ SE above 0.39), every size at or above 0.4005; S(0.4003) is 0.9997 of S(0.4005).
