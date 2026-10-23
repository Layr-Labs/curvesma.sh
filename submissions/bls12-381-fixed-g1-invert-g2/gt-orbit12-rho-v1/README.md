# gt-orbit12-rho-v1

Candidate for `bls12-381-fixed-g1-invert-g2`: fixed-argument pairing inversion
FAPI-1(P, z) -> Q with e(P,Q) = z.

* `proof.md` — reduction, the order-12 Frobenius class quotient, the
  four-Fp-multiplication canonicalisation, the lift/kappa walk rule, the pinned
  dictionary with its disjoint bit slices and hard caps, the collision solve,
  cycle escapes, the success argument and the exact ledger. Claim: log2 W <= 126.11.
* `submission.json` — the claim, score 12611/100, and the declared experiment.
* `experiment.py` — stdlib-only. Section A verifies every field identity at
  production scale (blst tower, ord_r(p) = 12, gamma1, zeta, and the twelve-key
  formula against brute-force Frobenius). Section B runs the complete algorithm
  at toy scale through the same code path, including the exact dictionary.
  Section C evaluates the ledger in exact integer / rational arithmetic.
* `results.json` — captured output of
  `echo '{"production_checks":2,"toy":true,"trials":16,"btrials":200,"seed":1}' | python3 -I experiment.py`.

Run: `echo '{"production_checks":2,"toy":true,"trials":16,"btrials":200,"seed":1}' | python3 -I experiment.py`
