# truncated-bsgs-gt

Scored candidate for `bls12-381-fixed-g2-invert-g1`.

Given a nonidentity `Q` in G2 and `z` in GT, the package reduces the instance to
the discrete logarithm of `z` to the base `g = e(P0,Q)` in the prime-order group
GT, and solves it with a truncated baby-step giant-step grid: `m` consecutive
baby steps `g^j` and `m` giant steps `z g^(-i m)` with
`m = 144825239754852386682327417126749682654`, the least integer with
`m^2 >= (2/5) r`. The grid covers exactly the exponents `[0, m^2)`, so the attack
succeeds with probability `m^2/r >= 0.4` and its work is the constant
`2m - 2` dense Fp12 multiplications plus the record, sort and merge work
itemised in `proof.md`.

Files:

* `submission.json` - the scored claim, score `127895611/1000000`, evidence
  `proof.md`.
* `proof.md` - reduction, algorithm, correctness, exact success count, complete
  work ledger and the peak-storage bound.
* `claim.json` - the same quantities in machine-readable form plus the SHA-256
  digests of the protected rules, profile, accounting policy and cost model this
  argument was written against.

No experiment is declared: the claim is an analytical instruction-equivalent
bound, which `tracks/docs/COST_MODEL.md` states is sufficient, so nothing is
executed at intake.
