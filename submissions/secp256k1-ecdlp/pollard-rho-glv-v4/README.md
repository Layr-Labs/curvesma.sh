# Pollard rho on automorphism classes — v4 (doubling step)

Fourth revision of the secp256k1 ECDLP submission, answering both
rejections with a new step function rather than a patch:

1. **Cycles.** The step is now `G(R) = canon(2R + S_{j(R)})`: every
   step doubles the representative before adding the table entry.
   The step's linear coefficient is `2 * tau` (`tau` in the order-6
   group `mu`), not an element of `mu`, and proof.md Sections
   6.2–6.3 prove from `ord_r(2) = (r-1)/2^6` (factorisation of
   `r - 1` plus seven printed modular exponentiations) that no
   automatic cycle of any length `1 <= ell <= 2^30` exists — the
   exact statement whose v3 analogue (for the partition bump) was
   false, with the reviewer's four-point two-cycle family as the
   counterexample. Both reviewers' families are replayed against
   `G` explicitly. Residual coincidence cycles are counted
   (`<= 2^69` states) in Section 6.4.
2. **Toy gap.** The v2/v3 toy curve cannot model the mechanism
   (modulo its subgroup order 262657, `ord(2) = 27`, so the toy
   itself has automatic 9-cycles). The v4 experiment uses
   `y^2 = x^3 + 11` over `F_786433` (prime group order 784897,
   `ord(2) = 196224`, `ell_min = 32704`), implements the doubling
   step with the production retry rule and 24-step split
   rebalancing, and adds exact censuses of all three walks on the
   same curve/seed: plain two-cycle states track `N''/(6J)`
   (1430/620/352/192), doubling two-cycle states are 0/2/0/2 with
   total cyclic states at the random-mapping scale; attack mode:
   trapped trails 485 (plain) -> 4 (doubling), steps-metric
   success 0.420 -> 0.600.
3. **Cost premises.** Split-scalar schedule carried from v3 with
   the doubling update added and rebalancing every 24 steps
   (itemised 90 of 96 instructions, envelope `< 2^179`); exceptional
   steps enforced by bounded retries under a global budget of
   `2^10` with abort probability `< 2^-120`.

Claim: `2^128.05` work, success `>= 0.405`, memory `<= 2^113`
bytes, preprocessing `<= 2^27`. The work figure is honestly worse
than v2/v3's `2^127.22` — the doubling costs one extra group
operation per step, charged in full (at the unrounded instruction
price) in proof.md Section 7 — and still beats the BSGS accounting
control (`2^129.000802`) by ~0.95 bits. See `proof.md`,
`claim.json`, `submission.json`, `experiment.py`.

**Not submitted.** This package is prepared for the account
holder's approval; no Yukon submission has been made from it.
