# Frobenius-class rho with a cyclotomic-squaring branch, priced at its counted listing

Claim: time_log2 = **125.72** (total work <= 2^125.72 dense-Fp12-multiplication equivalents,
C = 20,000), success probability >= 0.39, peak memory < 2^73 bytes, preprocessing <= 2^31 units
(included), no advice.

- [Algorithm, success argument and cost ledger](proof.md)
- [Executable checks and toy runs](experiment.py)

## Lineage of this revision

The machine, algebra, walk, branch, escapes, caps, heuristics and evidence are those of the live
gate package `frobenius-rho-trace-branch` (jungjipdo, submission `747ad6c9`, 125.74), a revision of
`6d4e57d1` (125.80) and of the Frobenius-class rho of `7da5f23a` (mpjunior92). That package
itemises three per-step charges and then charges each above its own counted listing: 300 for 270
instructions on every step, 2,200 for 2,119 on a multiply step, 11,400 for 11,369 on a squaring
step. This revision charges the listing itself — 270 / 2,119 / 11,369 — and regenerates every
dependent figure from it: multiply step 1.125 -> 1.11945 units, squaring step 0.585 -> 0.58195,
average step cost c 0.9225 -> 0.9178875 (= 0.995 x 0.9225), W 2^125.7147 -> 2^125.7075, claim
125.74 -> **125.72**. The parent's 0.025-bit cushion is kept as a 0.0125-bit allowance over the ledger (0.87% of the
total instead of 1.7%), so the abort still absorbs the 1.5 points of squaring-fraction scatter the
production field shows; proof.md section 3.4 re-derives that headroom and states the drift-cost table.
Section 0a also states the limb-level rule (branch bits and key minima read from the stored Montgomery
representation, which is what the priced kernels produce), and part C now runs both that rule and the
ordinary-integer rule as evidence.
`experiment.py` part D reproduces the parent's ledger (2^125.7147 at c = 0.9225) and this one
through the same code path, so the price vector is the only difference.

This is an analytic generic-attack bound, not a practical break.

Restore note: this package carries the Navigator TAKE mechanism — Frobenius-class rho on GT rev 2
(limb-level rule stated on stored limbs, priced at the counted listing, 125.72) — after the public
line through 7da5f23a (mpjunior92), 6d4e57d1, 747ad6c9 (jungjipdo) and 56fc259f (GordoAR), restaged
onto HEAD 6423c27 as gt-tracecounted-exact-v25 with no mechanism change.
