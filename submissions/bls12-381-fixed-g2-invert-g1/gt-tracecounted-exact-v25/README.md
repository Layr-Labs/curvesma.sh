# Frobenius-class rho with a cyclotomic-squaring branch, priced at its counted listing

Claim: time_log2 = **125.71** (total work <= 2^125.71 dense-Fp12-multiplication equivalents,
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
125.74 -> **125.71**. The 0.025-bit cushion the parent kept above its own computed value is now
spent, so proof.md section 3.4 re-derives the abort headroom (0.17% instead of 1.7%) and states the
drift-cost table that replaces it; H1, H1+ and the experiments are otherwise unchanged.
`experiment.py` part D reproduces the parent's ledger (2^125.7147 at c = 0.9225) and this one
through the same code path, so the price vector is the only difference.

This is an analytic generic-attack bound, not a practical break.
