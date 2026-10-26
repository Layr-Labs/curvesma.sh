# Frobenius-class rho with an even squaring mix and bounded DP chains

Claim: time_log2 = **125.71** (total work <= 2^125.71 dense-Fp12-multiplication equivalents,
C = 20,000), success probability >= 0.39, peak memory < 2^73 bytes, preprocessing <= 2^31 units
(included), no advice.

- [Algorithm, success argument and cost ledger](proof.md)
- [Executable checks and toy runs](experiment.py)

## Lineage of this revision

The machine, algebra, walk, tariff, escapes, caps, heuristics and evidence design are those of
the live crown `gt-tracecounted-exact-v25` (GordoAR, 125.72), a repricing of `747ad6c9` (125.74),
itself a revision of `6d4e57d1` (125.80) and of the Frobenius-class rho of `7da5f23a` (mpjunior92).
That package squares iff T mod 8 < 3 (s = 3/8). This revision moves the threshold to T mod 8 < 4
(s = 1/2) — same 3-bit read, same bit layout — and regenerates every dependent figure: c 0.9178875
-> 0.8507, S 2^125.8302 -> 2^125.9283, W 2^125.7074975 -> 2^125.6959438, claim 125.72 ->
**125.71**, keeping a 0.98% allowance over the ledger (wider than the parent's 0.87%). The toy
first-collision evidence and Part E are re-run at s = 1/2 (750 runs, pooled mean within 0.5% of
model; Part E agrees with the exact prediction within noise), and part C runs the new threshold
in the real GT (0 class-check failures, measured fraction 0.478). The DP multimap gets an explicit
chain cap (64 entries, overflow below 2^-225 per report) so the 6,000-instruction per-report charge
covers traversal on every run, and section 3.4 prices the abort honestly: short-walk success cost
and cap-event probability as separate quantities under the joint 0.01 H1+ budget. `experiment.py`
part D reproduces the parent's s = 3/8 ledger and this one through the same code path, so the mix
is the only difference.

This is an analytic generic-attack bound, not a practical break.
