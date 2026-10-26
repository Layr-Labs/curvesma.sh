# Frobenius-class rho with a doubled twisted table and an honest thin cap

Claim: time_log2 = **125.71** (total work <= 2^125.71 dense-Fp12-multiplication equivalents,
C = 20,000), success probability >= 0.39, peak memory < 2^73 bytes, preprocessing <= 2^31 units
(included), no advice.

- [Algorithm, success argument and cost ledger](proof.md)
- [Executable checks and toy runs](experiment.py)

## Lineage of this revision

The machine, algebra, walk, mix, tariff, escapes, caps, heuristics and evidence are those of
the live crown `gt-tracecounted-exact-v25` (GordoAR, 125.72), a repricing of `747ad6c9` (125.74),
itself a revision of `6d4e57d1` (125.80) and of the Frobenius-class rho of `7da5f23a` (mpjunior92).
That package runs R = 1024 twisted multipliers. This revision doubles the table to R = 2048
(24,576 entries, 14.2 MB, construction <= 2^21.1 units — inside the 2^31 preprocessing budget)
and regenerates every dependent figure: the conservative branch term halves to 0.140816, eps
falls 6.45e-4 -> 2.63e-4, S 2^125.8302 -> 2^125.8300, W 2^125.7074975 -> 2^125.7067861, claim
125.72 -> **125.71**. The claim keeps a thin but positive 0.22% allowance over the ledger,
with the declared run asserting slack and tolerated drop are both positive. Section 3.4 prices
the metered abort honestly: short-walk success cost and cap-event probability as separate
quantities under the joint 0.01 H1+ budget (no exact-3/8 assertion), and part C runs the
production rule at the declared R = 2048 (measured squaring fraction 0.376, 0 class-check
failures). `experiment.py` part D reproduces the parent's R = 1024 ledger and this one through
the same code path, so the table size is the only difference.

This is an analytic generic-attack bound, not a practical break.
