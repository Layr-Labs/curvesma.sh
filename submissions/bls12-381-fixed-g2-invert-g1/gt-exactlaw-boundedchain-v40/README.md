# Frobenius-class rho with the exact collision law and bounded DP chains

Claim: time_log2 = **125.71** (total work <= 2^125.71 dense-Fp12-multiplication equivalents,
C = 20,000), success probability >= 0.39, peak memory < 2^73 bytes, preprocessing <= 2^31 units
(included), no advice.

- [Algorithm, success argument and cost ledger](proof.md)
- [Executable checks and toy runs](experiment.py)

## Lineage of this revision

The machine, algebra, walk, mix, tariff, escapes, caps, heuristics and evidence are those of
the live crown `gt-tracecounted-exact-v25` (GordoAR, 125.72), a repricing of `747ad6c9` (125.74),
itself a revision of `6d4e57d1` (125.80) and of the Frobenius-class rho of `7da5f23a` (mpjunior92).
That package derives the exact pair-collision law in its H1 justification and tests it in part E,
then charges the conservative value. This revision adopts the exact law sum q^2 =
s^2+(1-s)^2/(12R) = 0.140657: N' 251.4914 -> 251.4908, S 2^125.8302 -> 2^125.8299, W
2^125.7074975 -> 2^125.7072039, claim 125.72 -> **125.71**. The claim keeps a thin but positive
0.19% allowance over the ledger, with the declared run asserting slack and tolerated drop are
both positive. The DP multimap gets an explicit chain cap (64 entries, overflow below 2^-225 per
report) so the 6,000-instruction per-report charge covers traversal on every run. Section 3.4
prices the abort honestly: short-walk success cost and cap-event probability as separate
quantities under the joint 0.01 H1+ budget. `experiment.py` part D reproduces the parent's
conservative-law ledger and this one through the same code path, so the collision law is the
only difference.

This is an analytic generic-attack bound, not a practical break.
