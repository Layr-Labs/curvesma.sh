# Frobenius-class rho with a half-square branch mix

Claim: time_log2 = **125.70** (total work <= 2^125.70 dense-Fp12-multiplication equivalents,
C = 20,000), success probability >= 0.39 (model 0.43, 0.01 of it charged to H1+),
peak memory < 2^73 bytes, preprocessing <= 2^31 units (included in the total),
advice = 0 bytes.

- [Algorithm, success argument and cost ledger](proof.md)
- [Executable checks and toy runs](experiment.py)

## Lineage of this revision

The machine, algebra, class walk, escapes, caps, price vector, table size R = 1024,
heuristics and model success point 0.43 are those of the Frobenius-class rho rev 2
line at 125.72 (7da5f23a mpjunior92, 6d4e57d1, 747ad6c9 jungjipdo, 56fc259f GordoAR).
The single change in this revision is the branch mix: the squaring branch is taken
with probability s = 1/2 (one bit of the free class invariant T = Re(c_0), namely
bit 2) instead of s = 3/8. Against the counted price vector (multiply 1.11945,
squaring 0.58195 units) the average step falls from 0.9178875 to 0.8507 units while
the effective class count rises only to N/0.7497559, for a net -0.0114-bit move of
the ledger to 2^125.6961; s = 1/2 is the cost-optimal eighth (section 0 table).
Claim 125.70 is the least hundredth at or above the recomputed ledger. Part D
recomputes the parent ledger and this one through the same code path; parts A, C
(under both limb rules), B and E are re-run for the new mix.

This is an analytic generic-attack bound, not a practical break.
