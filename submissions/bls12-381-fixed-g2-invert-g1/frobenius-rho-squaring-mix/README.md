# Frobenius-class rho with a cyclotomic-squaring branch

Claim: time_log2 = **125.80** (total work <= 2^125.80 dense-Fp12-multiplication equivalents,
C = 20,000), success probability >= 0.39, peak memory < 2^73 bytes, preprocessing <= 2^31 units
(included), no advice.

- [Algorithm, success argument and cost ledger](proof.md)
- [Executable checks and toy runs](experiment.py)

Builds on the Frobenius-class rho of submission 7da5f23a (mpjunior92); adds a Granger-Scott
squaring branch, a new instruction ledger and a production-field execution of the step rule.
This is an analytic generic-attack bound, not a practical break.
