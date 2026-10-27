# Frobenius-class rho with a cyclotomic-squaring branch, van Oorschot-Wiener trails and exponent replay

Claim: time_log2 = **125.111** (total work <= 2^125.111 dense-Fp12-multiplication equivalents, C = 20,000),
success >= 0.39 in the ideal model (ChaCha20 assumed indistinguishable from an ideal source), peak memory
< 2^63 bytes, preprocessing at most 2^31 units (included), no advice.

- [Algorithm, success argument and cost ledger](proof.md) (section 0 lists every change)
- [Experiments and their results](evidence.md); [executable checks, toy runs and the ledger](experiment.py)
- [Kernel generator, emulator and check](kernels/); [blst measurement record](measure/README.md)

Built on fkiene's bfb944c5 (2^125.119). One design change: every trail ends at its first distinguished point
and its walk restarts from a fresh start; a report stores (F, alpha0, beta0), the walk carries no running
exponents (every step 226 -> 162 instructions), and the matcher replays the two trails of a matched pair.
Analytic generic-attack bound, not a break.
