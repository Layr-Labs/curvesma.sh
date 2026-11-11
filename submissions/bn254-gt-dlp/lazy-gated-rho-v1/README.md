# bn254-gt-dlp: gated Frobenius-class rho with lazy-reduction kernels (expected work 2^122.85)

The gated, run-to-collision Frobenius-class rho of submissions 2ece9f04 and 8b901b3d, priced with
generated straight-line scalar AArch64 Fp12 routines that use lazy (once per output coefficient)
Montgomery reduction: multiply step 9,991 instructions, squaring step 4,627, average 0.36545 units.
Fixed step caps are enforced by 192-bit step counters (no instruction meter).

- [Algorithm, success and expected-work argument](proof.md)
- [Instruction ledger and generated routines](machine.md)
- `lazy.py`: routine generator, simulator, checks and the exact ledger; `experiment.py`: base package checks.

Exploratory analytic claim under disclosed heuristics (H1, H2, H-mean); not a practical break.
