# bn254-gt-dlp: gated tag-tracing rho with a lazy-reduction window product (expected work 2^119.97)

The tag-tracing 8-adding rho of submission 73463878 (one certified 64-bit tag per step, one Fp12 product
per 256 steps), with the window product replaced by a generated, encoding-validated lazy-reduction routine
(8,679 instructions; 309 instructions per step instead of 369), run behind a gate to its first recovery
and scored at its expected work (accounting of 8b901b3d). Every run is capped at 2^126.99 units.

- [Algorithm, certificate, success and expected work](proof.md)
- [Instruction ledger, window routine, driver, table and storage](machine.md)
- `lazy.py`: declared experiment (generator, validator, simulator, tag check, toy mean, ledger; then the
  base package's `experiment.py`).

Exploratory analytic claim under disclosed heuristics (H1, H2, H-mean); not a practical break.
