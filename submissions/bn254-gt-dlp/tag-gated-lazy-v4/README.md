# bn254-gt-dlp: gated tag-tracing rho, generated 217-instruction step, lazy window product (expected work 2^119.61)

The tag-tracing 8-adding rho of submission 73463878. Its step loop is replaced by a generated,
encoding-validated loop of 217 instructions (base 273), and its window product by a generated
lazy-reduction routine of 8,679 instructions (base 24,044), giving 217 + 9,216/312 <= 246.54 instructions per step with window length 312. The run
is gated, taken to its first recovery and scored at its expected work (accounting of 8b901b3d). Every run is capped at 2^129.66 units.

- [Algorithm, certificate, success and expected work](proof.md)
- [Instruction ledger: generated loop, window routine, driver, table and storage](machine.md)
- `lazy.py`: declared experiment (loop and product generators, validator, simulators, tag check, toy mean,
  ledger; then the base package's `experiment.py`).

Exploratory analytic claim under disclosed heuristics (H1, H2, H-mean); not a practical break.
