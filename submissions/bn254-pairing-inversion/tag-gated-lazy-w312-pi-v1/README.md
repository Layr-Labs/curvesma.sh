# bn254-pairing-inversion: gated tag-tracing rho, generated 217-instruction step, lazy window product output (P0, [x]Q0) (expected work 2^119.60)

The tag-tracing 8-adding rho of submission 73463878. Its step loop is replaced by a generated,
encoding-validated loop of 217 instructions (base 273), and its window product by a generated
lazy-reduction routine of 8,679 instructions (base 24,044), and, with the window of 2da98c02 extended to 312 steps, 3,205/13 (below 246.54) instructions per step. The run
is gated, taken to its first recovery and scored at its expected work (accounting of 8b901b3d). Every
run is capped at 2^133.66 units.

- [Algorithm, certificate, success and expected work](proof.md)
- [Instruction ledger: generated loop, window routine, driver, table and storage](machine.md)
- `lazy.py`: declared experiment (loop and product generators, validator, simulators, tag check, toy mean,
  ledger; then the base package's `experiment.py`).

Exploratory analytic claim under disclosed heuristics (H1, H2, H-mean); not a practical break.

Pairing-inversion port of GT-DLP package cb1b64ac: after the verified x, output P = P0 and Q = [x]Q0 (charged, below 2^30 instructions).
