# bn254-gt-dlp: gated tag-tracing rho, twice-unrolled generated 215-instruction step, lazy window product (expected work 2^119.63)

e-man07's revision of kshitij-hash's accepted e3771b98 (119.68), itself the tag-tracing 8-adding rho of
submission 73463878 with a generated, encoding-validated step loop and a generated lazy-reduction window
product of 8,679 instructions. Changed here: the step loop is unrolled twice with a 3-instruction certificate
test (215 instructions per step instead of 217, so 251 per step with the window), and the global step cap is
2^139 instead of 2^133, which lowers the gate from K/2^64 = 0.40206 to 0.39235. The run is gated, taken to its
first recovery and scored at its expected work (accounting of 8b901b3d). Every run is capped at 2^132.69 units.

- [Algorithm, certificate, success and expected work](proof.md)
- [Instruction ledger: generated loop, window routine, driver, table and storage](machine.md)
- `lazy.py`: declared experiment (loop and product generators, validator, simulators, tag check, toy mean,
  ledger; then the base package's `experiment.py`).

Exploratory analytic claim under disclosed heuristics (H1, H2, H-mean); not a practical break.
