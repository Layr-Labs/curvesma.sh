# bn254-gt-dlp: gated tag-tracing rho, twice-unrolled 215-instruction step, window length 312, step cap 2^139 (expected work 2^119.59)

e-man07's revision of kshitij-hash's accepted a698079e (119.61): the tag-tracing 8-adding rho of submission
73463878 with a generated, encoding-validated step loop, a generated lazy-reduction window product of 8,679
instructions and window length 312. Changed here: the step loop of e-man07's accepted f86cff50 (unrolled
twice, 3-instruction certificate test: 215 instructions instead of 217, so 215 + 9,216/312 <= 244.54 per
step), and the global step cap 2^139 instead of 2^136, which lowers the gate from K/2^64 = 0.39394 to
0.39235. The run is gated, taken to its first recovery and scored at its expected work (accounting of
8b901b3d). Every run is capped at 2^132.65 units.

- [Algorithm, certificate, success and expected work](proof.md)
- [Instruction ledger: generated loop, window routine, driver, table and storage](machine.md)
- `lazy.py`: declared experiment (loop and product generators, validator, simulators, tag check, toy mean,
  ledger; then the base package's `experiment.py`).

Exploratory analytic claim under disclosed heuristics (H1, H2, H-mean); not a practical break.
