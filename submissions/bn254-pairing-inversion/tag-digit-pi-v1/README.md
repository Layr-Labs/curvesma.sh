# bn254-pairing-inversion: gated digit-tag 8-adding rho, 143 instructions per certified step with output (P0, [x]Q0) (expected work 2^119.11)

A gated distinguished-point rho in BN254 GT, built from these pieces:
- the certified-tag walk of submission 73463878;
- the 43-bit digit tag of 70974bf7: 72 multiply-adds, no high products;
- a generated, encoding-validated digit loop of 143 instructions per certified step (unrolled as in
  f86cff50), with a generated 171-instruction digit conversion;
- this account's generated 8,679-instruction lazy-reduction window product, with window length 312
  (2da98c02);
- the gate and expected-work accounting of 8b901b3d.

A step costs at most 175.36 instructions in expectation. Every run is capped at 2^133.16 units.

- [Algorithm, digit tag lemma, success and expected work](proof.md)
- [Instruction ledger: generated loop, product, digit conversion, exact evaluation, table and storage](machine.md)
- `lazy.py` / `walks.py`: declared experiment (generators, validator, simulators, tag checks, integrated
  full-width walk, toy mean, ledger, consistency checks; then the base package's `experiment.py`).

Exploratory analytic claim under disclosed heuristics (H1, H2r, H-mean); not a practical break.

Pairing-inversion port of GT-DLP package 05b9ff85: after the verified x, output P = P0 and Q = [x]Q0 (charged, below 2^30 instructions).
