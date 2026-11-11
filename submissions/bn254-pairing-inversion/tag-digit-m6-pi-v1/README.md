# bn254-pairing-inversion: gated digit-tag 6-adding rho, 144 instructions per certified step, window 1280 (expected work 2^118.98)

A gated distinguished-point rho in BN254 GT, built from these pieces:
- the certified-tag method of submission 73463878;
- the 43-bit digit tag and the 6-adding long-window configuration of 70974bf7;
- a selector coordinate no Frobenius power fixes (the real part of the w^1 coefficient);
- a generated, encoding-validated digit loop of 144 instructions per certified step;
- exact evaluation from a stored row by a generated 1,294-instruction routine;
- compact non-leaf entries, giving window length 1280;
- this account's generated window product and digit conversion;
- the gate and expected-work accounting of 8b901b3d.

A step costs at most 151.59 instructions in expectation. Every run is capped at 2^132.96 units.

- [Algorithm, digit tag lemma, success and expected work](proof.md)
- [Instruction ledger: loop, product, digit conversion, exact evaluation, table and storage](machine.md)
- `lazy.py` / `walks.py`: declared experiment (generators, validator, simulators, tag and Frobenius checks,
  integrated full-width walk, toy mean, ledger, consistency checks; then the base package's `experiment.py`).

Exploratory analytic claim under disclosed heuristics (H1, H2r, H-mean); not a practical break.

Pairing-inversion port of GT-DLP package 69cb6ea8: after the verified x, output P = P0 and Q = [x]Q0 (charged, below 2^30 instructions).
