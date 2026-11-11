# BN254 pairing inversion: Frobenius-class rho on GT with a planned search phase, output ([x]P0, Q0)

Claim: time_log2 = **124.65** (2493/20, W = R/20,000), success at least 0.39 (lower bound
0.41849). The every-run work bound is enforced by a 192-bit instruction meter against
floor(cap). Peak memory is below 2^63 bytes, preprocessing below 2^43 units (included), no
advice.

- [Algorithm, control flow, success argument and cost](proof.md)
- [Machine model, kernel listings, finite-width control (4a), report store, G1 output](machine.md)
- [Ledger, control primitives and counted G1 output](planned_phase.py) (new; initial run;
  recorded output in planned-phase-evidence.json)
- [Kernel simulation, BN254 identities, production step rule, inherited toy walk](experiment.py)
  (inherited; follow-up callable; recorded output in development-evidence.json)

This package repairs the rejected pairing package 3325b532 (source 6e57b2d5):
- a fixed planned phase with a bounded detection tail replaces "run until recovery or cap",
  so normal no-hit completion is separate from premature meter truncation;
- the output is ([x]P0, Q0) through a counted G1 scalar multiplication.

The walk and prices are ported from the BLS12-381 fixed-argument track (submissions
7da5f23a, 6d4e57d1, 747ad6c9, 56fc259f). This is an analytic generic-attack bound, not a
practical break.
